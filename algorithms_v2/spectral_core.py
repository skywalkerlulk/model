# -*- coding: utf-8 -*-
"""
spectral_core.py — Chebyshev 谱配点核心 (2026 CUMCM A题 v2 独立方案)
==================================================================
空间离散: Chebyshev-Gauss-Lobatto 配点 x∈[−1,1], r = R(x+1)/2 映射到 [0,R]
    x=−1 ↔ r=0 (轴), x=+1 ↔ r=R (表面)
文献依据: Gasparin, Dutykh & Mendes (arXiv:1902.07775) — 谱方法解热湿耦合问题;
          Guo & Labrosse, Chebyshev-Spectral Method in Transport Phenomena.
边界处理: **谱边界消元 (spectral boundary elimination)** — 轴心与表面节点均不作为
    ODE 自由度, 由边界条件在每个右端求值中精确解出:
      轴心 (对称):  D_r[0]·C = 0 → C_0 = −(D0m·y + D0N·C_N)/D00    (线性)
      表面 (Robin): D(s)·C_r(s) + β(s−C_env) = 0 → 标量非线性根 (brentq)
      热量表面:     k·T_r + h(T_N−T_env) = 0 → 与轴心条件联立 2×2 线性系统
    必要性: 若轴心节点保留方程行 (2D·C_rr 形式) 而不施加对称条件, 半离散算子
    出现正实部特征值 (max ≈ +6.7e3 s⁻¹), 数值不稳定 — 消元后全部非正。
误差估计: 末位 Chebyshev 系数 |a_N| 为后验误差指示 (Gasparin 同款)。
状态约定: y = 内部节点 (x_1..x_{N-1}), 长度 N−1。
"""

import numpy as np
from scipy.fft import dct
from scipy.optimize import brentq


def cheb(N):
    """Chebyshev-Gauss-Lobatto 节点 x_j (j=0..N) 与微分矩阵 D (Trefethen, 2000)
    注意: c 含交替符号 c_j = 2δ_{j0}·(−1)^j (端点 ×2), 即
    D_{ij} = (c_i/c_j)/(x_i−x_j) = (−1)^{i+j}(…)/(x_i−x_j)"""
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.ones(N + 1)
    c[0] = c[-1] = 2.0
    c *= (-1.0) ** np.arange(N + 1)
    X = x[:, None] - x[None, :]
    X = X + np.eye(N + 1)
    D = (c[:, None] / c[None, :]) / X
    np.fill_diagonal(D, 0.0)
    D = D - np.diag(D.sum(axis=1))
    return x, D


def barycentric_eval(x_nodes, u, x_out):
    """重心插值: 配点值 u → 任意点 x_out (Berrut & Trefethen, 2004)"""
    n = len(x_nodes)
    w = np.ones(n)
    w[1::2] = -1.0
    w[0] /= 2.0; w[-1] /= 2.0
    x_out = np.atleast_1d(np.asarray(x_out, dtype=float))
    out = np.empty_like(x_out, dtype=float)
    for idx, xo in enumerate(x_out):
        dx = xo - x_nodes
        hit = np.abs(dx) < 1e-15
        if hit.any():
            out[idx] = u[hit][0]
        else:
            wj = w / dx
            out[idx] = (wj @ u) / wj.sum()
    return out


def cheb_coeffs(u):
    """配点值 → Chebyshev 系数 (DCT-I, 归一化)"""
    n = len(u) - 1
    a = dct(u, type=1) / n
    a[0] /= 2.0; a[-1] /= 2.0
    return a


def cc_weights(x):
    """Clenshaw-Curtis 积分权重: ∫_{−1}^{1} u(x) dx ≈ Σ w_j u(x_j)"""
    n = len(x) - 1
    theta = np.arccos(-x)
    w = np.zeros_like(x)
    for j in range(n + 1):
        s = 0.0
        for k in range(1, n, 2):
            s += np.cos(k * theta[j]) / (1.0 - k * k)
        w[j] = (2.0 / n) * (0.5 + s + 0.5 * np.cos(n * theta[j]) / (1.0 - n * n))
    w[0] /= 2.0; w[-1] /= 2.0
    return w


class RadialField:
    """圆柱轴对称谱配点场 (Q1-Q3; Q4 由 ALE 子类扩展)

    水分: ∂C/∂t = D(C,T)·[C_rr + C_r/r] + D'(C)·C_r²
    热量: ρcp·∂T/∂t = k(C)·[T_rr + T_r/r] + k'(C)·T_r·C_r
    状态 = 内部节点 x_1..x_{N-1} (长度 N−1); 轴心与表面由边界消元恢复。
    """

    def __init__(self, N, R, D_fun, dDdC_fun, dDdT_fun, beta, C_env_fun,
                 rho_fun, cp_fun, k_fun, dkdC_fun, h, T_env_fun):
        self.N = N
        self.R = R
        self.x, self.Dmat = cheb(N)
        # CGL 顺序: x_0=+1 … x_N=−1 → r = R(1−x)/2 使 r_0=0 (轴), r_N=R (表面)
        self.r = R * (1.0 - self.x) / 2.0
        self.Dr = -(2.0 / R) * self.Dmat           # d/dr = −(2/R)·d/dx
        self.Drr = self.Dr @ self.Dr               # d²/dr²
        # 行分块: 0=轴, 1..N-1=内部, N=表面
        self.D00, self.D0N = self.Dr[0, 0], self.Dr[0, -1]
        self.D0m = self.Dr[0, 1:-1]
        self.DN0, self.DNN = self.Dr[-1, 0], self.Dr[-1, -1]
        self.DNm = self.Dr[-1, 1:-1]
        self.D_fun, self.dDdC_fun, self.dDdT_fun = D_fun, dDdC_fun, dDdT_fun
        self.rho_fun, self.cp_fun = rho_fun, cp_fun
        self.k_fun, self.dkdC_fun = k_fun, dkdC_fun
        self.beta, self.h = beta, h
        self.C_env_fun, self.T_env_fun = C_env_fun, T_env_fun
        self._surf_C = None
        self._surf_T = None

    # ---------------- 边界消元 ----------------
    def _axis_value(self, y, s):
        """轴心对称条件 D_r[0]·C = 0 → C_0 (线性)"""
        return -(self.D0m @ y + self.D0N * s) / self.D00

    def _surf_moisture(self, t, y, T_surf=None):
        """表面 Robin 标量根: D(s)·[A + B·s] + β(s − C_env) = 0
        A, B: 表面梯度 Cr(s) 的常数项与斜率 (轴心值已按对称条件代入)。

        注意 F(s) 有两个零点: 物理支 (s ≈ 内部含水率量级, 梯度连续) 与
        "干支" (s ≈ C_env, 此时 D(s)→0 使 β 项退化出伪根)。用 Newton 连续法
        从上一时刻的表面值出发沿物理支跟踪, 避免二分支歧义。"""
        c_env = self.C_env_fun(t)
        A = self.DN0 * (-self.D0m @ y / self.D00) + self.DNm @ y
        B = self.DN0 * (-self.D0N / self.D00) + self.DNN
        Ts = None if T_surf is None else np.array([T_surf])
        def F(s):
            Ds = self.D_fun(np.array([s]), Ts)[0]
            return Ds * (A + B * s) + self.beta * (s - c_env)
        def dF(s):
            Ds = self.D_fun(np.array([s]), Ts)[0]
            dD = self.dDdC_fun(np.array([s]), Ts)[0]
            return dD * (A + B * s) + Ds * B + self.beta
        # Newton 续延: 从上一次的表面值出发沿物理支跟踪 (积分中为前一步值,
        # 重建中为前一个输出时刻的值 — 调用方按时间顺序维护 self._surf_C)
        s = self._surf_C if self._surf_C is not None else float(np.max(y))
        s_max = 3.0 * max(np.max(y), 0.5)
        for _ in range(40):
            Fs = F(s)
            ds = -Fs / dF(s)
            s_new = s + ds
            if s_new <= 0.0 or s_new > s_max:           # 步出物理区间 → 折半重试
                ds *= 0.5
                s_new = s + ds
            s = s_new
            if abs(ds) < 1e-13 * max(1.0, abs(s)):
                return s
        raise RuntimeError('表面水分边界 Newton 未收敛 (t=%.4g)' % t)

    def _surf_heat(self, t, y, C_full):
        """表面 Robin (k 冻结于 C_full) 与轴心对称条件联立 2×2:
        [D00        D0N       ] [T0]   [ −D0m·y               ]
        [DN0   DNN + h/k_N    ] [TN] = [ h·T_env/k_N − DNm·y ]"""
        k_N = self.k_fun(C_full)[-1]
        M = np.array([[self.D00, self.D0N],
                      [self.DN0, self.DNN + self.h / k_N]])
        b = np.array([-self.D0m @ y,
                      self.h * self.T_env_fun(t) / k_N - self.DNm @ y])
        T0, TN = np.linalg.solve(M, b)
        return T0, TN

    def _full_field(self, y, u0, s):
        return np.concatenate([[u0], y, [s]])

    # ---------------- 单场 RHS (乘性守恒形式) ----------------
    # 采用 (1/r)·∂_r(r·q) 的乘积链式离散 (Trefethen 极坐标谱方法标准做法):
    #   展开形式 C_rr + C_r/r 在轴附近 (r→0) 的谱离散病态 (产生正实部伪模态),
    #   乘性形式中 r·q 在轴处精确为零, 算子尺度保持正则。
    def moisture_rhs(self, t, y, T_full=None):
        """dC/dt = (1/r)·∂_r(r·D(C)·C_r), y = 内部节点 (x_1..x_{N-1})"""
        T_surf = None if T_full is None else T_full[-1]
        s = self._surf_moisture(t, y, T_surf)
        u0 = self._axis_value(y, s)
        self._surf_C = s
        C = self._full_field(y, u0, s)
        D = self.D_fun(C, T_full)
        flux = self.r * D * (self.Dr @ C)          # r·q, 轴处 = 0 精确
        f = np.divide(self.Dr @ flux, self.r, out=np.zeros_like(self.r),
                      where=self.r > 0)
        return f[1:-1]

    def heat_rhs(self, t, y, C_full=None):
        """dT/dt = (1/(ρcp·r))·∂_r(r·k(C)·T_r), y = 内部节点 (x_1..x_{N-1})"""
        Cf = C_full if C_full is not None else np.full(self.N + 1, 2.55)
        T0, TN = self._surf_heat(t, y, Cf)
        self._surf_T = TN
        T = self._full_field(y, T0, TN)
        k = self.k_fun(Cf)
        rho_cp = self.rho_fun(Cf) * self.cp_fun(Cf)
        flux = self.r * k * (self.Dr @ T)          # r·q_T, 轴处 = 0 精确
        div = np.divide(self.Dr @ flux, self.r, out=np.zeros_like(self.r),
                        where=self.r > 0)
        return np.divide(div, rho_cp, out=np.zeros_like(div), where=rho_cp > 0)[1:-1]

    # ---------------- Q2+ 耦合系统 ----------------
    def coupled_rhs(self, t, y):
        """y = [T_1..T_{N-1}; C_1..C_{N-1}] → dy/dt
        表面边界值交叉耦合 (k(C_N) 与 D(C_N,T_N)): 定点扫描 4 次"""
        n = self.N - 1
        yT, yC = y[:n], y[n:]
        sC = self._surf_moisture(t, yC, self.T_env_fun(t))
        C0 = self._axis_value(yC, sC)
        T0, sT = self._surf_heat(t, yT, self._full_field(yC, C0, sC))
        for _ in range(4):
            C_full = self._full_field(yC, C0, sC)
            T0, sT = self._surf_heat(t, yT, C_full)
            sC = self._surf_moisture(t, yC, sT)
            C0 = self._axis_value(yC, sC)
        self._surf_C, self._surf_T = sC, sT
        C_full = self._full_field(yC, C0, sC)
        T_full = self._full_field(yT, T0, sT)
        return np.concatenate([self.heat_rhs(t, yT, C_full),
                               self.moisture_rhs(t, yC, T_full)])

    # ---------------- 恢复与诊断 ----------------
    def full_C(self, t, y):
        s = self._surf_moisture(t, y)
        return self._full_field(y, self._axis_value(y, s), s)

    def full_T(self, t, y, C_full=None):
        Cf = C_full if C_full is not None else np.full(self.N + 1, 2.55)
        T0, TN = self._surf_heat(t, y, Cf)
        return self._full_field(y, T0, TN)

    def integral_radial(self, u):
        """∫₀^R u(r)·2πr dr (Clenshaw-Curtis + Jacobian)"""
        w = cc_weights(self.x)
        return np.sum(w * u * 2.0 * np.pi * self.r) * (self.R / 2.0)

    def tail_coeff(self, u):
        """末位 Chebyshev 系数 (后验误差指示)"""
        return abs(cheb_coeffs(u)[-1])
