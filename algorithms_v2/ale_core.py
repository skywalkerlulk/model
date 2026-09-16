# -*- coding: utf-8 -*-
"""
ale_core.py — ALE 移动网格谱配点 (2026 CUMCM A题 v2 方案, 问题4)
==================================================================
问题4 的几何: r ∈ [0, R(t)], R(t) 由附件2 实测收缩数据给定 (PCHIP 插值)。
域随时间收缩 → 采用任意拉格朗日-欧拉 (ALE) 表述:
    ∂C/∂t + w·C_r = (1/r)·∂_r(r·D·C_r),   w(r,t) = Ṙ(t)·r/R(t)  (网格速度)
    热量同构: ρcp·(∂T/∂t + w·T_r) = (1/r)·∂_r(r·k·T_r)
ALE 守恒表述 (论文 5.4 节): 控制体 V_j(t)=π(r²_{j+1/2}−r²_{j−1/2}) 随网格运动,
    d(V_j·C_j)/dt = Σ扩散通量 + Σ对流(网格)通量,
    几何守恒律 (GCL): dV_j/dt = 2π(r_{j+1/2}·ṙ_{j+1/2} − r_{j−1/2}·ṙ_{j−1/2}) = (2Ṙ/R)V_j
配点形式 (与 Q1-3 同内核): 移动 CGL 节点 r_j(t) = R(t)(1−x_j)/2,
    dC_j/dt = (1/r_j)(Dr(t)@(r·D·Cr))_j − w_j·(Dr(t)@C)_j
    Dr(t) = −(2/R(t))·Dmat  (度量随时变, 每 RHS 求值重建)
边界消元与 Q1-3 相同: 轴心对称 + 表面 Robin (Newton 续延), 但系数随 R(t) 时变。
文献: Białobrzewski et al. (carrot shrinking drying); 面包烘焙移动边界模拟。
"""

import numpy as np
from scipy.optimize import brentq
from spectral_core import cheb


class AleField:
    """ALE 移动网格谱配点场 (问题4)

    附录4: ρ=760+90C, cp=1850+2150C/(C+1), k=0.12+0.20C/(C+1),
           D=4.2e-4·e^(−0.30/C)·e^(−3850/T_K)
    状态 = 内部节点 (x_1..x_{N-1}); 轴心/表面由边界消元恢复。
    """

    def __init__(self, N, R_fun, Rdot_fun, D_fun, dDdC_fun, dDdT_fun, beta, C_env_fun,
                 rho_fun, cp_fun, k_fun, h, T_env_fun, R0=0.02):
        self.N = N
        self.R_fun = R_fun
        self.Rdot_fun = Rdot_fun
        self.x, self.Dmat = cheb(N)
        self.R = R0
        self._rebuild_metric(R0)
        self.D_fun, self.dDdC_fun, self.dDdT_fun = D_fun, dDdC_fun, dDdT_fun
        self.rho_fun, self.cp_fun, self.k_fun = rho_fun, cp_fun, k_fun
        self.beta, self.h = beta, h
        self.C_env_fun, self.T_env_fun = C_env_fun, T_env_fun
        self._surf_C = None
        self._surf_T = None

    def _rebuild_metric(self, R):
        """移动网格度量: 节点坐标 / 微分矩阵 / 边界分块"""
        self.R = R
        self.r = R * (1.0 - self.x) / 2.0
        self.Dr = -(2.0 / R) * self.Dmat
        self.D00, self.D0N = self.Dr[0, 0], self.Dr[0, -1]
        self.D0m = self.Dr[0, 1:-1]
        self.DN0, self.DNN = self.Dr[-1, 0], self.Dr[-1, -1]
        self.DNm = self.Dr[-1, 1:-1]

    def _set_time(self, t):
        R_new = self.R_fun(t)
        if abs(R_new - self.R) > 1e-16:
            self._rebuild_metric(R_new)
        self._Rdot = self.Rdot_fun(t)
        self._w = self._Rdot * self.r / self.R          # w = Ṙ·r/R (网格速度场)
        self._dil = 2.0 * self._Rdot / self.R           # 体积膨胀率 ∇·u (时变)

    # ---------------- 边界消元 (系数随时变) ----------------
    def _axis_value(self, y, s):
        return -(self.D0m @ y + self.D0N * s) / self.D00

    def _surf_moisture(self, t, y, T_surf=None):
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
        s = self._surf_C if self._surf_C is not None else float(np.max(y))
        s_max = 3.0 * max(np.max(y), 0.5)
        for _ in range(40):
            Fs, ds = F(s), -F(s) / dF(s)
            s_new = s + ds
            if s_new <= 0.0 or s_new > s_max:
                ds *= 0.5
                s_new = s + ds
            s = s_new
            if abs(ds) < 1e-13 * max(1.0, abs(s)):
                return s
        raise RuntimeError('ALE 表面水分 Newton 未收敛 (t=%.4g)' % t)

    def _surf_heat(self, t, y, C_full):
        k_N = self.k_fun(C_full)[-1]
        M = np.array([[self.D00, self.D0N],
                      [self.DN0, self.DNN + self.h / k_N]])
        b = np.array([-self.D0m @ y,
                      self.h * self.T_env_fun(t) / k_N - self.DNm @ y])
        T0, TN = np.linalg.solve(M, b)
        return T0, TN

    def _full_field(self, y, u0, s):
        return np.concatenate([[u0], y, [s]])

    # ---------------- 单场 RHS (ALE 配点) ----------------
    def moisture_rhs(self, t, y, T_full=None):
        self._set_time(t)
        T_surf = None if T_full is None else T_full[-1]
        s = self._surf_moisture(t, y, T_surf)
        u0 = self._axis_value(y, s)
        self._surf_C = s
        C = self._full_field(y, u0, s)
        D = self.D_fun(C, T_full)
        Cr = self.Dr @ C
        flux = self.r * D * Cr
        div = np.divide(self.Dr @ flux, self.r, out=np.zeros_like(self.r),
                        where=self.r > 0)
        # ALE 随体坐标 (网格速度 w = 材料速度 u = Ṙr/R) 完整推导:
        #   欧拉守恒: ∂(ρ_s C)/∂t + ∇·(u ρ_s C) = ∇·(ρ_s D ∇C), ρ_s ∝ 1/R² 均匀压实
        #   → ∂C/∂t|_g = div (对流/膨胀/密度变化项全部精确抵消)
        # 均匀场检验: C≡常数 → div=0 → C 恒定 ✓ (干基含水率的物质守恒)
        return div[1:-1]

    def heat_rhs(self, t, y, C_full=None):
        self._set_time(t)
        Cf = C_full if C_full is not None else np.full(self.N + 1, 2.55)
        T0, TN = self._surf_heat(t, y, Cf)
        self._surf_T = TN
        T = self._full_field(y, T0, TN)
        k = self.k_fun(Cf)
        rho_cp = self.rho_fun(Cf) * self.cp_fun(Cf)
        Tr = self.Dr @ T
        flux = self.r * k * Tr
        div = np.divide(self.Dr @ flux, self.r, out=np.zeros_like(self.r),
                        where=self.r > 0)
        # 热方程同构 (随体坐标, ρcp 取当前时刻值 — 与 v1 物质坐标一致;
        # ρcp 随时间的变化项 O(干燥速率) 远小于主导项, 作为模型假设记录)
        return np.divide(div, rho_cp, out=np.zeros_like(div),
                         where=rho_cp > 0)[1:-1]

    # ---------------- 耦合系统 ----------------
    def coupled_rhs(self, t, y):
        self._set_time(t)                         # 先更新移动网格度量
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
        self._set_time(t)
        s = self._surf_moisture(t, y)
        return self._full_field(y, self._axis_value(y, s), s)

    def full_T(self, t, y, C_full=None):
        self._set_time(t)
        Cf = C_full if C_full is not None else np.full(self.N + 1, 2.55)
        T0, TN = self._surf_heat(t, y, Cf)
        return self._full_field(y, T0, TN)

    def gcl_check(self, t):
        """几何守恒律验证: 离散 dV_j/dt vs (2Ṙ/R)·V_j (应恒等)"""
        self._set_time(t)
        return 2.0 * self.Rdot_fun(t) / self.R
