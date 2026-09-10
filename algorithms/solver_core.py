# -*- coding: utf-8 -*-
"""
solver_core.py — 统一数值内核 (2026 CUMCM A题 问题1-4 共用)
============================================================
圆柱轴对称 (一维径向) 瞬态传热传质模型的有限体积法求解器。

统一模型族 (各子问题只是物性与边界条件的切换):
    rho(C)cp(C) dT/dt = (1/r) d/dr( r·k(C) dT/dr )          [热]
    dC/dt             = (1/r) d/dr( r·D(C,T) dC/dr )         [质]
    边界 r=R:  -k dT/dr = h(T - T_env(t)),  -D dC/dr = beta(C - C_env(t))
    边界 r=0:  对称 (零通量)
    Q1: 常物性 (附录2); Q2/Q3: 变物性 (附录3); Q4: Landau 变换 (附录4, 后续扩展)

数值格式:
    空间: 顶点中心有限体积法, 均匀网格 r_j = j·dr (j = 0..N), 守恒型;
          输出点 r = 0, 0.1, 0.5 cm ... 与节点对齐 → 零插值误差
    时间: BDF2 (第一步后退欧拉起步), 固定步长, 全隐式无条件稳定
    非线性: Newton 迭代, 解析 Jacobian (面系数取调和平均, 含 dD/dC 全导数)
验证链 (validate_solver.py):
    ① 常物性阶跃边界 vs Bessel 级数解析解
    ② 网格收敛 (N 倍增, Richardson 对比)
    ③ 时间步长收敛 (dt 减半)
    ④ 水分/能量守恒残差

作者: 建模组 | 用途: 论文附录 A
"""

import numpy as np
from scipy.linalg import solve_banded

# ----------------------------------------------------------------------
# 基础工具
# ----------------------------------------------------------------------

def harmonic_mean(a_left, a_right):
    """面处物性的调和平均 H(a,b)=2ab/(a+b) (经典 FVM, 保证通量连续性)"""
    return 2.0 * a_left * a_right / (a_left + a_right)


def harmonic_mean_deriv(a_left, a_right):
    """H(a,b) 对 a, b 的偏导数"""
    denom = (a_left + a_right) ** 2
    return 2.0 * a_right ** 2 / denom, 2.0 * a_left ** 2 / denom


def thomas_solve(lo, di, up, rhs):
    """Thomas 算法求解三对角方程组 diag(di) + sub(lo) + sup(up)"""
    n = len(di)
    d = np.array(di, dtype=float)
    r = np.array(rhs, dtype=float)
    l = np.array(lo, dtype=float)
    u = np.array(up, dtype=float)
    cp = np.zeros(n - 1)
    dp = np.zeros(n)
    cp[0] = u[0] / d[0]
    dp[0] = r[0] / d[0]
    for i in range(1, n):
        denom = d[i] - l[i - 1] * cp[i - 1]
        dp[i] = (r[i] - l[i - 1] * dp[i - 1]) / denom
        if i < n - 1:
            cp[i] = u[i] / denom
    x = np.zeros(n)
    x[-1] = dp[-1]
    for i in range(n - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]
    return x


# ----------------------------------------------------------------------
# 圆柱轴对称 FVM 网格 (顶点中心)
# ----------------------------------------------------------------------

class CylinderMesh:
    """节点 r_j = j·dr (j=0..N); 控制体体积与界面面积 (按单位长度计)"""

    def __init__(self, R, N):
        self.R, self.N = R, N
        self.dr = R / N
        self.r = np.arange(N + 1) * self.dr             # 节点坐标 (m)
        rf = (np.arange(N) + 0.5) * self.dr             # 内部界面 r_{j+1/2}
        # 控制体体积 V_j = pi*(r_{j+1/2}^2 - r_{j-1/2}^2), 两端为半控制体
        self.V = np.zeros(N + 1)
        self.V[0] = np.pi * rf[0] ** 2
        self.V[1:-1] = np.pi * (rf[1:] ** 2 - rf[:-1] ** 2)
        self.V[N] = np.pi * (R ** 2 - (R - 0.5 * self.dr) ** 2)
        # 界面面积 A_{j+1/2} = 2*pi*r_{j+1/2}; 外表面面积 A_R = 2*pi*R
        self.Af = np.zeros(N + 1)
        self.Af[1:] = 2.0 * np.pi * rf
        self.AR = 2.0 * np.pi * R


# ----------------------------------------------------------------------
# 扩散-对流边值算子的统一组装
# ----------------------------------------------------------------------

def assemble_diffusion(mesh, D_node, dDdu_node, bc_coef, bc_env):
    """组装 F(u) = Σ界面扩散通量 + 边界项 及其解析 Jacobian

    离散方程:  V_j·du_j/dt = g_{j+1/2}(u_{j+1}-u_j) - g_{j-1/2}(u_j-u_{j-1})
                其中 g = D_face·A_face/dr, D_face = H(D_j, D_{j+1}) (调和平均)
    边界 j=N:   -D dU/dr = bc_coef·(u - bc_env)   ⇒  F_N 增补 -bc_coef·A_R·(u_N-bc_env)
    返回: (lo, di, up, dlo, ddi, dup, bc_vec)
          L(u) = lo·u[j-1] + di·u[j] + up·u[j+1]; J = dF/du 三对角
          (J 中含 dg/du·Δu 的全导数修正由 _apply_face_deriv 叠加)
    """
    m = mesh
    n = m.N + 1
    Df = np.zeros(n)
    Df[1:] = harmonic_mean(D_node[:-1], D_node[1:])
    g = Df * m.Af / m.dr                        # 界面传导系数 (含面积)
    # ---- 线性部分 L(u) ----
    lo = np.zeros(n)
    up = np.zeros(n)
    lo[1:] = g[1:]                              # lo[j] 乘 u[j-1], j=1..N
    up[:-1] = g[1:]                             # up[j] 乘 u[j+1], j=0..N-1
    di = -np.r_[g[1], g[1:-1] + g[2:], g[-1]]
    di[-1] -= bc_coef * m.AR                    # 边界对流项并入
    bc_vec = np.zeros(n)
    bc_vec[-1] = bc_coef * m.AR * bc_env
    dlo, ddi, dup = lo.copy(), di.copy(), up.copy()
    return lo, di, up, dlo, ddi, dup, bc_vec


def _apply_face_deriv(dlo, ddi, dup, mesh, D_node, dDdu_node, u):
    """将 dg/du·(Δu) 的贡献叠加进 Jacobian 三对角 (F 对 u 的全导数修正)

    界面 j-1/2 (j=1..N) 连接节点 j-1, j; g 随两端物性变化:
        dg_l = dH/da·D'(u_{j-1})·A/dr, dg_r = dH/db·D'(u_j)·A/dr
    F_j 含 +g(u_j-u_{j-1}) → dF_j/du_{j-1} += dg_l·Δu, dF_j/du_j += dg_r·Δu
    F_{j-1} 含 -g(u_j-u_{j-1}) → dF_{j-1}/du_{j-1} -= dg_l·Δu, dF_{j-1}/du_j -= dg_r·Δu
    """
    m = mesh
    j = np.arange(1, m.N + 1)
    Ha, Hb = harmonic_mean_deriv(D_node[j - 1], D_node[j])
    dg_l = Ha * dDdu_node[j - 1] * m.Af[j] / m.dr
    dg_r = Hb * dDdu_node[j] * m.Af[j] / m.dr
    du_face = u[j] - u[j - 1]
    ddi[j] += dg_r * du_face
    ddi[j - 1] += -dg_l * du_face
    dlo[j] += dg_l * du_face
    dup[j - 1] += -dg_r * du_face
    return dlo, ddi, dup


# ----------------------------------------------------------------------
# 变物性耦合传热传质求解器 (Q1-Q3 用; Q4 由 Landau 变换扩展)
# ----------------------------------------------------------------------

class CylinderDryer:
    """FVM + BDF2 + Newton 求解器

    参数:
        R: 药材半径 (m); N: 网格数
        rho_fun(C), cp_fun(C), k_fun(C): 热物性
        D_fun(C,T), dDdC_fun(C,T), dDdT_fun(C,T): 扩散系数及偏导
        h: 对流换热系数 (W/(m2·K)); beta: 对流传质系数 (m/s)
        T_env_fun(t), C_env_fun(t): 烘房环境
        T0, C0: 初始条件
    """

    def __init__(self, R, N, rho_fun, cp_fun, k_fun,
                 D_fun, dDdC_fun, dDdT_fun,
                 h, beta, T_env_fun, C_env_fun, T0, C0):
        self.mesh = CylinderMesh(R, N)
        self.rho_fun, self.cp_fun, self.k_fun = rho_fun, cp_fun, k_fun
        self.D_fun, self.dDdC_fun, self.dDdT_fun = D_fun, dDdC_fun, dDdT_fun
        self.h, self.beta = h, beta
        self.T_env_fun, self.C_env_fun = T_env_fun, C_env_fun
        self.T = np.full(N + 1, T0, dtype=float)
        self.C = np.full(N + 1, C0, dtype=float)
        self.t = 0.0
        self._T_prev, self._C_prev = None, None
        self._T_prev2, self._C_prev2 = None, None

    # ---------------- 单方程 BDF2 + Newton 推进 ----------------
    def _advance_field(self, u_init, u_hist_now, u_hist_prev2, dt, t_new, mass, op_fun):
        """mass·du/dt = F(u,t); BDF2 (首步/重启时退化为 BE); 返回 u^{n+1}

        u_init:      Newton 迭代初值 (GS 迭代中为上一轮扫描值)
        u_hist_now:  BDF2 历史项 u^n (本步入口的收敛状态, 迭代期间不变)
        u_hist_prev2: BDF2 历史项 u^{n-1} (None → 后退欧拉)
        """
        u = u_init.copy()
        if u_hist_prev2 is None:                 # 第一步或步长切换: 后退欧拉
            coef = mass / dt
            hist = mass * u_hist_now / dt
        else:
            coef = 3.0 * mass / (2.0 * dt)
            hist = mass * (4.0 * u_hist_now - u_hist_prev2) / (2.0 * dt)
        for _ in range(40):
            lo, di, up, dlo, ddi, dup, bc = op_fun(u, t_new)
            F = np.zeros_like(u)
            F[:-1] += up[:-1] * u[1:]
            F[1:] += lo[1:] * u[:-1]
            F += di * u
            F += bc
            G = coef * u - hist - F
            J_lo = -dlo
            J_di = coef - ddi
            J_up = -dup
            ab = np.zeros((3, len(u)))
            ab[0, 1:] = J_up[:-1]
            ab[1, :] = J_di
            ab[2, :-1] = J_lo[1:]
            delta = solve_banded((1, 1), ab, -G)
            u = u + delta
            tol = 1e-12 * max(1.0, np.max(np.abs(u)))   # 相对容差 (场值量级)
            if np.max(np.abs(delta)) < tol:
                return u
        raise RuntimeError("Newton 未收敛 (dt=%.4g)" % dt)

    # ---------------- 热方程算子: F_T, 质量阵 ρcp·V ----------------
    def _heat_op(self, T, t_new, C=None):
        m = self.mesh
        Cc = self.C if C is None else C
        k = self.k_fun(Cc)
        lo, di, up, dlo, ddi, dup, bc = assemble_diffusion(
            m, k, np.zeros_like(k), self.h, self.T_env_fun(t_new))
        mass = self.rho_fun(Cc) * self.cp_fun(Cc) * m.V
        return lo, di, up, dlo, ddi, dup, bc, mass

    # ---------------- 水分方程算子: F_C, 质量阵 V ----------------
    def _moisture_op(self, C, t_new, T=None):
        m = self.mesh
        Tc = self.T if T is None else T
        D = self.D_fun(C, Tc)
        dDdC = self.dDdC_fun(C, Tc)
        lo, di, up, dlo, ddi, dup, bc = assemble_diffusion(
            m, D, dDdC, self.beta, self.C_env_fun(t_new))
        dlo, ddi, dup = _apply_face_deriv(dlo, ddi, dup, m, D, dDdC, C)
        mass = m.V
        return lo, di, up, dlo, ddi, dup, bc, mass

    # ---------------- 时间推进 ----------------
    def step(self, dt, force_restart=False):
        """推进一个时间步, 块 Gauss-Seidel 交替迭代 (Q1 时 1 次扫描即收敛)"""
        t_new = self.t + dt
        T_now, C_now = self.T.copy(), self.C.copy()      # u^n: BDF2 历史基准
        T_prev2, C_prev2 = self._T_prev2, self._C_prev2  # u^{n-1}
        if force_restart:                        # 步长切换: 重置为后退欧拉
            T_prev2 = C_prev2 = None
        T_it, C_it = T_now, C_now
        for sweep in range(12):
            # 1) 热场: k(C), ρcp(C) 在 C_it 下冻结 → 线性三对角求解
            op = lambda T, tt: self._heat_op(T, tt, C=C_it)[:7]
            mass_T = self._heat_op(T_it, t_new, C=C_it)[7]
            T_new = self._advance_field(T_it, T_now, T_prev2, dt, t_new, mass_T, op)
            # 2) 水分场: D(C,T_new) 中 T 冻结 → Newton 迭代
            op = lambda C, tt: self._moisture_op(C, tt, T=T_new)[:7]
            mass_C = self._moisture_op(C_it, t_new, T=T_new)[7]
            C_new = self._advance_field(C_it, C_now, C_prev2, dt, t_new, mass_C, op)
            dT = np.max(np.abs(T_new - T_it)); dC = np.max(np.abs(C_new - C_it))
            T_it, C_it = T_new, C_new
            tol_sweep = 1e-12 * max(1.0, np.max(np.abs(T_it)), np.max(np.abs(C_it)))
            if max(dT, dC) < tol_sweep:
                break
        self.T, self.C = T_it, C_it
        # 历史状态滚动
        self._T_prev2, self._C_prev2 = self._T_prev, self._C_prev
        self._T_prev, self._C_prev = self.T.copy(), self.C.copy()
        self.t = t_new

    def run(self, t_end, dt, t_out=None):
        """推进到 t_end; t_out (list) 中时刻的状态随推进保存"""
        snapshots = []
        out_idx = 0
        while self.t < t_end - 1e-12:
            dt_eff = min(dt, t_end - self.t)
            self.step(dt_eff)
            if t_out is not None and out_idx < len(t_out) and self.t >= t_out[out_idx] - 1e-12:
                snapshots.append((self.t, self.T.copy(), self.C.copy()))
                out_idx += 1
        return snapshots

    # ---------------- 守恒校验 ----------------
    def moisture_conservation(self, t0, t1, C_hist):
        """d/dt ∫C dV = -A_R·β(C_R - C_env) 两端的相对偏差"""
        m = self.mesh
        C0, C1 = C_hist
        lhs = np.sum(m.V * (C1 - C0)) / (t1 - t0)
        rhs = -m.AR * self.beta * ((C0[-1] + C1[-1]) / 2.0 - self.C_env_fun((t0 + t1) / 2.0))
        return abs(lhs - rhs) / max(abs(lhs), abs(rhs), 1e-20)


# ----------------------------------------------------------------------
# 问题4: 收缩干燥 (移动边界, ξ = r/R(t) 仿射收缩物质坐标)
# ----------------------------------------------------------------------
#
# 推导要点 (论文 5.4.1 节): 在仿射收缩假设下 ξ 为物质坐标, 干基含水率 C 依附于
# 干物质, 控制体 (固定 ξ 区间) 内干物质质量守恒, 水分仅通过扩散相对干物质迁移:
#     d/dt(V_j(t)·C_j) = Σ g·(ΔC) ,   V_j(t) = πR²(t)(ξ²_{j+1/2}-ξ²_{j-1/2})
#     g_{j} = D_face·2πξ_{j-1/2}/dξ   (R 在扩散通量中恰好抵消)
#     边界 ξ=1: -D/R·∂C/∂ξ = β(C - C_env) → 通量 β·2πR(t)·(C - C_env)
# 等价 PDE:  ∂C/∂t = 1/(R²ξ) ∂/∂ξ( ξ·D·∂C/∂ξ )
# 热方程同构: 热质量 ρ(C)cp(C)·V(t), 边界 h·2πR(t)·(T - T_env)

class XiMesh:
    """ξ 网格 (0..1), 界面"面积" Af = 2πξ_{j-1/2} (R 已并入时变体积/边界)"""

    def __init__(self, N):
        self.N = N
        self.dr = 1.0 / N                       # dξ (复用 CylinderMesh 接口名)
        xi_f = (np.arange(N) + 0.5) * self.dr
        self.Af = np.zeros(N + 1)
        self.Af[1:] = 2.0 * np.pi * xi_f        # 无量纲界面周长
        self.AR = 2.0 * np.pi                    # 由 step() 更新为 2πR(t)


class ShrinkDryer:
    """收缩干燥求解器: ξ 物质坐标 + FVM + BDF2 + 块 GS 迭代

    方程 (物质坐标, 仿射收缩):  V_j(t)·dC_j/dt = Σ D_face·2πξ/dξ·ΔC
      V_j(t) = πR²(t)(ξ²_{j+1/2}-ξ²_{j-1/2}) 仅为时变系数 (不在时间导数内);
      等价 PDE: ∂C/∂t = 1/(R²ξ)·∂/∂ξ(ξD∂C/∂ξ), 边界 -D/R·∂C/∂ξ = β(C-C_env)
    热方程同构: ρcp·V(t)·dT/dt = Σ k_face·2πξ/dξ·ΔT + h·2πR(t)·(T_env-T)
    """

    def __init__(self, N, rho_fun, cp_fun, k_fun, D_fun, dDdC_fun, dDdT_fun,
                 h, beta, T_env_fun, C_env_fun, R_fun, T0, C0, R0=0.02):
        self.mesh = XiMesh(N)
        self.N = N
        self.rho_fun, self.cp_fun, self.k_fun = rho_fun, cp_fun, k_fun
        self.D_fun, self.dDdC_fun, self.dDdT_fun = D_fun, dDdC_fun, dDdT_fun
        self.h, self.beta = h, beta
        self.T_env_fun, self.C_env_fun = T_env_fun, C_env_fun
        self.R_fun = R_fun
        self.T = np.full(N + 1, T0, dtype=float)
        self.C = np.full(N + 1, C0, dtype=float)
        self.t = 0.0
        self.R = R0
        self._xi2 = self._xi2_factors()        # ξ²_{j+1/2}-ξ²_{j-1/2} 几何因子
        self._T_prev, self._C_prev = None, None
        self._T_prev2, self._C_prev2 = None, None

    def _xi2_factors(self):
        m = self.mesh
        xi_f = (np.arange(m.N) + 0.5) * m.dr
        vol = np.zeros(m.N + 1)
        vol[0] = xi_f[0] ** 2
        vol[1:-1] = xi_f[1:] ** 2 - xi_f[:-1] ** 2
        vol[-1] = 1.0 - xi_f[-1] ** 2
        return vol

    def _V_of_R(self, R):
        return np.pi * R ** 2 * self._xi2

    # ---------------- 算子 (ξ 坐标) ----------------
    def _heat_op(self, T, t_new, C=None):
        m = self.mesh
        Cc = self.C if C is None else C
        k = self.k_fun(Cc)
        lo, di, up, dlo, ddi, dup, bc = assemble_diffusion(
            m, k, np.zeros_like(k), self.h, self.T_env_fun(t_new))
        mass = self.rho_fun(Cc) * self.cp_fun(Cc) * self._V_of_R(self.R)
        return lo, di, up, dlo, ddi, dup, bc, mass

    def _moisture_op(self, C, t_new, T=None):
        m = self.mesh
        Tc = self.T if T is None else T
        D = self.D_fun(C, Tc)
        dDdC = self.dDdC_fun(C, Tc)
        lo, di, up, dlo, ddi, dup, bc = assemble_diffusion(
            m, D, dDdC, self.beta, self.C_env_fun(t_new))
        dlo, ddi, dup = _apply_face_deriv(dlo, ddi, dup, m, D, dDdC, C)
        mass = self._V_of_R(self.R)
        return lo, di, up, dlo, ddi, dup, bc, mass

    # ---------------- 单方程 BDF2 + Newton (复用 CylinderDryer 同构逻辑) ----------------
    def _advance_field(self, u_init, u_hist_now, u_hist_prev2, dt, t_new, mass, op_fun):
        u = u_init.copy()
        if u_hist_prev2 is None:
            coef = mass / dt
            hist = mass * u_hist_now / dt
        else:
            coef = 3.0 * mass / (2.0 * dt)
            hist = mass * (4.0 * u_hist_now - u_hist_prev2) / (2.0 * dt)
        for _ in range(40):
            lo, di, up, dlo, ddi, dup, bc = op_fun(u, t_new)
            F = np.zeros_like(u)
            F[:-1] += up[:-1] * u[1:]
            F[1:] += lo[1:] * u[:-1]
            F += di * u
            F += bc
            G = coef * u - hist - F
            ab = np.zeros((3, len(u)))
            ab[0, 1:] = -dup[:-1]
            ab[1, :] = coef - ddi
            ab[2, :-1] = -dlo[1:]
            delta = solve_banded((1, 1), ab, -G)
            u = u + delta
            tol = 1e-12 * max(1.0, np.max(np.abs(u)))   # 相对容差 (场值量级)
            if np.max(np.abs(delta)) < tol:
                return u
        raise RuntimeError("Newton 未收敛 (dt=%.4g)" % dt)

    # ---------------- 时间推进 ----------------
    def step(self, dt, force_restart=False):
        t_new = self.t + dt
        R_new = self.R_fun(t_new)
        self.R = R_new                            # 算子按 t^{n+1} 的体积/面积
        self.mesh.AR = 2.0 * np.pi * R_new        # 时变表面积
        C_now, T_now = self.C.copy(), self.T.copy()
        T_prev2, C_prev2 = self._T_prev2, self._C_prev2
        if force_restart:                        # 步长切换: 重置为后退欧拉
            T_prev2 = C_prev2 = None
        T_it, C_it = T_now, C_now
        for sweep in range(12):
            op = lambda T, tt: self._heat_op(T, tt, C=C_it)[:7]
            mass_T = self._heat_op(T_it, t_new, C=C_it)[7]
            T_new = self._advance_field(T_it, T_now, T_prev2, dt, t_new, mass_T, op)
            op = lambda C, tt: self._moisture_op(C, tt, T=T_new)[:7]
            mass_C = self._moisture_op(C_it, t_new, T=T_new)[7]
            C_new = self._advance_field(C_it, C_now, C_prev2, dt, t_new, mass_C, op)
            dT = np.max(np.abs(T_new - T_it)); dC = np.max(np.abs(C_new - C_it))
            T_it, C_it = T_new, C_new
            tol_sweep = 1e-12 * max(1.0, np.max(np.abs(T_it)), np.max(np.abs(C_it)))
            if max(dT, dC) < tol_sweep:
                break
        self.T, self.C = T_it, C_it
        self._T_prev2, self._C_prev2 = self._T_prev, self._C_prev
        self._T_prev, self._C_prev = self.T.copy(), self.C.copy()
        self.t = t_new
