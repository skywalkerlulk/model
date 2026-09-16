# -*- coding: utf-8 -*-
"""
q1_analytic.py — 问题1 温度场的 Duhamel 解析解 (2026 CUMCM A题 v2 独立方案)
==================================================================
结构: 圆柱 Robin 边界 (对流换热) 的特征函数展开 + Duhamel 原理。

环境分解 (env_model.py): T_env(t) = [T_∞ − ΔT·e^(−t/τ)]   (指数主项, 解析可积)
                               + δ(t)                     (残差, PCHIP 光滑)

主项闭式解 (对指数边界的 Duhamel 卷积显式求值):
    T(r,t) = T_fit(t) − ΔT·Σ_n A_n·J₀(λ_n r/R)·g_n(t)
    g_n(t) = [e^(−t/τ) − e^(−λ_n²αt/R²)] / [τλ_n²α/R² − 1]
    A_n = 2Bi / [(λ_n²+Bi²)·J₀(λ_n)],   λ_n: λJ₁(λ) = Bi·J₀(λ) 的正根
极限自检: τ→0⁺ 恢复阶跃响应; t=0 恢复初始值 (完备性 ΣA_nJ₀=1)。

残差贡献: 指数滤波递归卷积 (梯形增量, Δt=1 s):
    I_n(t+Δt) = e^(−λ_n²αΔt/R²)·I_n(t) + (Δt/2)·[δ'(t)·e^(−λ_n²αΔt/R²) + δ'(t+Δt)]
    T_res(r,t) = Σ_n A_n·J₀(λ_n r/R)·I_n(t)

文献: Carslaw & Jaeger, Conduction of Heat in Solids (2nd ed.), Ch. VII;
      Duhamel 原理 (标准热传导理论)。
"""

import numpy as np
from scipy.special import j0, j1, jn_zeros
from scipy.optimize import brentq


def robin_roots(Bi, n_roots):
    """λ·J₁(λ) = Bi·J₀(λ) 的前 n 个正根 (与 J₀ 零点交错, 逐区间二分)"""
    z = jn_zeros(0, n_roots + 2)
    roots = []
    prev = 1e-10
    for zk in z:
        x = np.linspace(prev, zk, 4001)
        f = x * j1(x) - Bi * j0(x)
        for i in range(len(x) - 1):
            if f[i] * f[i + 1] < 0:
                roots.append(brentq(lambda l: l * j1(l) - Bi * j0(l),
                                    x[i], x[i + 1], xtol=1e-14))
                break
        prev = zk
        if len(roots) >= n_roots:
            break
    return np.array(roots[:n_roots])


class DuhamelTemperature:
    """Q1 温度场解析解: T(r, t) 在任意 (r, t) 网格上的精确求值"""

    def __init__(self, R, rho, cp, k, h, T0, env_model, n_terms=80):
        self.R = R
        self.alpha = k / (rho * cp)
        self.Bi = h * R / k
        self.lam = robin_roots(self.Bi, n_terms)
        self.A = 2.0 * self.Bi / ((self.lam ** 2 + self.Bi ** 2) * j0(self.lam))
        self.T0 = T0
        self.env = env_model
        self.dT = env_model.y_inf - env_model.y0       # ΔT = T_∞ − T_0
        self.tau = env_model.tau

    def _main_term(self, r, t):
        """指数边界主项的闭式解 (t: 标量或数组)"""
        t = np.atleast_1d(np.asarray(t, dtype=float))
        fit = self.env.y_inf - self.dT * np.exp(-t / self.tau)
        out = np.zeros((len(t), len(r)))
        decay = self.lam ** 2 * self.alpha / self.R ** 2
        denom = self.tau * decay - 1.0
        for n in range(len(self.lam)):
            g = (np.exp(-t / self.tau) - np.exp(-decay[n] * t)) / denom[n]
            out += -self.dT * self.A[n] * np.outer(g, j0(self.lam[n] * r / self.R))
        return out + fit[:, None]

    def _residual_term(self, r, t, dt=1.0):
        """残差 δ(t) = env(t) − fit(t) 的响应: δ(t) − Σ A·J₀·I_n^res(t)

        注意必须包含 δ(t) 本身 (均匀项): 残差边界对内部场的贡献 =
        瞬时跟随项 δ(t) 减去扩散滞后项。t=0 时两者相消给出 T(r,0)=T_env(0)=T0。"""
        t = np.atleast_1d(np.asarray(t, dtype=float))
        decay = self.lam ** 2 * self.alpha / self.R ** 2
        main_deriv = lambda tt: self.dT / self.tau * np.exp(-tt / self.tau)
        fit = lambda tt: self.env.y_inf - self.dT * np.exp(-tt / self.tau)
        t_grid = np.arange(0.0, t[-1] + dt / 2, dt)
        dp = np.array([self.env.deriv(tt) - main_deriv(tt) for tt in t_grid])
        I = np.zeros(len(self.lam))
        out = np.zeros((len(t), len(r)))
        k = 0
        for i in range(1, len(t_grid)):
            w = np.exp(-decay * dt)
            I = w * I + 0.5 * dt * (dp[i - 1] * w + dp[i])
            while k < len(t) and abs(t_grid[i] - t[k]) < 0.5 * dt:
                delta = self.env(t_grid[i]) - fit(t_grid[i])   # δ(t) 均匀项
                out[k] = delta - np.sum(self.A[:, None]
                                        * j0(np.outer(self.lam, r / self.R))
                                        * I[:, None], axis=0)
                k += 1
        return out

    def solve(self, r, t):
        """T(r, t): r (m, 数组), t (s, 数组) → (len(t), len(r))
        主项已含初始值 (t=0 时 T_fit=T0, 级数项为零), 不可再加 T0"""
        return self._main_term(r, t) + self._residual_term(r, t)


def bessel_series_cyl(r, t, R, D, beta, u0, u_env, n_terms=80):
    """常数 D 水分扩散的 Bessel 级数解 (阶跃边界, 验证用):
    C(r,t) = u_env + (u0−u_env)·Σ A_n J₀(λ_n r/R) e^(−λ_n²Dt/R²)
    A_n = 2Bi_m/((λ_n²+Bi_m²)J₀(λ_n)), λ_n 为 λJ₁ = Bi_m·J₀ 之根"""
    Bi = beta * R / D
    lam = robin_roots(Bi, n_terms)
    A = 2.0 * Bi / ((lam ** 2 + Bi ** 2) * j0(lam))
    t = np.atleast_1d(np.asarray(t, dtype=float))
    val = np.zeros((len(t), len(r)))
    for n in range(n_terms):
        val += (u0 - u_env) * np.outer(np.exp(-lam[n] ** 2 * D * t / R ** 2),
                                       A[n] * j0(lam[n] * r / R))
    return val + u_env
