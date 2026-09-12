# -*- coding: utf-8 -*-
"""
control_q4.py — 问题四控制变量实验 (附录4物性 + 固定半径 R=2cm)
================================================================
目的: 分离"物性体系变化"与"纯几何收缩效应":
  基准A: 附录4物性 + R=2cm 固定        → t_fixed,App4
  模型B: 附录4物性 + 附件2 R(t) 收缩    → t_shrink,App4 (= 51.15 h, 已有)
  Δt_shrink = t_fixed,App4 − t_shrink,App4  = 纯几何收缩效应
  总差异   = t_Q3(57.83) − t_Q4(51.15)   = 6.68 h (问题设置变化的总差异,
                                           含物性体系变化 + 收缩两重因素)
数值: 谱配点 N=96, 0→48h 粗(1h) + 48h→72h 细(60s), 判据 C(0,t)<0.15
"""

import sys, os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
sys.path.insert(0, os.path.join(_ROOT, 'algorithms_v2'))
os.chdir(_ROOT)

import numpy as np
from scipy.integrate import solve_ivp
from spectral_core import RadialField
from env_model import load_env_models

N = 96
R = 0.02
T_FULL = 540000.0   # 150h (120h 仍未达标, 再延长)
THRESHOLD = 0.15


def main():
    T_env, C_env, _ = load_env_models()
    fld = RadialField(N, R,
        D_fun=lambda C, T: 4.2e-4 * np.exp(-0.30 / np.maximum(C, 1e-6))
                            * np.exp(-3850.0 / (T + 273.15)),
        dDdC_fun=lambda C, T: 4.2e-4 * np.exp(-0.30 / np.maximum(C, 1e-6))
                               * np.exp(-3850.0 / (T + 273.15))
                               * 0.30 / np.maximum(C, 1e-6) ** 2,
        dDdT_fun=lambda C, T: 4.2e-4 * np.exp(-0.30 / np.maximum(C, 1e-6))
                               * np.exp(-3850.0 / (T + 273.15))
                               * 3850.0 / (T + 273.15) ** 2,
        beta=8e-7, C_env_fun=C_env,
        rho_fun=lambda C: 760.0 + 90.0 * C,
        cp_fun=lambda C: 1850.0 + 2150.0 * C / (C + 1.0),
        k_fun=lambda C: 0.12 + 0.20 * C / (C + 1.0),
        dkdC_fun=lambda C: 0.20 / (C + 1.0) ** 2,
        h=25.0, T_env_fun=T_env)
    n = N - 1
    y0 = np.concatenate([np.full(n, 28.0), np.full(n, 2.55)])
    t_coarse = np.arange(3600.0, 48 * 3600.0 + 1, 3600.0)
    sol = solve_ivp(fld.coupled_rhs, (0.0, 48 * 3600.0), y0, method='BDF',
                    t_eval=t_coarse, rtol=1e-8, atol=1e-8)
    i48 = int(np.where(np.isclose(sol.t, 48 * 3600.0))[0][0])
    t_fine = np.arange(48 * 3600.0, T_FULL + 1, 60.0)
    sol2 = solve_ivp(fld.coupled_rhs, (48 * 3600.0, T_FULL), sol.y[:, i48],
                     method='BDF', t_eval=t_fine, rtol=1e-8, atol=1e-8)
    # 中心值扫描
    fld._surf_C = None
    fld._surf_T = None
    c48 = None
    for i, tt in enumerate(sol.t[:-1]):
        y = sol.y[:, i]
        yC = y[n:]
        Ts = fld._surf_T if fld._surf_T is not None else T_env(tt)
        sC = fld._surf_moisture(tt, yC, Ts)
        C_full = fld._full_field(yC, fld._axis_value(yC, sC), sC)
        T0v, sT = fld._surf_heat(tt, y[:n], C_full)
        fld._surf_C, fld._surf_T = sC, sT
        if abs(tt - 48 * 3600.0) < 1:
            c48 = float(C_full[0])
        if C_full[0] < THRESHOLD:
            print(f't_fixed,App4 = {tt:.0f} s = {tt/3600:.4f} h'); return
    for i, tt in enumerate(sol2.t):
        y = sol2.y[:, i]
        yC = y[n:]
        if abs(tt - 72*3600.0) < 1:
            print(f'中心含水率 @72h = {C_full[0]:.6f}', flush=True)
        Ts = fld._surf_T if fld._surf_T is not None else T_env(tt)
        sC = fld._surf_moisture(tt, yC, Ts)
        C_full = fld._full_field(yC, fld._axis_value(yC, sC), sC)
        T0v, sT = fld._surf_heat(tt, y[:n], C_full)
        fld._surf_C, fld._surf_T = sC, sT
        if C_full[0] < THRESHOLD:
            print(f't_fixed,App4 = {tt:.0f} s = {tt/3600:.4f} h')
            print(f'中心含水率 @48h = {c48:.6f}')
            return
    print('72h 内未达标')


if __name__ == '__main__':
    main()
