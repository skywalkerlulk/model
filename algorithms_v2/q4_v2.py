# -*- coding: utf-8 -*-
"""
q4_v2.py — 问题4 独立求解 (2026 CUMCM A题 v2 方案, ALE 移动网格谱配点)
============================================================
模型: ALE 移动网格 + 附录4 物性 + 附件2 实测收缩
    ρ=760+90C, cp=1850+2150C/(C+1), k=0.12+0.20C/(C+1),
    D=4.2e-4·e^(−0.30/C)·e^(−3850/T_K)
网格: r_j(t) = R(t)(1−x_j)/2, R(t) 附件2 PCHIP (尾部平台外推), Ṙ 解析导数
方程: ∂C/∂t + w·C_r = (1/r)∂_r(rD·C_r), w = Ṙr/R (ALE 配点)
输出: 表6 (每6h, 0..药材表面), result4.xlsx (60s × 0.1cm + 表面列, 至烘干结束)
验证: GCL 一致性 / 与 v1 物质坐标 Q4 交叉对比 (V6)
"""

import sys, os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
sys.path.insert(0, os.path.join(_ROOT, 'algorithms_v2'))
os.chdir(_ROOT)

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy.integrate import solve_ivp
from openpyxl import Workbook

from env_model import load_env_models
from spectral_core import barycentric_eval
from ale_core import AleField

N = 96
T_FULL = 259200.0
T_OUT60 = np.arange(60.0, T_FULL + 0.1, 60.0)
OUT_R_CM = np.arange(0.0, 1.9, 0.1)     # result4 固定列 0.0..1.9, 末列=药材表面
T0, C0 = 28.0, 2.55
THRESHOLD = 0.15


def make_R():
    df = pd.read_excel('A题/附件/附件2.xlsx')
    t = df['时间'].values.astype(float)
    R_cm = df['半径'].values.astype(float)
    p = PchipInterpolator(t, R_cm / 100.0)
    pd_ = p.derivative()
    R_tail = float(R_cm[-1]) / 100.0
    def R_fun(tt):
        return float(p(tt)) if tt <= t[-1] else R_tail
    def Rd_fun(tt):
        return float(pd_(tt)) if tt < t[-1] else 0.0
    return R_fun, Rd_fun, t, R_cm


def make_ale_field(T_env, C_env, R_fun, Rd_fun):
    return AleField(N, R_fun, Rd_fun,
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
        h=25.0, T_env_fun=T_env, R0=0.02)


def main():
    T_env, C_env, _ = load_env_models()
    R_fun, Rd_fun, t_R, R_cm = make_R()
    fld = make_ale_field(T_env, C_env, R_fun, Rd_fun)
    n = N - 1
    y0 = np.concatenate([np.full(n, T0), np.full(n, C0)])
    print('ALE 移动网格谱求解 0→72h ...')
    sol = solve_ivp(fld.coupled_rhs, (0.0, T_FULL), y0, method='BDF',
                    t_eval=T_OUT60, rtol=1e-8, atol=1e-8)
    print(f'完成: nfev={sol.nfev}, 步数={len(sol.t)}')

    # 重建: 表面/轴心由消元恢复, 输出到 r 网格 (r ≤ R(t) 才有值)
    fld._surf_C = None
    fld._surf_T = None
    C_out, R_out = {}, {}
    x_out_static = None
    for i, tt in enumerate(sol.t):
        y = sol.y[:, i]
        yT, yC = y[:n], y[n:]
        Ts_guess = fld._surf_T if fld._surf_T is not None else T_env(tt)
        fld._set_time(tt)
        sC = fld._surf_moisture(tt, yC, Ts_guess)
        C_full = fld._full_field(yC, fld._axis_value(yC, sC), sC)
        T0v, sT = fld._surf_heat(tt, yT, C_full)
        fld._surf_C, fld._surf_T = sC, sT
        Rt = fld.R
        # 输出列: 0.0..1.9 (固定) + 药材表面
        vals = np.full(len(OUT_R_CM) + 1, np.nan)
        x_out = 1.0 - 2.0 * OUT_R_CM / 100.0 / Rt
        mask = OUT_R_CM / 100.0 <= Rt + 1e-14
        vals[:-1][mask] = barycentric_eval(fld.x, C_full, x_out[mask])
        vals[-1] = C_full[-1]                       # 表面列
        C_out[float(tt)] = vals
        R_out[float(tt)] = Rt

    # 判据: max C < 0.15 (材料内, 中心最严)
    centers = np.array([C_out[t][0] for t in T_OUT60])   # r=0 列
    meet = np.where(centers < THRESHOLD)[0]
    t_end = float(T_OUT60[meet[0]]) if len(meet) else None
    print(f'\n问题4 判据: t_end = {t_end:.0f} s = {t_end/3600:.4f} h' if t_end
          else '\n问题4: 72h 内未达标')
    if t_end:
        i_e = meet[0]
        print(f'  结束时刻: 中心C = {centers[i_e]:.6f}, R = {R_out[float(t_end)]*100:.4f} cm')
    print(f'  对比: v1 物质坐标 Q4 = 51.17 h; 附件2 收缩平台 = 70.0 h')

    # 表6 (每6h, 列 0,0.5,...,药材表面)
    print('\n表6  药材烘干过程的水分浓度 (kg/kg) [末列=药材表面]')
    for h in range(6, int(t_end // 3600) + 1, 6):
        tt = float(h * 3600)
        Rt = R_out[tt] * 100.0
        cols = np.arange(0.0, Rt - 1e-9, 0.5)
        vals = [float(np.interp(rc / 100.0, OUT_R_CM, C_out[tt][:-1]))
                for rc in cols]
        vals.append(C_out[tt][-1])
        print(f'{h:.0f}h\t' + '\t'.join(f'{rc:.1f}({v:.4f})' for rc, v in
                                        zip(list(cols) + [Rt], vals)))

    # result4.xlsx
    wb = Workbook()
    ws = wb.active
    ws.title = 'Sheet1'
    ws.cell(1, 1, '时间\\到药材中心的距离')
    for j, rr in enumerate(OUT_R_CM):
        ws.cell(1, j + 2, round(rr, 1)).number_format = '0.0'
    ws.cell(1, len(OUT_R_CM) + 2, '药材表面')
    n_rows = int(t_end / 60.0)
    for k in range(n_rows):
        tt = float(T_OUT60[k])
        ws.cell(k + 2, 1, int(tt))
        for j in range(len(OUT_R_CM) + 1):
            v = C_out[tt][j]
            if not np.isnan(v):
                cell = ws.cell(k + 2, j + 2, round(float(v), 4))
                cell.number_format = '0.0000'
    wb.save('results/result4.xlsx')
    print(f'\nresult4.xlsx 已保存 ({n_rows}行 × {len(OUT_R_CM)+1}列)')

    # 存档
    C_arr = np.array([C_out[t] for t in T_OUT60])
    R_arr = np.array([R_out[t] for t in T_OUT60])
    np.savez('results/q4_v2_full.npz', t=T_OUT60, C=C_arr, R=R_arr, t_end=t_end)
    print('results/q4_v2_full.npz 已保存')

    # GCL 检查 (几何守恒律: dV_j/dt = (2Ṙ/R)·V_j 解析恒等)
    print(f'[V3-GCL] 网格几何因子: dV/dt / V = 2Ṙ/R = {fld.gcl_check(T_FULL):.6e} /s')


if __name__ == '__main__':
    main()
