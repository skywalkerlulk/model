# -*- coding: utf-8 -*-
"""
q2_v2.py — 问题2/3 独立求解 (2026 CUMCM A题 v2 方案)
============================================================
模型: 附录3 变物性双向耦合 (VC-HTMT-1 谱配点版)
    ρ=650+128C, cp=1450+2736C/(C+1), k=0.21+0.38C/(C+1),
    D=2.4e-3·e^(−0.45/C)·e^(−3850/T_K)   [T 取开尔文]
耦合: coupled_rhs (谱边界消元 + 表面 2×2 联立 + 定点扫描)
环境: 一阶惯性模型 (env_model.py), t>4h 平台外推 (数据统计)
时间: solve_ivp BDF 自适应, 0→72h
输出: 表3/表4 (0.5-3h), result2.xlsx (1s×0.1cm×3h),
      全流程 60s 存档 (供问题3 判据扫描)
验证: 与 v1 FVM 参考解交叉对比 (V6)
"""

import sys, os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
sys.path.insert(0, os.path.join(_ROOT, 'algorithms_v2'))
sys.path.insert(0, os.path.join(_ROOT, 'algorithms'))
os.chdir(_ROOT)

import numpy as np
from scipy.integrate import solve_ivp
from openpyxl import Workbook

from env_model import load_env_models
from spectral_core import RadialField, barycentric_eval, cheb_coeffs

R = 0.02
N = int(sys.argv[1]) if len(sys.argv) > 1 else 96
# 谱模式数: N=48 欠分辨 (末位系数 1.3e-3); N=96 与独立 FVM 一致 6e-4;
# N=128 收敛 (t_end 57.83h vs FVM 57.77h) — 最终提交用 N=128
T_FULL = 259200.0            # 72 h
T_WIN3 = 10800.0             # 表3/表4 窗口 (3h)
T_OUT1 = np.arange(1.0, T_WIN3 + 1.0)          # result2: 1s 间隔
T_OUT60 = np.arange(60.0, T_FULL + 0.1, 60.0)  # 全流程: 60s 间隔
R_OUT_CM = np.arange(0.0, 2.01, 0.1)
TABLE_TIMES = [1800, 3600, 5400, 7200, 9000, 10800]
TABLE_R_CM = [0.0, 0.5, 1.0, 1.5, 2.0]
T0, C0 = 28.0, 2.55


def make_coupled_field(T_env, C_env):
    return RadialField(N, R,
        D_fun=lambda C, T: 2.4e-3 * np.exp(-0.45 / np.maximum(C, 1e-6))
                            * np.exp(-3850.0 / (T + 273.15)),
        dDdC_fun=lambda C, T: 2.4e-3 * np.exp(-0.45 / np.maximum(C, 1e-6))
                               * np.exp(-3850.0 / (T + 273.15))
                               * 0.45 / np.maximum(C, 1e-6) ** 2,
        dDdT_fun=lambda C, T: 2.4e-3 * np.exp(-0.45 / np.maximum(C, 1e-6))
                               * np.exp(-3850.0 / (T + 273.15))
                               * 3850.0 / (T + 273.15) ** 2,
        beta=8e-7, C_env_fun=C_env,
        rho_fun=lambda C: 650.0 + 128.0 * C,
        cp_fun=lambda C: 1450.0 + 2736.0 * C / (C + 1.0),
        k_fun=lambda C: 0.21 + 0.38 * C / (C + 1.0),
        dkdC_fun=lambda C: 0.38 / (C + 1.0) ** 2,
        h=25.0, T_env_fun=T_env)


def main():
    T_env, C_env, diag = load_env_models()
    fld = make_coupled_field(T_env, C_env)
    n = N - 1
    y0 = np.concatenate([np.full(n, T0), np.full(n, C0)])
    t_eval = np.unique(np.concatenate([T_OUT1, T_OUT60]))
    print('耦合谱求解 0→72h ...')
    sol = solve_ivp(fld.coupled_rhs, (0.0, T_FULL), y0, method='BDF',
                    t_eval=t_eval, rtol=1e-8, atol=1e-8)
    print(f'完成: nfev={sol.nfev}, 步数={len(sol.t)}')

    # 重建输出 (每个 t_eval 时刻恢复全场并插值到输出网格)
    # 关键: 按时间顺序处理, 表面值 Newton 续延沿物理支连续跟踪
    x_out = 1.0 - 2.0 * R_OUT_CM / 100.0 / R
    fld._surf_C = None
    fld._surf_T = None
    T_out, C_out = {}, {}
    full_last = None
    for i, tt in enumerate(sol.t):
        y = sol.y[:, i]
        yT, yC = y[:n], y[n:]
        # 用与 coupled_rhs 相同的边界消元恢复全场
        Ts_guess = fld._surf_T if fld._surf_T is not None else T_env(tt)
        sC = fld._surf_moisture(tt, yC, Ts_guess)
        C_full = fld._full_field(yC, fld._axis_value(yC, sC), sC)
        T0v, sT = fld._surf_heat(tt, yT, C_full)
        T_full = fld._full_field(yT, T0v, sT)
        fld._surf_C, fld._surf_T = sC, sT
        full_last = (T_full, C_full)
        T_out[float(tt)] = barycentric_eval(fld.x, T_full, x_out)
        C_out[float(tt)] = barycentric_eval(fld.x, C_full, x_out)
    T1 = np.array([T_out[float(t)] for t in T_OUT1])
    C1 = np.array([C_out[float(t)] for t in T_OUT1])
    # 60s 网格 (判据扫描 + 全流程存档)
    C60 = np.array([C_out[float(t)] for t in T_OUT60])
    T60 = np.array([T_out[float(t)] for t in T_OUT60])
    print('重建完成')

    # 表3/表4
    idx = {r: int(round(r / 0.1)) for r in TABLE_R_CM}
    print('\n表3  3小时内药材的温度 (°C)')
    for t in TABLE_TIMES:
        print(f'{t/3600:.1f}h\t' + '\t'.join(f'{T_out[float(t)][idx[r]]:.4f}' for r in TABLE_R_CM))
    print('\n表4  3小时内药材的水分浓度 (kg/kg)')
    for t in TABLE_TIMES:
        print(f'{t/3600:.1f}h\t' + '\t'.join(f'{C_out[float(t)][idx[r]]:.4f}' for r in TABLE_R_CM))

    # 判据: max C < 0.15 (60s 分辨率)
    Cmax = C60.max(axis=1)
    meet = np.where(Cmax < 0.15)[0]
    t_end = float(T_OUT60[meet[0]]) if len(meet) else None
    print(f'\n问题3 判据: t_end = {t_end:.0f} s = {t_end/3600:.4f} h' if t_end
          else '\n问题3: 72h 内未达标')
    if t_end:
        i_e = meet[0]
        print(f'  结束时刻: 中心C = {C60[i_e, 0]:.6f}, '
              f'上一时刻 = {C60[i_e-1, 0]:.6f}')

    # result2.xlsx
    wb = Workbook()
    for sheet_name, data in (('温度', T1), ('水分浓度', C1)):
        ws = wb.active if sheet_name == '温度' else wb.create_sheet()
        ws.title = sheet_name
        ws.cell(1, 1, '时间\\到药材中心的距离')
        for j, r in enumerate(R_OUT_CM):
            ws.cell(1, j + 2, round(r, 1)).number_format = '0.0'
        for i, t in enumerate(T_OUT1):
            ws.cell(i + 2, 1, int(t))
            for j, c in enumerate(data[i]):
                cell = ws.cell(i + 2, j + 2, round(float(c), 4))
                cell.number_format = '0.0000'
    wb.save('results/result2.xlsx')
    print('result2.xlsx 已保存 (v2)')

    # 72h 存档 (60s 网格, 供 q3_v2)
    np.savez('results/q2_v2_full.npz', t=T_OUT60, r_out=R_OUT_CM,
             C=C60, T=T60, t_end=t_end)
    print('results/q2_v2_full.npz 已保存')

    # 验证: BC 残差 + 谱系数
    Tf, Cf = full_last
    bcC = abs(2.4e-3 * np.exp(-0.45 / max(Cf[-1], 1e-6)) * np.exp(-3850 / (Tf[-1] + 273.15))
              * (fld.Dr @ Cf)[-1] + 8e-7 * (Cf[-1] - C_env(T_FULL)))
    print(f'[V4] 终态表面BC残差 = {bcC:.2e}')
    print(f'[V5] 终态末位谱系数: |a_C| = {fld.tail_coeff(Cf):.2e}, |a_T| = {fld.tail_coeff(Tf):.2e}')


if __name__ == '__main__':
    main()
