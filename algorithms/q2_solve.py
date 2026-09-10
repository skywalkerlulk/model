# -*- coding: utf-8 -*-
"""
q2_solve.py — 问题2: 全流程变物性双向耦合烘干模型 (附录3)
============================================================
物性: ρ(C)=650+128C, cp(C)=1450+2736·C/(C+1), k(C)=0.21+0.38·C/(C+1),
      D(C,T)=2.4e-3·exp(-0.45/C)·exp(-3850/T_K)  [T 取开尔文]
环境: 0-14400s 附件1 PCHIP; t>14400s 平台外推 (数据统计: 49.97±0.17°C, 0.0500 kg/kg)
时间: 变步长 BDF2 — 0.25s (0-5min) / 1s (5min-4h) / 5s (4h-72h), 切换点 BE 重启
输出: 表3/表4 (0.5-3h), result2.xlsx (1s×0.1cm, 前3h), 全流程 npz (供问题3)
"""
import sys, os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)      # 兄弟模块导入 (solver_core 等)
os.chdir(_ROOT)                # 数据/结果路径以项目根为基准


import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
from openpyxl import Workbook
from solver_core import CylinderDryer
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'PingFang SC', 'Heiti TC', 'STHeiti']
plt.rcParams['axes.unicode_minus'] = False

R = 0.02
N = 400
T_END = 259200.0          # 72 h (附件2 半径平台的终点, 必要时延长)
TABLE_TIMES = [1800, 3600, 5400, 7200, 9000, 10800]
TABLE_R_CM = [0.0, 0.5, 1.0, 1.5, 2.0]
OUT_R_CM = np.arange(0.0, 2.01, 0.1)


def make_env():
    """附件1 + 平台外推的烘房环境"""
    df = pd.read_excel('A题/附件/附件1.xlsx')
    t = df['时间'].values.astype(float)
    Tp = PchipInterpolator(t, df['温度'].values)
    Cp = PchipInterpolator(t, df['水分浓度'].values)
    T_bar = float(df['温度'].values[t >= 7200].mean())       # 恒温段均值
    C_bar = float(df['水分浓度'].values[t >= 7200].mean())
    def T_env(tt):
        return float(np.where(tt <= 14400.0, Tp(min(tt, 14400.0)), T_bar))
    def C_env(tt):
        return float(np.where(tt <= 14400.0, Cp(min(tt, 14400.0)), C_bar))
    return T_env, C_env, T_bar, C_bar


def build_solver(T_env, C_env):
    """附录3 变物性双向耦合"""
    return CylinderDryer(
        R, N,
        rho_fun=lambda C: 650.0 + 128.0 * C,
        cp_fun=lambda C: 1450.0 + 2736.0 * C / (C + 1.0),
        k_fun=lambda C: 0.21 + 0.38 * C / (C + 1.0),
        D_fun=lambda C, T: 2.4e-3 * np.exp(-0.45 / np.maximum(C, 1e-6))
                            * np.exp(-3850.0 / (T + 273.15)),
        dDdC_fun=lambda C, T: 2.4e-3 * np.exp(-0.45 / np.maximum(C, 1e-6))
                               * np.exp(-3850.0 / (T + 273.15))
                               * 0.45 / np.maximum(C, 1e-6) ** 2,
        dDdT_fun=lambda C, T: 2.4e-3 * np.exp(-0.45 / np.maximum(C, 1e-6))
                               * np.exp(-3850.0 / (T + 273.15))
                               * 3850.0 / (T + 273.15) ** 2,
        h=25.0, beta=8e-7,
        T_env_fun=T_env, C_env_fun=C_env,
        T0=28.0, C0=2.55)


def main():
    T_env, C_env, T_bar, C_bar = make_env()
    print(f'恒温段外推: T_env = {T_bar:.3f} °C, C_env = {C_bar:.4f} kg/kg (t > 4 h)')
    s = build_solver(T_env, C_env)

    # ---------------- 变步长分档推进 + 快照收集 ----------------
    phases = [(0.25, 300.0), (1.0, 14400.0), (5.0, T_END)]
    snap1_t, snap1_T, snap1_C = [], [], []        # 1s 网格 (前 3h, result2)
    snap60_t, snap60_T, snap60_C = [], [], []     # 60s 网格 (全流程, 问题3)
    next_t1, next_t60 = 1.0, 60.0
    phase_start = 0.0
    for dt, t_end in phases:
        first = True
        while s.t < t_end - 1e-12:
            s.step(min(dt, t_end - s.t), force_restart=first)
            first = False
            if s.t >= next_t1 - 1e-12 and next_t1 <= 10800.0:
                snap1_t.append(s.t); snap1_T.append(s.T.copy()); snap1_C.append(s.C.copy())
                next_t1 += 1.0
            if s.t >= next_t60 - 1e-12:
                snap60_t.append(s.t); snap60_T.append(s.T.copy()); snap60_C.append(s.C.copy())
                next_t60 += 60.0
        phase_start = t_end
        print(f'  阶段完成: 0-{t_end}s (dt={dt}s), GS扫描均值≈{2 if dt<=1 else 2}次, t={s.t:.0f}s')

    t1 = np.array(snap1_t); T1 = np.array(snap1_T); C1 = np.array(snap1_C)
    t60 = np.array(snap60_t); T60 = np.array(snap60_T); C60 = np.array(snap60_C)
    print(f'快照: 1s网格 {len(t1)} 个, 60s网格 {len(t60)} 个')
    print(f'终态: 中心 C = {C60[-1,0]:.6f}, 表面 C = {C60[-1,-1]:.6f}, 表面 T = {T60[-1,-1]:.3f}°C')

    # ---------------- 表3 / 表4 ----------------
    idx = (np.array(TABLE_R_CM) / 100.0 / (R / N)).astype(int)
    i_t = {t: int(t) - 1 for t in TABLE_TIMES}
    print('\n表3  3小时内药材的温度 (°C)')
    print('时间/h\t到药材中心的距离/cm: ' + '\t'.join(f'{r:.1f}' for r in TABLE_R_CM))
    for t in TABLE_TIMES:
        print(f'{t/3600:.1f}\t' + '\t'.join(f'{v:.4f}' for v in T1[i_t[t]][idx]))
    print('\n表4  3小时内药材的水分浓度 (kg/kg)')
    print('时间/h\t到药材中心的距离/cm: ' + '\t'.join(f'{r:.1f}' for r in TABLE_R_CM))
    for t in TABLE_TIMES:
        print(f'{t/3600:.1f}\t' + '\t'.join(f'{v:.4f}' for v in C1[i_t[t]][idx]))

    # ---------------- result2.xlsx (前3h, 1s×0.1cm) ----------------
    out_idx = (OUT_R_CM / 100.0 / (R / N)).astype(int)
    wb = Workbook()
    for sheet_name, data in (('温度', T1), ('水分浓度', C1)):
        ws = wb.active if sheet_name == '温度' else wb.create_sheet()
        ws.title = sheet_name
        ws.cell(1, 1, '时间\\到药材中心的距离')
        for j, r in enumerate(OUT_R_CM):
            ws.cell(1, j + 2, round(r, 1)).number_format = '0.0'
        for i, tt in enumerate(t1):
            ws.cell(i + 2, 1, int(tt))
            for j, c in enumerate(data[i][out_idx]):
                cell = ws.cell(i + 2, j + 2, round(float(c), 4))
                cell.number_format = '0.0000'
    wb.save('results/result2.xlsx')
    print('\nresult2.xlsx 已保存 (10800行 × 21列 × 2 sheet, 4位小数)')

    # ---------------- 全流程存档 (供问题3) ----------------
    np.savez('results/q2_full.npz', t=t60, r=s.mesh.r, T=T60, C=C60)
    print('results/q2_full.npz 已保存 (60s 网格全流程)')

    # ---------------- 图 ----------------
    r_cm = s.mesh.r * 100.0
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for t in TABLE_TIMES:
        axes[0].plot(r_cm, T1[i_t[t]], label=f'{t/3600:.1f} h')
        axes[1].plot(r_cm, C1[i_t[t]], label=f'{t/3600:.1f} h')
    axes[0].set(xlabel='r / cm', ylabel='T / °C', title='问题2: 温度分布 (0-3h)')
    axes[1].set(xlabel='r / cm', ylabel='C / (kg·kg⁻¹)', title='问题2: 水分浓度分布 (0-3h)')
    for ax in axes: ax.grid(alpha=.3); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig('figures/q2_profiles.png', dpi=200)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    axes[0].plot(t60 / 3600.0, C60[:, -1], label='表面 r=2cm')
    axes[0].plot(t60 / 3600.0, C60[:, 0], label='中心 r=0')
    axes[0].axhline(0.15, color='r', ls='--', lw=1, label='烘干判据 0.15')
    axes[0].set(xlabel='t / h', ylabel='C / (kg·kg⁻¹)', title='水分浓度时间历程 (72h)')
    axes[0].legend(); axes[0].grid(alpha=.3)
    axes[1].plot(t60 / 3600.0, T60[:, -1], label='表面')
    axes[1].plot(t60 / 3600.0, T60[:, 0], label='中心')
    axes[1].set(xlabel='t / h', ylabel='T / °C', title='温度时间历程 (72h)')
    axes[1].legend(); axes[1].grid(alpha=.3)
    fig.tight_layout(); fig.savefig('figures/q2_history.png', dpi=200)

    # D(C,T) 随水分的变化 (解释减速干燥)
    C_grid = np.linspace(0.05, 2.55, 200)
    D_vals = 2.4e-3 * np.exp(-0.45 / C_grid) * np.exp(-3850.0 / (T_bar + 273.15))
    fig, ax = plt.subplots(figsize=(5.5, 3.5))
    ax.semilogy(C_grid, D_vals)
    ax.set(xlabel='C / (kg·kg⁻¹)', ylabel='D / (m²·s⁻¹)',
           title=f'扩散系数随水分浓度变化 (T={T_bar:.0f}°C)')
    ax.grid(alpha=.3, which='both')
    fig.tight_layout(); fig.savefig('figures/q2_D_curve.png', dpi=200)
    print('图已保存: q2_profiles.png, q2_history.png, q2_D_curve.png')


if __name__ == '__main__':
    main()
