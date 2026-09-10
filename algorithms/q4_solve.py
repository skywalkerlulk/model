# -*- coding: utf-8 -*-
"""
q4_solve.py — 问题4: 考虑收缩的移动边界烘干模型 (附录4 + 附件2)
============================================================
坐标: ξ = r/R(t) (仿射收缩物质坐标, R(t) 取自附件2 PCHIP 插值)
物性: 附录4  ρ=760+90C, cp=1850+2150·C/(C+1), k=0.12+0.20·C/(C+1),
      D=4.2e-4·exp(-0.30/C)·exp(-3850/T_K)
方程: d(V(t)·C)/dt = Σ D_face·2πξ/dξ·ΔC (乘积型 BDF2, 见 solver_core.ShinkDryer)
      d(ρcpV·T)/dt = Σ k_face·2πξ/dξ·ΔT + h·2πR(t)·(T_env-T)
输出: 表6 (每6h, 0..药材表面), result4.xlsx (60s × 0.1cm + 表面列, 至烘干结束)
"""
import sys, os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)      # 兄弟模块导入 (solver_core 等)
os.chdir(_ROOT)                # 数据/结果路径以项目根为基准


import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
from openpyxl import Workbook
from solver_core import ShrinkDryer
from q2_solve import make_env
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'PingFang SC', 'Heiti TC', 'STHeiti']
plt.rcParams['axes.unicode_minus'] = False

N = 400
T_END = 259200.0          # 72h; 若未达标则延长
THRESHOLD = 0.15
OUT_R_CM = np.arange(0.0, 1.9, 0.1)   # result4 固定列 0.0..1.9, 末列=药材表面


def make_R():
    """附件2 半径 (cm→m), PCHIP 插值 + 尾部外推保持常数"""
    df = pd.read_excel('A题/附件/附件2.xlsx')
    t = df['时间'].values.astype(float)
    R_cm = df['半径'].values.astype(float)
    p = PchipInterpolator(t, R_cm)
    R_tail = float(R_cm[-1])
    def R_fun(tt):
        if tt <= t[-1]:
            return float(p(tt)) / 100.0
        return R_tail / 100.0
    return R_fun, t, R_cm


def build_shrink_solver(T_env, C_env, R_fun):
    return ShrinkDryer(
        N,
        rho_fun=lambda C: 760.0 + 90.0 * C,
        cp_fun=lambda C: 1850.0 + 2150.0 * C / (C + 1.0),
        k_fun=lambda C: 0.12 + 0.20 * C / (C + 1.0),
        D_fun=lambda C, T: 4.2e-4 * np.exp(-0.30 / np.maximum(C, 1e-6))
                            * np.exp(-3850.0 / (T + 273.15)),
        dDdC_fun=lambda C, T: 4.2e-4 * np.exp(-0.30 / np.maximum(C, 1e-6))
                               * np.exp(-3850.0 / (T + 273.15))
                               * 0.30 / np.maximum(C, 1e-6) ** 2,
        dDdT_fun=lambda C, T: 4.2e-4 * np.exp(-0.30 / np.maximum(C, 1e-6))
                               * np.exp(-3850.0 / (T + 273.15))
                               * 3850.0 / (T + 273.15) ** 2,
        h=25.0, beta=8e-7,
        T_env_fun=T_env, C_env_fun=C_env, R_fun=R_fun,
        T0=28.0, C0=2.55, R0=0.02)


def main():
    T_env, C_env, _, _ = make_env()
    R_fun, t_R, R_cm = make_R()
    s = build_shrink_solver(T_env, C_env, R_fun)

    # ---------------- 分档推进 + 60s 快照 ----------------
    phases = [(0.25, 300.0), (1.0, 14400.0), (5.0, T_END)]
    snap_t, snap_C, snap_R, snap_T = [], [], [], []
    next_t60 = 60.0
    for dt, t_end in phases:
        first = True
        while s.t < t_end - 1e-12:
            s.step(min(dt, t_end - s.t), force_restart=first)
            first = False
            if s.t >= next_t60 - 1e-12:
                snap_t.append(s.t); snap_C.append(s.C.copy())
                snap_R.append(s.R); snap_T.append(s.T.copy())
                next_t60 += 60.0
    t = np.array(snap_t); C = np.array(snap_C); Rt = np.array(snap_R)
    print(f'60s 快照 {len(t)} 个; 终态 R = {s.R*100:.3f} cm, 中心 C = {C[-1,0]:.6f}')

    # ---- 判据: max_ξ C < 0.15 (最严点为 ξ=0) ----
    Cmax = C.max(axis=1)
    meet = np.where(Cmax < THRESHOLD)[0]
    if len(meet) == 0:
        print('72h 内未达标 — 需延长模拟!')
        return
    i_end = meet[0]
    t_end = float(t[i_end])
    print(f'烘干结束判定: t_end = {t_end:.0f} s = {t_end/3600:.4f} h')
    print(f'  对比: 问题3 (固定半径) = 57.77 h; 附件2 收缩平台 = 70.0 h')
    print(f'  结束时 R = {Rt[i_end]*100:.4f} cm, 中心 C = {C[i_end,0]:.6f}')

    # ---------------- 表6: 每6h, 列 0,0.5,...,药材表面 ----------------
    xi = np.linspace(0.0, 1.0, N + 1)
    print('\n表6  药材烘干过程的水分浓度 (kg/kg) [列: 0, 0.5, ... cm, 末列=药材表面]')
    for tt in range(6 * 3600, int(t_end) + 1, 6 * 3600):
        i = int(tt / 60.0) - 1
        R_cm_t = Rt[i] * 100.0
        cols = np.arange(0.0, R_cm_t, 0.5)          # 0, 0.5, ..., < R(t)
        vals = []
        for rc in cols:
            vals.append(float(np.interp(rc / 100.0 / Rt[i], xi, C[i])))
        vals.append(C[i, -1])                       # 表面
        cols = list(cols) + [R_cm_t]
        print(f'{tt/3600:.1f}h\t' + '\t'.join(
            [f'{rc:.1f}({v:.4f})' for rc, v in zip(cols, vals)]))

    # ---------------- result4.xlsx: 60s × (0.0..1.9 + 药材表面) ----------------
    wb = Workbook()
    ws = wb.active
    ws.title = 'Sheet1'
    ws.cell(1, 1, '时间\\到药材中心的距离')
    for j, rr in enumerate(OUT_R_CM):
        ws.cell(1, j + 2, round(rr, 1)).number_format = '0.0'
    ws.cell(1, len(OUT_R_CM) + 2, '药材表面')
    for k in range(i_end + 1):
        ws.cell(k + 2, 1, int(t[k]))
        R_k = Rt[k]
        for j, rr in enumerate(OUT_R_CM):
            r_m = rr / 100.0
            if r_m <= R_k + 1e-14:
                v = float(np.interp(r_m / R_k, xi, C[k]))
                cell = ws.cell(k + 2, j + 2, round(v, 4))
                cell.number_format = '0.0000'
            # r > R(t) 处无材料 → 留空
        cell = ws.cell(k + 2, len(OUT_R_CM) + 2, round(float(C[k, -1]), 4))
        cell.number_format = '0.0000'
    wb.save('results/result4.xlsx')
    print(f'\nresult4.xlsx 已保存 ({i_end+1}行 × {len(OUT_R_CM)+1}列, 60s间隔至烘干结束)')

    # ---------------- 图 ----------------
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    axes[0].plot(t_R / 3600.0, R_cm, '.-', ms=4, label='附件2 实测')
    axes[0].axvline(252000 / 3600.0, color='g', ls=':', lw=1.2, label='收缩平台 70h')
    axes[0].axvline(t_end / 3600.0, color='r', ls='--', lw=1.2, label=f'模型烘干结束 {t_end/3600:.1f}h')
    axes[0].set(xlabel='t / h', ylabel='R / cm', title='药材半径收缩曲线')
    axes[0].legend(); axes[0].grid(alpha=.3)
    axes[1].plot(t / 3600.0, C[:, -1], label='表面')
    axes[1].plot(t / 3600.0, C[:, 0], label='中心')
    axes[1].axhline(0.15, color='r', ls='--', lw=1, label='判据 0.15')
    axes[1].axvline(t_end / 3600.0, color='g', ls=':', lw=1.2, label=f'结束 {t_end/3600:.1f}h')
    axes[1].set(xlabel='t / h', ylabel='C / (kg·kg⁻¹)', title='水分浓度时间历程 (收缩模型)')
    axes[1].legend(); axes[1].grid(alpha=.3)
    fig.tight_layout(); fig.savefig('figures/q4_shrink.png', dpi=200)

    # 与问题3 对比
    q2 = np.load('results/q2_full.npz')
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(q2['t'] / 3600.0, q2['C'][:, 0], label='问题3 (固定半径) 中心')
    ax.plot(t / 3600.0, C[:, 0], label='问题4 (收缩) 中心')
    ax.axhline(0.15, color='r', ls='--', lw=1)
    ax.set(xlabel='t / h', ylabel='中心水分浓度 C(0,t) / (kg·kg⁻¹)',
           title='问题3 vs 问题4: 收缩对烘干时长的影响')
    ax.legend(); ax.grid(alpha=.3)
    fig.tight_layout(); fig.savefig('figures/q4_vs_q3.png', dpi=200)

    # 收缩-失水一致性检验 (理想收缩假设 vs 附件2 数据)
    xi2 = s._xi2
    V0 = np.sum(xi2) * np.pi * 0.02 ** 2
    Mw = np.array([np.sum(xi2 * C[k]) for k in range(len(t))]) * np.pi * 0.02 ** 2  # 水含量(干基参照)
    wl = 1.0 - Mw / Mw[0]                                # 模型失水分数
    vs = 1.0 - (Rt / 0.02) ** 2                          # 体积收缩分数
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(wl, vs, lw=1.5, label='模型: 体积收缩 vs 失水')
    ax.plot([0, 1], [0, 1], 'k--', lw=0.8, label='理想收缩线 (1:1)')
    ax.set(xlabel='模型失水分数 1-Mw/Mw0', ylabel='附件2 体积收缩分数 1-(R/R0)²',
           title='问题4: 收缩与失水的一致性检验')
    ax.legend(); ax.grid(alpha=.3)
    fig.tight_layout(); fig.savefig('figures/q4_consistency.png', dpi=200)
    print('图已保存: q4_shrink.png, q4_vs_q3.png, q4_consistency.png')
    np.savez('results/q4_full.npz', t=t, xi=xi, C=C, R=Rt)
    print('results/q4_full.npz 已保存')


if __name__ == '__main__':
    main()
