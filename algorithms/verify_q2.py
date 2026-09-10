# -*- coding: utf-8 -*-
"""
verify_q2.py — 问题2/3 的数值稳健性验证 (论文 6 节素材)
=========================================================
A. 时间步长减半 (0.125/0.5/2.5s 分档) → 对比 t_end 与 60s 网格解
B. 网格加倍 (N=800) → 对比 t_end 与 60s 网格解
"""
import sys, os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)      # 兄弟模块导入 (solver_core 等)
os.chdir(_ROOT)                # 数据/结果路径以项目根为基准


import numpy as np
from q2_solve import build_solver, make_env, T_END
from solver_core import CylinderDryer

R = 0.02
THRESHOLD = 0.15


def run_variant(dt_scale=1.0, N=400, T_end=T_END):
    T_env, C_env, _, _ = make_env()
    s = build_solver(T_env, C_env)
    s2 = CylinderDryer(R, N,
        rho_fun=s.rho_fun, cp_fun=s.cp_fun, k_fun=s.k_fun,
        D_fun=s.D_fun, dDdC_fun=s.dDdC_fun, dDdT_fun=s.dDdT_fun,
        h=s.h, beta=s.beta, T_env_fun=s.T_env_fun, C_env_fun=s.C_env_fun,
        T0=28.0, C0=2.55)
    phases = [(0.25 * dt_scale, 300.0), (1.0 * dt_scale, 14400.0), (5.0 * dt_scale, T_end)]
    snap60_t, snap60_C = [], []
    next_t60 = 60.0
    for dt, t_end in phases:
        first = True
        while s2.t < t_end - 1e-12:
            s2.step(min(dt, t_end - s2.t), force_restart=first)
            first = False
            if s2.t >= next_t60 - 1e-12:
                snap60_t.append(s2.t); snap60_C.append(s2.C.copy())
                next_t60 += 60.0
    t = np.array(snap60_t); C = np.array(snap60_C)
    Cmax = C.max(axis=1)
    meet = np.where(Cmax < THRESHOLD)[0]
    if len(meet) == 0:
        return None, C, None
    return float(t[meet[0]]), C, None


def main():
    base = np.load('results/q2_full.npz')
    Cbase, tbase = base['C'], base['t']
    # 基准 60s 网格上 N=800 解的插值对比 (N=800 节点更密, 用其自身节点对比中心/表面即可)
    print('基准 (dt=0.25/1/5s, N=400): t_end = 207960 s = 57.7667 h')
    for label, kw in (('dt减半', dict(dt_scale=0.5)), ('N=800', dict(N=800))):
        t_end, C, _ = run_variant(**kw)
        if t_end is None:
            print(f'{label}: 72h 内未达标')
            continue
        print(f'{label}: t_end = {t_end:.0f} s = {t_end/3600:.4f} h  '
              f'(Δ = {(t_end - 207960.0)/3600:.3f} h)')
        i24 = int(24 * 3600 / 60) - 1
        i48 = int(48 * 3600 / 60) - 1
        d24 = abs(C[i24, 0] - Cbase[i24, 0])
        d48 = abs(C[i48, 0] - Cbase[i48, 0])
        dS = abs(C[i24, -1] - Cbase[i24, -1])
        print(f'  表5对比 (中心): |ΔC| @24h = {d24:.2e}, @48h = {d48:.2e}, 表面@24h = {dS:.2e}')


if __name__ == '__main__':
    main()
