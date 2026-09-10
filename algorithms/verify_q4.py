# -*- coding: utf-8 -*-
"""
verify_q4.py — 问题4 收缩模型的数值稳健性验证
============================================
A. 时间步长减半 (0.125/0.5/2.5s 分档) → 对比 t_end
B. 网格加倍 (N=800) → 对比 t_end 与 表6 中心值
"""
import sys, os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)      # 兄弟模块导入 (solver_core 等)
os.chdir(_ROOT)                # 数据/结果路径以项目根为基准


import numpy as np
from q4_solve import make_R, build_shrink_solver, make_env, THRESHOLD

T_END = 259200.0


def run_variant(dt_scale=1.0, N=400):
    T_env, C_env, _, _ = make_env()
    R_fun, _, _ = make_R()
    s = build_shrink_solver(T_env, C_env, R_fun)
    if N != 400:
        from solver_core import ShrinkDryer
        s = ShrinkDryer(N,
            rho_fun=s.rho_fun, cp_fun=s.cp_fun, k_fun=s.k_fun,
            D_fun=s.D_fun, dDdC_fun=s.dDdC_fun, dDdT_fun=s.dDdT_fun,
            h=25.0, beta=8e-7, T_env_fun=T_env, C_env_fun=C_env,
            R_fun=R_fun, T0=28.0, C0=2.55, R0=0.02)
    phases = [(0.25 * dt_scale, 300.0), (1.0 * dt_scale, 14400.0), (5.0 * dt_scale, T_END)]
    t_hist, cmax_hist, c_hist = [], [], []
    for dt, t_end_ph in phases:
        first = True
        while s.t < t_end_ph - 1e-12:
            s.step(min(dt, t_end_ph - s.t), force_restart=first)
            first = False
            if s.t >= (len(t_hist) + 1) * 60.0 - 1e-12:
                t_hist.append(s.t)
                cmax_hist.append(float(s.C.max()))
                c_hist.append(float(s.C[0]))       # 中心
    for i, cm in enumerate(cmax_hist):
        if cm < THRESHOLD:
            return t_hist[i], t_hist, c_hist
    return None, t_hist, c_hist


def main():
    base_t, base_tv, base_c = run_variant(1.0, 400)
    print(f'基准 (dt=0.25/1/5s, N=400): t_end = {base_t:.0f} s = {base_t/3600:.4f} h')
    for label, kw in (('dt减半', dict(dt_scale=0.5)), ('N=800', dict(N=800))):
        te, tv, cv = run_variant(**kw)
        print(f'{label}: t_end = {te:.0f} s = {te/3600:.4f} h  (Δ = {(te-base_t)/3600:.3f} h)')
        i24 = int(24 * 3600 / 60) - 1
        i48 = int(48 * 3600 / 60) - 1
        print(f'  中心C对比: |Δ| @24h = {abs(cv[i24]-base_c[i24]):.2e}, @48h = {abs(cv[i48]-base_c[i48]):.2e}')


if __name__ == '__main__':
    main()
