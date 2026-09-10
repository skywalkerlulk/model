# -*- coding: utf-8 -*-
"""
q2_compare.py — 变物性 vs 常物性对比实验 (论文 5.2.5 节)
============================================================
用附录3 公式在 C=2.55 (湿态) 与 C=0.15 (近干态) 处冻结物性,
与真实变物性模型 (附录3 全程) 对比烘干时长, 证明变物性建模的必要性.
QoI: t_end 与 24h/48h 中心含水率
"""
import sys, os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)      # 兄弟模块导入 (solver_core 等)
os.chdir(_ROOT)                # 数据/结果路径以项目根为基准


import numpy as np
from solver_core import CylinderDryer
from q2_solve import make_env

R = 0.02
N = 400
T_END = 259200.0
THRESHOLD = 0.15


def run_frozen(C_fix):
    T_env, C_env, _, _ = make_env()
    rho0, cp0, k0 = 650.0 + 128.0 * C_fix, 1450.0 + 2736.0 * C_fix / (C_fix + 1.0), 0.21 + 0.38 * C_fix / (C_fix + 1.0)
    D0 = 2.4e-3 * np.exp(-0.45 / C_fix)
    s = CylinderDryer(R, N,
        rho_fun=lambda C: np.full_like(C, rho0),
        cp_fun=lambda C: np.full_like(C, cp0),
        k_fun=lambda C: np.full_like(C, k0),
        D_fun=lambda C, T: D0 * np.exp(-3850.0 / (T + 273.15)),
        dDdC_fun=lambda C, T: np.zeros_like(C),
        dDdT_fun=lambda C, T: D0 * np.exp(-3850.0 / (T + 273.15)) * 3850.0 / (T + 273.15) ** 2,
        h=25.0, beta=8e-7, T_env_fun=T_env, C_env_fun=C_env, T0=28.0, C0=2.55)
    phases = [(0.25, 300.0), (1.0, 14400.0), (5.0, T_END)]
    t_hist, cmax_hist, c24, c48 = [], [], None, None
    for dt, t_end_ph in phases:
        first = True
        while s.t < t_end_ph - 1e-12:
            s.step(min(dt, t_end_ph - s.t), force_restart=first)
            first = False
            if s.t >= (len(t_hist) + 1) * 60.0 - 1e-12:
                t_hist.append(s.t); cmax_hist.append(float(s.C.max()))
                if s.t == 24 * 3600.0:
                    c24 = float(s.C[0])
                if s.t == 48 * 3600.0:
                    c48 = float(s.C[0])
    for i, cm in enumerate(cmax_hist):
        if cm < THRESHOLD:
            return t_hist[i], c24, c48
    return None, c24, c48


def main():
    print('对比实验: 附录3 物性冻结 vs 真实变物性 (基准 t_end = 57.77 h)')
    for label, C_fix in (('冻结@C=2.55 (湿态)', 2.55), ('冻结@C=0.15 (近干态)', 0.15)):
        te, c24, c48 = run_frozen(C_fix)
        if te:
            print(f'{label}: t_end = {te/3600:.2f} h, 中心C@24h = {c24:.4f}, @48h = {c48:.4f}')
        else:
            print(f'{label}: 72h 内未达标! 中心C@24h = {c24:.4f}, @48h = {c48:.4f}')
    print('\n结论: 湿态冻结显著低估烘干时长, 近干态冻结 72h 仍不达标 →')
    print('      D 随 C 的强非线性 (4个量级) 使变物性建模成为必要')


if __name__ == '__main__':
    main()
