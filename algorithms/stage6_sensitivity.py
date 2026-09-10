# -*- coding: utf-8 -*-
"""
stage6_sensitivity.py — 阶段6: 全局灵敏度分析 (问题2/3 模型, QoI = 烘干时长 t_end)
============================================================
方法: ① LHS 多变量联合扰动 (±10%, 10 样本) → t_end 置信区间
     ② OAT 单参数 ±10% → 主导参数排序 (tornado)
参数: h(25), β(8e-7), D前置因子(2.4e-3), 恒温段 T_bar(49.969), C_bar(0.0500)
"""
import sys, os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)      # 兄弟模块导入 (solver_core 等)
os.chdir(_ROOT)                # 数据/结果路径以项目根为基准


import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy.stats import qmc
from solver_core import CylinderDryer

R = 0.02
N = 400
T_END = 259200.0
THRESHOLD = 0.15

BASE = {'h': 25.0, 'beta': 8e-7, 'Dpre': 2.4e-3, 'Tbar': 49.969, 'Cbar': 0.0500}


def make_env_perturbed(Tbar, Cbar):
    df = pd.read_excel('A题/附件/附件1.xlsx')
    t = df['时间'].values.astype(float)
    Tp = PchipInterpolator(t, df['温度'].values)
    Cp = PchipInterpolator(t, df['水分浓度'].values)
    def T_env(tt):
        return float(np.where(tt <= 14400.0, Tp(min(tt, 14400.0)), Tbar))
    def C_env(tt):
        return float(np.where(tt <= 14400.0, Cp(min(tt, 14400.0)), Cbar))
    return T_env, C_env


def run_t_end(h, beta, Dpre, Tbar, Cbar):
    """给定参数运行全流程, 返回 t_end (60s 分辨率, 未达标返回 None)"""
    T_env, C_env = make_env_perturbed(Tbar, Cbar)
    s = CylinderDryer(R, N,
        rho_fun=lambda C: 650.0 + 128.0 * C,
        cp_fun=lambda C: 1450.0 + 2736.0 * C / (C + 1.0),
        k_fun=lambda C: 0.21 + 0.38 * C / (C + 1.0),
        D_fun=lambda C, T: Dpre * np.exp(-0.45 / np.maximum(C, 1e-6))
                            * np.exp(-3850.0 / (T + 273.15)),
        dDdC_fun=lambda C, T: Dpre * np.exp(-0.45 / np.maximum(C, 1e-6))
                               * np.exp(-3850.0 / (T + 273.15))
                               * 0.45 / np.maximum(C, 1e-6) ** 2,
        dDdT_fun=lambda C, T: Dpre * np.exp(-0.45 / np.maximum(C, 1e-6))
                               * np.exp(-3850.0 / (T + 273.15))
                               * 3850.0 / (T + 273.15) ** 2,
        h=h, beta=beta, T_env_fun=T_env, C_env_fun=C_env, T0=28.0, C0=2.55)
    phases = [(0.25, 300.0), (1.0, 14400.0), (5.0, T_END)]
    Cmax_hist = []
    t_hist = []
    for dt, t_end_ph in phases:
        first = True
        while s.t < t_end_ph - 1e-12:
            s.step(min(dt, t_end_ph - s.t), force_restart=first)
            first = False
            if s.t >= (len(t_hist) + 1) * 60.0 - 1e-12:
                t_hist.append(s.t); Cmax_hist.append(float(s.C.max()))
    for i, cm in enumerate(Cmax_hist):
        if cm < THRESHOLD:
            return t_hist[i]
    return None


def main():
    # ① OAT ±10%
    print('① 单参数扰动 (OAT ±10%), QoI = t_end')
    results = []
    base_t = 207960.0
    for p, v in BASE.items():
        lo = run_t_end(**{**BASE, p: v * 0.9})
        hi = run_t_end(**{**BASE, p: v * 1.1})
        dlo = (lo - base_t) / 3600.0 if lo else np.nan
        dhi = (hi - base_t) / 3600.0 if hi else np.nan
        results.append((p, v, dlo, dhi))
        print(f'  {p}: -10% → Δt_end = {dlo:+.2f} h, +10% → Δt_end = {dhi:+.2f} h')
    print('\n主导参数排序 (按 |Δ| 之和):')
    for p, v, dlo, dhi in sorted(results, key=lambda x: -abs(x[2]) - abs(x[3])):
        print(f'  {p}: |Δ(-10%)|+|Δ(+10%)| = {abs(dlo)+abs(dhi):.2f} h')

    # ② LHS 联合扰动 ±10%, 10 样本
    print('\n② LHS 多变量联合扰动 (±10%, 10 样本)')
    sampler = qmc.LatinHypercube(d=5, seed=42)
    unit = sampler.random(10)
    lo_b = np.array([0.9 * BASE[p] for p in ('h', 'beta', 'Dpre', 'Tbar', 'Cbar')])
    hi_b = np.array([1.1 * BASE[p] for p in ('h', 'beta', 'Dpre', 'Tbar', 'Cbar')])
    samples = qmc.scale(unit, lo_b, hi_b)
    t_ends = []
    for k, srow in enumerate(samples):
        te = run_t_end(h=srow[0], beta=srow[1], Dpre=srow[2], Tbar=srow[3], Cbar=srow[4])
        t_ends.append(te)
        print(f'  样本{k+1:2d}: t_end = {te/3600 if te else np.nan:8.2f} h'
              f'  (h={srow[0]:.1f}, β={srow[1]:.2e}, Dpre={srow[2]:.2e}, Tbar={srow[3]:.1f}, Cbar={srow[4]:.4f})')
    t_ends = np.array(t_ends) / 3600.0
    print(f'\n  LHS 统计: t_end ∈ [{t_ends.min():.2f}, {t_ends.max():.2f}] h, '
          f'均值 {t_ends.mean():.2f} h, 标准差 {t_ends.std():.2f} h')
    print(f'  结论: 烘干时长对参数联合扰动±10%的稳健区间 = {t_ends.mean():.2f} ± {2*t_ends.std():.2f} h (95%)')


if __name__ == '__main__':
    main()
