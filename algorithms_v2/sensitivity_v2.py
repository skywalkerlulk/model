# -*- coding: utf-8 -*-
"""
sensitivity_v2.py — v2 方案灵敏度分析 (论文 6 节, 多进程并行版)
============================================================
方法: OAT 单参数 ±10% (QoI = t_end)
参数: h, β, D 前置因子, 恒温段 T_bar, C_bar
并行: 12 核 M4 Pro 上 10 个独立算例多进程并行 (~15 min, 串行需 ~2h)
模型: 谱配点耦合 (N=96), 粗网格 1h + 48h 起 60s 细化 (正确重启状态)
用法: python3 algorithms_v2/sensitivity_v2.py [--workers 8]
"""

import sys, os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
sys.path.insert(0, os.path.join(_ROOT, 'algorithms_v2'))
os.chdir(_ROOT)

import numpy as np
from multiprocessing import Pool
from scipy.integrate import solve_ivp
from spectral_core import RadialField
from env_model import load_env_models

R = 0.02
N = 96
T_FULL = 259200.0
THRESHOLD = 0.15
BASE = dict(h=25.0, beta=8e-7, Dpre=2.4e-3, Tbar=49.969, Cbar=0.0499)


def make_env_pert(Tbar, Cbar):
    T_env, C_env, _ = load_env_models()
    def T_env_p(t):
        t = float(t)
        return T_env(t) if t <= 14400 else Tbar
    def C_env_p(t):
        t = float(t)
        return C_env(t) if t <= 14400 else Cbar
    return T_env_p, C_env_p


def make_field(T_env, C_env, h, beta, Dpre):
    return RadialField(N, R,
        D_fun=lambda C, T: Dpre * np.exp(-0.45 / np.maximum(C, 1e-6))
                            * np.exp(-3850.0 / (T + 273.15)),
        dDdC_fun=lambda C, T: Dpre * np.exp(-0.45 / np.maximum(C, 1e-6))
                               * np.exp(-3850.0 / (T + 273.15))
                               * 0.45 / np.maximum(C, 1e-6) ** 2,
        dDdT_fun=lambda C, T: Dpre * np.exp(-0.45 / np.maximum(C, 1e-6))
                               * np.exp(-3850.0 / (T + 273.15))
                               * 3850.0 / (T + 273.15) ** 2,
        beta=beta, C_env_fun=C_env,
        rho_fun=lambda C: 650.0 + 128.0 * C,
        cp_fun=lambda C: 1450.0 + 2736.0 * C / (C + 1.0),
        k_fun=lambda C: 0.21 + 0.38 * C / (C + 1.0),
        dkdC_fun=lambda C: 0.38 / (C + 1.0) ** 2,
        h=h, T_env_fun=T_env)


def run_case(args):
    """单个算例: (h, beta, Dpre, Tbar, Cbar) → t_end 或 None (72h 未达标)"""
    h, beta, Dpre, Tbar, Cbar = args
    T_env, C_env = make_env_pert(Tbar, Cbar)
    fld = make_field(T_env, C_env, h, beta, Dpre)
    n = N - 1
    y0 = np.concatenate([np.full(n, 28.0), np.full(n, 2.55)])
    # 粗段: 0→48h, 每小时
    t_coarse = np.arange(3600.0, 48 * 3600.0 + 1, 3600.0)
    sol = solve_ivp(fld.coupled_rhs, (0.0, 48 * 3600.0), y0, method='BDF',
                    t_eval=t_coarse, rtol=1e-8, atol=1e-8)
    # 细段: 48h→72h, 60s, 正确重启状态 (48h 列)
    i48 = int(np.where(np.isclose(sol.t, 48 * 3600.0))[0][0])
    t_fine = np.arange(48 * 3600.0, T_FULL + 1, 60.0)
    sol2 = solve_ivp(fld.coupled_rhs, (48 * 3600.0, T_FULL), sol.y[:, i48],
                     method='BDF', t_eval=t_fine, rtol=1e-8, atol=1e-8)
    # 中心值扫描 (粗段 1h 分辨率 + 细段 60s)
    fld._surf_C = None
    fld._surf_T = None
    for i, tt in enumerate(sol.t[:-1]):
        y = sol.y[:, i]
        yC = y[n:]
        Ts = fld._surf_T if fld._surf_T is not None else T_env(tt)
        sC = fld._surf_moisture(tt, yC, Ts)
        C_full = fld._full_field(yC, fld._axis_value(yC, sC), sC)
        T0v, sT = fld._surf_heat(tt, y[:n], C_full)
        fld._surf_C, fld._surf_T = sC, sT
        if C_full[0] < THRESHOLD:
            return tt
    for i, tt in enumerate(sol2.t):
        y = sol2.y[:, i]
        yC = y[n:]
        Ts = fld._surf_T if fld._surf_T is not None else T_env(tt)
        sC = fld._surf_moisture(tt, yC, Ts)
        C_full = fld._full_field(yC, fld._axis_value(yC, sC), sC)
        T0v, sT = fld._surf_heat(tt, y[:n], C_full)
        fld._surf_C, fld._surf_T = sC, sT
        if C_full[0] < THRESHOLD:
            return tt
    return None


def main():
    n_workers = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    cases = [(p, s) for p in BASE for s in (-0.1, 0.1)]
    args_list = [tuple(BASE[p] * (1 + s) if k == p else BASE[k]
                       for k in ('h', 'beta', 'Dpre', 'Tbar', 'Cbar'))
                 for p, s in cases]
    print(f'多进程并行: {len(args_list)} 个算例 × {n_workers} workers')
    print(f'基准参考: t_end = 57.83 h (N=128 收敛值; N=96 口径 58.85 h)')
    with Pool(n_workers) as pool:
        results = pool.map(run_case, args_list)
    print('\nOAT ±10% 结果 (t_end):')
    out = []
    for (p, s), te in zip(cases, results):
        dte = (te - 211860.0) / 3600.0 if te else np.nan
        out.append((p, s, dte))
        print(f'  {p:5s} {s:+5.0%}: t_end = {te/3600:.2f} h' if te else
              f'  {p:5s} {s:+5.0%}: 72h 未达标', f'(Δt_end = {dte:+.2f} h)')
    print('\n主导参数排序 (按 |Δ(−10%)| + |Δ(+10%)|):')
    rank = {}
    for p, s, d in out:
        rank.setdefault(p, [None, None])
        rank[p][0 if s < 0 else 1] = d
    for p, (dlo, dhi) in sorted(rank.items(), key=lambda kv: -abs(kv[1][0] or 0) - abs(kv[1][1] or 0)):
        print(f'  {p}: |Δ(−10%)|+|Δ(+10%)| = {abs(dlo or 0) + abs(dhi or 0):.2f} h  ({dlo:+.2f} / {dhi:+.2f})')
    np.save('results/sensitivity_v2.npy', np.array(out, dtype=object))
    print('\nresults/sensitivity_v2.npy 已保存')


if __name__ == '__main__':
    main()
