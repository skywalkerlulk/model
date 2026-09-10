# -*- coding: utf-8 -*-
"""validation_figs.py — 论文 5.1.3 节验证图 (解析解对比 + 收敛曲线)"""
import sys, os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)      # 兄弟模块导入 (solver_core 等)
os.chdir(_ROOT)                # 数据/结果路径以项目根为基准

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import j0, j1, jn_zeros
from scipy.optimize import brentq
from validate_solver import series_cyl, robin_roots, q1_solver
from q1_solve import T_END

plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'PingFang SC', 'Heiti TC', 'STHeiti']
plt.rcParams['axes.unicode_minus'] = False

R = 0.02
k, rho, cp, h = 0.36, 820.0, 2600.0, 25.0
alpha = k / (rho * cp)

# ---- 图1: 数值解 vs 解析解 (阶跃边界, 常物性) ----
from solver_core import CylinderDryer
s = CylinderDryer(R, 400,
    rho_fun=lambda C: np.full_like(C, rho), cp_fun=lambda C: np.full_like(C, cp),
    k_fun=lambda C: np.full_like(C, k),
    D_fun=lambda C, T: np.full_like(C, 5e-9), dDdC_fun=lambda C, T: np.zeros_like(C),
    dDdT_fun=lambda C, T: np.zeros_like(C),
    h=h, beta=8e-7, T_env_fun=lambda t: 50.0, C_env_fun=lambda t: 0.02,
    T0=28.0, C0=2.55)
t_targets = [100.0, 600.0, 1800.0]
snaps = []
for tt in t_targets:
    while s.t < tt - 1e-12:
        s.step(min(0.25, tt - s.t), force_restart=False)
    snaps.append((s.t, s.T.copy()))
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
for (tt, Tnum) in snaps:
    Tref = series_cyl(s.mesh.r, tt, R, alpha, k, h, 28.0, 50.0)
    axes[0].plot(s.mesh.r * 100, Tnum, 'o', ms=3, label=f'数值 t={tt:.0f}s')
    axes[0].plot(s.mesh.r * 100, Tref, '-', lw=1, label=f'解析 t={tt:.0f}s')
    axes[1].semilogy(s.mesh.r * 100, np.abs(Tnum - Tref), label=f't={tt:.0f}s')
axes[0].set(xlabel='r / cm', ylabel='T / °C', title='数值解与解析解对比 (阶跃边界)')
axes[0].legend(fontsize=7); axes[0].grid(alpha=.3)
axes[1].set(xlabel='r / cm', ylabel='|误差| / °C', title='误差分布 (对数坐标)')
axes[1].legend(fontsize=7); axes[1].grid(alpha=.3, which='both')
fig.tight_layout(); fig.savefig('figures/q1_validation.png', dpi=200)

# ---- 图2: 网格收敛曲线 (真实 Q1 参数, t=1800s) ----
r_out = np.array([0.0, 0.5, 1.0, 1.5, 2.0]) / 100.0
errs_T, errs_C = [], []
Ns = [200, 400, 800, 1600]
prev_T, prev_C = None, None
for N in Ns:
    sq = q1_solver(N, 0.25)
    sq.run(T_END, 0.25)
    idx = (r_out / (0.02 / N)).astype(int)
    Tq, Cq = sq.T[idx], sq.C[idx]
    if prev_T is not None:
        errs_T.append(np.max(np.abs(Tq - prev_T)))
        errs_C.append(np.max(np.abs(Cq - prev_C)))
    prev_T, prev_C = Tq.copy(), Cq.copy()
fig, ax = plt.subplots(figsize=(5.5, 3.8))
ax.loglog(Ns[1:], errs_T, 'o-', label='温度')
ax.loglog(Ns[1:], errs_C, 's-', label='水分浓度')
ax.loglog(Ns[1:], [1e-4 * (400 / N) ** 2 for N in Ns[1:]], 'k--', lw=0.8, label='二阶参考线')
ax.set(xlabel='网格数 N', ylabel='相邻网格最大差', title='网格收敛性 (t=1800s)')
ax.legend(); ax.grid(alpha=.3, which='both')
fig.tight_layout(); fig.savefig('figures/q1_convergence.png', dpi=200)
print('图已保存: q1_validation.png, q1_convergence.png')
