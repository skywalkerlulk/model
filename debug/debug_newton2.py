import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'algorithms'))  # 兄弟模块
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # 项目根
# -*- coding: utf-8 -*-
"""跟踪 N=800, t=14400, dt=5 (BE 重启) 的热场 Newton 迭代 (按 step() 实际路径)"""
import numpy as np
from q2_solve import build_solver, make_env
from solver_core import CylinderDryer, assemble_diffusion
from scipy.linalg import solve_banded

T_env, C_env, _, _ = make_env()
s = build_solver(T_env, C_env)
s2 = CylinderDryer(0.02, 800,
    rho_fun=s.rho_fun, cp_fun=s.cp_fun, k_fun=s.k_fun,
    D_fun=s.D_fun, dDdC_fun=s.dDdC_fun, dDdT_fun=s.dDdT_fun,
    h=25.0, beta=8e-7, T_env_fun=T_env, C_env_fun=C_env, T0=28.0, C0=2.55)
phases = [(0.25, 300.0), (1.0, 14400.0)]
for dt, t_end in phases:
    first = True
    while s2.t < t_end - 1e-12:
        s2.step(min(dt, t_end - s2.t), force_restart=first)
        first = False

m = s2.mesh
dt = 5.0
t_new = s2.t + dt
T = s2.T.copy()
C = s2.C.copy()
k = s2.k_fun(C)
rho, cp = s2.rho_fun(C), s2.cp_fun(C)
mass = rho * cp * m.V
coef = mass / dt
hist = mass * T / dt
u = T.copy()
print(f'T 场范围: [{T.min():.6f}, {T.max():.6f}], k 范围: [{k.min():.4f}, {k.max():.4f}]')
for it in range(50):
    lo, di, up, dlo, ddi, dup, bc = assemble_diffusion(m, k, np.zeros_like(k), s2.h, s2.T_env_fun(t_new))
    F = np.zeros_like(u)
    F[:-1] += up[:-1] * u[1:]
    F[1:] += lo[1:] * u[:-1]
    F += di * u
    F += bc
    G = coef * u - hist - F
    ab = np.zeros((3, len(u)))
    ab[0, 1:] = -dup[:-1]
    ab[1, :] = coef - ddi
    ab[2, :-1] = -dlo[1:]
    delta = solve_banded((1, 1), ab, -G)
    u = u + delta
    print(f'it={it:2d}  max|G|={np.max(np.abs(G)):.4e}  max|δ|={np.max(np.abs(delta)):.4e}  '
          f'T表面={u[-1]:.8f}')
    if np.max(np.abs(delta)) < 1e-12:
        print('收敛')
        break
