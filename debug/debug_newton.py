import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'algorithms'))  # 兄弟模块
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # 项目根
# -*- coding: utf-8 -*-
"""解剖 N=800, t=14400, dt=5 (BE 重启) 的水分 Newton 迭代"""
import numpy as np
from q2_solve import build_solver, make_env
from solver_core import CylinderDryer, assemble_diffusion, _apply_face_deriv
from scipy.linalg import solve_banded

T_env, C_env, _, _ = make_env()
s = build_solver(T_env, C_env)

# 先用 N=800 跑到 t=14400 (与失败路径一致)
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

print(f'入口状态 (t={s2.t}s): C 表面={s2.C[-1]:.6f}, 中心={s2.C[0]:.6f}, T 表面={s2.T[-1]:.4f}')
m = s2.mesh
dt = 5.0
t_new = s2.t + dt
C = s2.C.copy()
D = s2.D_fun(C, s2.T)
dDdC = s2.dDdC_fun(C, s2.T)
print(f'D 范围: [{D.min():.3e}, {D.max():.3e}]')
print(f'表面 2 节点: C[-1]={C[-1]:.6f} (D={D[-1]:.3e}), C[-2]={C[-2]:.6f} (D={D[-2]:.3e})')

# 手动 BE Newton
coef = m.V / dt
hist = m.V * C / dt
u = C.copy()
for it in range(50):
    D = s2.D_fun(u, s2.T)
    dDdC = s2.dDdC_fun(u, s2.T)
    lo, di, up, dlo, ddi, dup, bc = assemble_diffusion(m, D, dDdC, s2.beta, s2.C_env_fun(t_new))
    dlo, ddi, dup = _apply_face_deriv(dlo, ddi, dup, m, D, dDdC, u)
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
    u_new = u + delta
    print(f'it={it:2d}  max|G|={np.max(np.abs(G)):.4e}  max|δ|={np.max(np.abs(delta)):.4e}  '
          f'C表面={u[-1]:.6f}')
    u = u_new
    if np.max(np.abs(delta)) < 1e-12:
        print('收敛')
        break
