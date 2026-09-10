import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'algorithms'))  # 兄弟模块
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # 项目根
# -*- coding: utf-8 -*-
"""定位 N=800 时 Newton 不收敛的时刻与原因"""
import numpy as np
from q2_solve import build_solver, make_env
from solver_core import CylinderDryer

T_env, C_env, _, _ = make_env()
s = build_solver(T_env, C_env)
s2 = CylinderDryer(0.02, 800,
    rho_fun=s.rho_fun, cp_fun=s.cp_fun, k_fun=s.k_fun,
    D_fun=s.D_fun, dDdC_fun=s.dDdC_fun, dDdT_fun=s.dDdT_fun,
    h=25.0, beta=8e-7, T_env_fun=T_env, C_env_fun=C_env, T0=28.0, C0=2.55)

phases = [(0.25, 300.0), (1.0, 14400.0), (5.0, 259200.0)]
for dt, t_end in phases:
    first = True
    while s2.t < t_end - 1e-12:
        try:
            s2.step(min(dt, t_end - s2.t), force_restart=first)
        except RuntimeError as e:
            print(f'FAIL: t={s2.t:.0f}s, dt={dt}, first={first}, {e}')
            print(f'  C: min={s2.C.min():.6e} max={s2.C.max():.6e}, 表面={s2.C[-1]:.6e}')
            print(f'  T: min={s2.T.min():.4f} max={s2.T.max():.4f}')
            raise SystemExit
        first = False
        if s2.t % 36000 == 0 and s2.t > 0:
            print(f't={s2.t/3600:.0f}h, C表面={s2.C[-1]:.5f}, C中心={s2.C[0]:.5f}')
print('N=800 全流程无失败')
