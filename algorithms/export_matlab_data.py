# -*- coding: utf-8 -*-
"""
export_matlab_data.py — 导出 MATLAB 绘图所需数据 (.mat)
============================================================
产出 results/matlab/ 下:
  q1.mat        r_cm, t1, T1, C1 (1s×401), 表1/表2 时刻与半径
  q2.mat        t60/T60/C60 (60s×401, 72h), t3/T3/C3 (3h@1s×0.1cm), r_out_cm
  q4.mat        t/C/R (60s, 收缩模型), xi, 附件2 半径数据
  env.mat       附件1 烘房环境 + 平台期统计
  validation.mat 解析解对比 (3 时刻) + 网格收敛数据
"""

import sys, os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
os.chdir(_ROOT)

import numpy as np
import pandas as pd
from scipy.io import savemat
from scipy.interpolate import PchipInterpolator
from validate_solver import series_cyl, q1_solver
from solver_core import CylinderDryer

OUT = 'results/matlab'
os.makedirs(OUT, exist_ok=True)

R = 0.02
N = 400
TABLE_TIMES = [100, 300, 600, 900, 1200, 1500, 1800]
TABLE_R_CM = [0.0, 0.5, 1.0, 1.5, 2.0]


def export_q1():
    d = np.load('results/q1_solution.npz')
    r_cm = d['r'] * 100.0                      # (401,)
    t1 = d['t']                                # (1800,)  t=1..1800 s
    T1 = d['T']; C1 = d['C']                   # (1800, 401)
    savemat(f'{OUT}/q1.mat', {
        'r_cm': r_cm, 't1': t1, 'T1': T1, 'C1': C1,
        'TABLE_TIMES': np.array(TABLE_TIMES, float),
        'TABLE_R_CM': np.array(TABLE_R_CM, float)})
    print('q1.mat 导出完成', T1.shape)


def export_q2():
    d = np.load('results/q2_full.npz')
    t60 = d['t']; T60 = d['T']; C60 = d['C']   # (4320, 401)
    # 3h 内 1s×0.1cm 数据 (从 result2.xlsx 读取, 与提交文件一致)
    xl = pd.ExcelFile('results/result2.xlsx')
    T3 = pd.read_excel(xl, '温度', header=0, index_col=0).values   # (10800, 21)
    C3 = pd.read_excel(xl, '水分浓度', header=0, index_col=0).values
    t3 = pd.read_excel(xl, '温度', header=0, index_col=0).index.values.astype(float)
    r_out = pd.read_excel(xl, '温度', header=0, index_col=0).columns.values.astype(float)
    savemat(f'{OUT}/q2.mat', {
        't60': t60, 'T60': T60, 'C60': C60, 'r_cm': d['r'] * 100.0,
        't3': t3, 'T3': T3, 'C3': C3, 'r_out_cm': r_out})
    print('q2.mat 导出完成', T60.shape, T3.shape)


def export_q4():
    d = np.load('results/q4_full.npz')         # t, xi, C, R (4320, 401)
    df = pd.read_excel('A题/附件/附件2.xlsx')
    savemat(f'{OUT}/q4.mat', {
        't': d['t'], 'xi': d['xi'], 'C': d['C'], 'R': d['R'],
        't_R': df['时间'].values.astype(float), 'R_cm': df['半径'].values.astype(float)})
    print('q4.mat 导出完成', d['C'].shape)


def export_env():
    df = pd.read_excel('A题/附件/附件1.xlsx')
    t = df['时间'].values.astype(float)
    T, C = df['温度'].values, df['水分浓度'].values
    savemat(f'{OUT}/env.mat', {
        't': t, 'T_env': T, 'C_env': C,
        'T_bar': float(T[t >= 7200].mean()), 'T_std': float(T[t >= 7200].std()),
        'C_bar': float(C[t >= 7200].mean()), 'C_std': float(C[t >= 7200].std())})
    print('env.mat 导出完成', len(t))


def export_validation():
    # ① 阶跃边界数值 vs Bessel 级数 (3 时刻, 热场)
    k, rho, cp, h = 0.36, 820.0, 2600.0, 25.0
    alpha = k / (rho * cp)
    s = CylinderDryer(R, N,
        rho_fun=lambda C: np.full_like(C, rho), cp_fun=lambda C: np.full_like(C, cp),
        k_fun=lambda C: np.full_like(C, k),
        D_fun=lambda C, T: np.full_like(C, 5e-9), dDdC_fun=lambda C, T: np.zeros_like(C),
        dDdT_fun=lambda C, T: np.zeros_like(C),
        h=h, beta=8e-7, T_env_fun=lambda t: 50.0, C_env_fun=lambda t: 0.02,
        T0=28.0, C0=2.55)
    t_targets = [100.0, 600.0, 1800.0]
    T_num = []; T_ref = []
    for tt in t_targets:
        while s.t < tt - 1e-12:
            s.step(min(0.25, tt - s.t))
        T_num.append(s.T.copy())
        T_ref.append(series_cyl(s.mesh.r, tt, R, alpha, k, h, 28.0, 50.0))
    # ② 网格收敛 (真实 Q1 参数, t=1800s)
    Ns = [200, 400, 800, 1600]
    err_T, err_C = [], []
    prev_T = prev_C = None
    for Nn in Ns:
        sq = q1_solver(Nn, 0.25)
        sq.run(1800.0, 0.25)
        idx = (np.array([0.0, 0.5, 1.0, 1.5, 2.0]) / 100.0 / (R / Nn)).astype(int)
        if prev_T is not None:
            err_T.append(float(np.max(np.abs(sq.T[idx] - prev_T))))
            err_C.append(float(np.max(np.abs(sq.C[idx] - prev_C))))
        prev_T, prev_C = sq.T[idx].copy(), sq.C[idx].copy()
    savemat(f'{OUT}/validation.mat', {
        'r_cm': s.mesh.r * 100.0, 't_targets': np.array(t_targets),
        'T_num': np.array(T_num), 'T_ref': np.array(T_ref),
        'Ns': np.array(Ns[1:], float), 'err_T': np.array(err_T), 'err_C': np.array(err_C)})
    print('validation.mat 导出完成')


if __name__ == '__main__':
    export_q1()
    export_q2()
    export_q4()
    export_env()
    export_validation()
    print('全部 .mat 导出至', OUT)
