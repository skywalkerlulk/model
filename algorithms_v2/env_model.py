# -*- coding: utf-8 -*-
"""
env_model.py — 烘房环境的一阶惯性模型 (2026 CUMCM A题 v2 独立方案)
==================================================================
物理依据: 烘房升温由加热系统惯性决定, 温度/湿度呈指数松弛逼近平台
    T_env(t) = T_∞ − (T_∞ − T_0)·e^(−t/τ)        (一阶惯性环节)
附件1 拟合: R²_T = 0.9973 (τ=1776 s), R²_C = 0.9944 (τ=2522 s)
残差 (数据 − 拟合) 由 PCHIP 光滑并叠加, 保证完全忠实于数据;
t > 4 h 平台外推 (数据统计: 49.97±0.17 °C / 0.0500 kg/kg)。
用法: T_env(t) 可调用; .deriv(t) 供 Duhamel 卷积; .tau/.T_inf 供解析主项。
"""

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.interpolate import PchipInterpolator

T_PLATEAU = 7200.0        # 平台期起点 (s, 数据统计)
T_DATA_END = 14400.0      # 附件1 数据终点 (s)


def _exp_relax(t, y_inf, y0, tau):
    return y_inf - (y_inf - y0) * np.exp(-t / tau)


class EnvModel:
    """单个环境量 (温度或水分浓度) 的一阶惯性 + 残差模型"""

    def __init__(self, t, y):
        idx = t <= T_DATA_END
        t_fit, y_fit = t[idx], y[idx]
        popt, _ = curve_fit(_exp_relax, t_fit, y_fit,
                            p0=(y_fit[-1], y_fit[0], 2000.0), maxfev=40000)
        self.y_inf, self.y0, self.tau = popt
        fit = _exp_relax(t_fit, *popt)
        resid = y_fit - fit
        self.r2 = 1.0 - np.sum(resid ** 2) / np.sum((y_fit - y_fit.mean()) ** 2)
        # 平台期统计 (外推依据)
        plat = y[t >= T_PLATEAU]
        self.plateau_mean = float(plat.mean())
        self.plateau_std = float(plat.std())
        self._resid = PchipInterpolator(t_fit, resid)   # 残差光滑 (单调保形)
        self._t_max = t_fit[-1]

    def __call__(self, t):
        """T_env(t): 拟合主项 + 残差; t>4h 平台外推"""
        t = np.asarray(t, dtype=float)
        main = np.where(t <= self._t_max, _exp_relax(t, self.y_inf, self.y0, self.tau),
                        self.plateau_mean)
        res = np.zeros_like(t)
        m = t <= self._t_max
        res[m] = self._resid(t[m])
        return main + res

    def deriv(self, t):
        """dT_env/dt: 主项解析导数 + 残差 PCHIP 导数"""
        t = np.asarray(t, dtype=float)
        d = np.zeros_like(t)
        m = t <= self._t_max
        d[m] = (self.y_inf - self.y0) / self.tau * np.exp(-t[m] / self.tau)
        d[m] += self._resid.derivative()(t[m])
        return d


def load_env_models():
    """返回 (T_env, C_env, 拟合诊断字典)"""
    df = pd.read_excel('A题/附件/附件1.xlsx')
    t = df['时间'].values.astype(float)
    T_model = EnvModel(t, df['温度'].values)
    C_model = EnvModel(t, df['水分浓度'].values)
    diag = {
        'T': dict(R2=T_model.r2, tau=T_model.tau, T_inf=T_model.y_inf, T0=T_model.y0,
                  plateau=(T_model.plateau_mean, T_model.plateau_std)),
        'C': dict(R2=C_model.r2, tau=C_model.tau, C_inf=C_model.y_inf, C0=C_model.y0,
                  plateau=(C_model.plateau_mean, C_model.plateau_std)),
    }
    return T_model, C_model, diag


if __name__ == '__main__':
    Tm, Cm, diag = load_env_models()
    for k, v in diag.items():
        print(f'{k}: R2={v["R2"]:.6f}, tau={v["tau"]:.1f}s, '
              f'平台均值±std = {v["plateau"][0]:.4f} ± {v["plateau"][1]:.4f}')
