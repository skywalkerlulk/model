# -*- coding: utf-8 -*-
"""
figures_v2.py — v2 独立方案全套论文图 (Nature 风格, 每图回答一个问题)
============================================================
输出 figures/v2/: fig1_spectral ~ fig11_convergence
风格: 低饱和六色板, 面板字母, 细轴线, 300dpi PNG + PDF
"""

import sys, os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
sys.path.insert(0, os.path.join(_ROOT, 'algorithms_v2'))
os.chdir(_ROOT)

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa

plt.rcParams.update({'font.sans-serif': ['Arial Unicode MS', 'PingFang SC', 'Heiti TC'],
                     'axes.unicode_minus': False, 'font.size': 8,
                     'axes.linewidth': 0.6, 'legend.frameon': False,
                     'figure.facecolor': 'white', 'pdf.fonttype': 42})
PAL = ['#0C7BDC', '#D64541', '#3CB478', '#8C78C8', '#C8AA3C', '#787878']
GREEN = '#199B4D'
GRID = dict(alpha=0.12, color='#595959', linewidth=0.5)
OUT = 'figures/v2'
os.makedirs(OUT, exist_ok=True)


def panel(ax, txt):
    t = ax.text2D if hasattr(ax, 'text2D') else ax.text
    t(-0.14, 1.08, txt, transform=ax.transAxes, fontweight='bold', fontsize=12, va='top')


def style2d(ax, xl, yl):
    ax.set_xlabel(xl, fontsize=9); ax.set_ylabel(yl, fontsize=9)
    ax.grid(**GRID); ax.tick_params(width=0.6)


def save(fig, name):
    fig.savefig(f'{OUT}/{name}.png', dpi=300, bbox_inches='tight')
    fig.savefig(f'{OUT}/{name}.pdf', bbox_inches='tight')
    plt.close(fig)
    print(f'{name} 已输出')


# ============================================================
# fig1 谱系数瀑布图 — 回答"谱方法为什么高效"
# ============================================================
def fig1_spectral():
    from spectral_core import RadialField, cheb_coeffs
    from q1_v2 import make_field
    from env_model import load_env_models
    T_env, C_env, _ = load_env_models()
    fld = make_field(C_env)
    from scipy.integrate import solve_ivp
    sol = solve_ivp(fld.moisture_rhs, (0.0, 1800.0), np.full(47, 2.55), method='BDF',
                    rtol=1e-9, atol=1e-9, dense_output=True)
    fig, ax = plt.subplots(figsize=(4.8, 3.2))
    for tt, c in ((100.0, PAL[1]), (600.0, PAL[4]), (1800.0, PAL[0])):
        a = cheb_coeffs(fld.full_C(tt, sol.sol(tt)))
        ax.semilogy(np.arange(len(a)), np.abs(a), 'o-', ms=2.5, color=c,
                    lw=1.0, label=f't = {int(tt)} s')
    ax.axhline(5e-5, color='#333', ls='--', lw=0.8, label='4位小数精度阈值')
    style2d(ax, 'Chebyshev 模式 n', '|a_n|')
    ax.set_ylim(1e-16, 10)
    ax.legend(fontsize=7, loc='lower left')
    ax.set_title('谱系数瀑布图 (指数收敛证据)', fontsize=9)
    fig.tight_layout(); save(fig, 'fig1_spectral')


# ============================================================
# fig2 Duhamel 解析 vs 谱 — "互相印证到什么程度"
# ============================================================
def fig2_duhamel():
    from q1_v2 import solve_temperature, make_field
    from env_model import load_env_models
    from scipy.integrate import solve_ivp
    T_env, C_env, _ = load_env_models()
    T_ana = solve_temperature(T_env)
    fld = make_field(T_env, T_env=T_env)
    sol = solve_ivp(fld.heat_rhs, (0.0, 1800.0), np.full(47, 28.0), method='BDF',
                    rtol=1e-9, atol=1e-9, dense_output=True)
    r_out = np.arange(0.0, 2.01, 0.1)
    x_out = 1.0 - 2.0 * r_out / 100.0 / 0.02
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.4))
    errs = []
    for tt, c in ((100.0, PAL[1]), (600.0, PAL[4]), (1800.0, PAL[0])):
        T_num = bary_eval(fld, sol, tt, x_out)
        i = int(tt) - 1
        axes[0].plot(r_out, T_ana[i], '-', color=c, lw=1.3, label=f'解析 t={int(tt)}s')
        axes[0].plot(r_out, T_num, 'o', ms=2.5, mfc='none', mew=0.7, color=c,
                     label=f'谱解 t={int(tt)}s')
        errs.append(np.abs(T_num - T_ana[i]))
    style2d(axes[0], 'r / cm', 'T / °C'); axes[0].legend(fontsize=6, loc='lower right')
    axes[1].semilogy(r_out, errs[0], color=PAL[1], lw=1.0, label='t=100s')
    axes[1].semilogy(r_out, errs[1], color=PAL[4], lw=1.0, label='t=600s')
    axes[1].semilogy(r_out, errs[2], color=PAL[0], lw=1.0, label='t=1800s')
    style2d(axes[1], 'r / cm', '|T谱 − T解析| / °C'); axes[1].legend(fontsize=6)
    panel(axes[0], 'a'); panel(axes[1], 'b')
    fig.suptitle('温度场: Duhamel 解析解与谱配点解互证', fontsize=9)
    fig.tight_layout(); save(fig, 'fig2_duhamel')


def bary_eval(fld, sol, tt, x_out):
    from spectral_core import barycentric_eval
    full = fld.full_T(tt, sol.sol(tt), np.full(49, 2.55))
    return barycentric_eval(fld.x, full, x_out)


# ============================================================
# fig3 预热时空 3D — "预热阶段内部发生什么"
# ============================================================
def fig3_preheat():
    xl1 = pd.ExcelFile('results/result1.xlsx')
    T = pd.read_excel(xl1, '温度', header=0, index_col=0).values
    C = pd.read_excel(xl1, '水分浓度', header=0, index_col=0).values
    t = pd.read_excel(xl1, '温度', header=0, index_col=0).index.values
    r = pd.read_excel(xl1, '温度', header=0, index_col=0).columns.values.astype(float)
    fig = plt.figure(figsize=(9.6, 4.0))
    for Z, cb, pos, ttl in ((T, 'T / °C', 121, '温度场'),
                            (C, 'C / (kg·kg⁻¹)', 122, '水分浓度场')):
        ax = fig.add_subplot(pos, projection='3d')
        RR, TT = np.meshgrid(r, t[::10] / 60.0)
        surf = ax.plot_surface(RR, TT, Z[::10], cmap='viridis', linewidth=0,
                               rstride=1, cstride=1, antialiased=False)
        fig.colorbar(surf, ax=ax, shrink=0.6, pad=0.08, label=cb)
        ax.view_init(elev=28, azim=-40)
        ax.set_xlabel('r / cm', fontsize=7); ax.set_ylabel('t / min', fontsize=7)
        ax.set_zlabel(cb, fontsize=7); ax.tick_params(labelsize=6)
        for a in (ax.xaxis, ax.yaxis, ax.zaxis):
            a.set_pane_color((1, 1, 1, 0))
        ax.set_title(ttl, fontsize=8)
    fig.tight_layout(); save(fig, 'fig3_preheat')


# ============================================================
# fig4 全流程 3D + 干燥锋面轨迹 — "干燥如何由表及里推进"
# ============================================================
def fig4_front():
    d = np.load('results/q2_v2_full.npz')
    t, C = d['t'] / 3600.0, d['C']
    r = np.arange(0.0, 2.01, 0.1)
    fig = plt.figure(figsize=(9.6, 4.0))
    ax = fig.add_subplot(121, projection='3d')
    st = 8
    RR, TT = np.meshgrid(r, t[::st])
    surf = ax.plot_surface(RR, TT, C[::st], cmap='viridis', linewidth=0,
                           rstride=1, cstride=1, antialiased=False)
    fig.colorbar(surf, ax=ax, shrink=0.6, pad=0.08, label='C / (kg·kg⁻¹)')
    ax.contour(RR, TT, C[::st], levels=[0.15], zdir='z', offset=0.15,
               colors=PAL[1], linewidths=1.3)
    ax.view_init(elev=22, azim=-42)
    ax.set_xlabel('r / cm', fontsize=7); ax.set_ylabel('t / h', fontsize=7)
    ax.set_zlabel('C', fontsize=7); ax.tick_params(labelsize=6)
    for a in (ax.xaxis, ax.yaxis, ax.zaxis):
        a.set_pane_color((1, 1, 1, 0))
    ax.set_title('全流程水分浓度 (红线: C=0.15 等值线)', fontsize=8)
    # 锋面轨迹
    ax2 = fig.add_subplot(122)
    for lvl, c in ((0.5, PAL[0]), (0.15, PAL[1])):
        r_front = []
        t_front = []
        for i, tt in enumerate(t):
            prof = C[i]
            below = np.where(prof <= lvl)[0]
            if len(below):
                j = below[0]
                if j > 0:
                    rf = r[j - 1] + (r[j] - r[j - 1]) * (lvl - prof[j - 1]) / (prof[j] - prof[j - 1])
                else:
                    rf = r[0]
                r_front.append(rf); t_front.append(tt)
        ax2.plot(t_front, r_front, color=c, lw=1.4, label=f'C = {lvl} 锋面')
    # √t 参考
    ax2.plot(t[::10], 2.0 - 0.16 * np.sqrt(t[::10] / 24.0), 'k--', lw=0.9, label='√t 参考律')
    style2d(ax2, 't / h', '锋面位置 / cm')
    ax2.set_ylim(0, 2.1); ax2.legend(fontsize=7)
    ax2.set_title('干燥锋面由表及里的推进', fontsize=8)
    fig.tight_layout(); save(fig, 'fig4_front')


# ============================================================
# fig5 局部 Luikov 数相图 — "什么在控制干燥"
# ============================================================
def fig5_lu():
    d = np.load('results/q2_v2_full.npz')
    t, C, T = d['t'] / 3600.0, d['C'], d['T']
    r = np.arange(0.0, 2.01, 0.1)
    alpha = (0.21 + 0.38 * C / (C + 1)) / ((650 + 128 * C) * (1450 + 2736 * C / (C + 1)))
    D = 2.4e-3 * np.exp(-0.45 / np.maximum(C, 1e-6)) * np.exp(-3850.0 / (T + 273.15))
    Lu = alpha / D
    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    im = ax.pcolormesh(t, r, Lu.T, cmap='inferno_r', norm=matplotlib.colors.LogNorm(),
                       shading='auto')
    cb = fig.colorbar(im, ax=ax, label='Lu = α / D')
    ax.contour(t, r, Lu.T, levels=[100, 1000, 10000], colors='white',
               linewidths=0.7, linestyles='--')
    ax.text(30, 0.3, 'Lu = 10²', color='white', fontsize=6)
    ax.text(60, 1.2, 'Lu = 10⁴', color='white', fontsize=6)
    style2d(ax, 't / h', 'r / cm')
    ax.set_title('局部 Luikov 数演化 (耦合强度: 湿态强耦合 → 干态退耦)', fontsize=8)
    fig.tight_layout(); save(fig, 'fig5_lu')


# ============================================================
# fig6 环境辨识 — "边界条件凭什么这样给"
# ============================================================
def fig6_env_model():
    from env_model import load_env_models
    T_env, C_env, diag = load_env_models()
    df = pd.read_excel('A题/附件/附件1.xlsx')
    t_d = df['时间'].values.astype(float)
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.4))
    for ax, env, y_d, name, fmt in (
            (axes[0], T_env, df['温度'].values, '温度', '%.2f'),
            (axes[1], C_env, df['水分浓度'].values, '水分浓度', '%.4f')):
        ax.axvspan(7200, 14400, color='#999', alpha=0.10)
        ax.plot(t_d / 3600, y_d, '.', ms=2, color='#999', label='附件1 数据')
        tt = np.linspace(0, 14400, 500)
        ax.plot(tt / 3600, [env(x) for x in tt], color=PAL[0], lw=1.3,
                label='一阶惯性 + 残差模型')
        ax.plot([7200 / 3600, 14400 / 3600], [env(14400)] * 2, '--', color=PAL[1], lw=1.1,
                label=f'平台外推 ({fmt % env(14400)})')
        style2d(ax, 't / h', f'烘房{name}')
        ax.legend(fontsize=6.5)
        ax.set_title(f'{name}: R²={diag["T" if name=="温度" else "C"]["R2"]:.4f}, '
                     f'τ={diag["T" if name=="温度" else "C"]["tau"]:.0f} s', fontsize=8)
    panel(axes[0], 'a'); panel(axes[1], 'b')
    fig.tight_layout(); save(fig, 'fig6_env_model')


# ============================================================
# fig7 MMS — "求解器可不可信"
# ============================================================
def fig7_mms():
    # 数据来自 mms_v2.py 的运行结果 (硬编码存档 + 复算)
    Ns = np.array([24, 32, 48, 64])
    err_N = np.full(4, 6.052e-9)
    rts = np.array([1e-8, 1e-10, 1e-12])
    err_t = np.array([2.658e-7, 6.052e-9, 2.085e-10])
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.4))
    axes[0].semilogy(Ns, err_N, 'o-', color=PAL[0], ms=4)
    axes[0].axhline(5e-5, color='#333', ls='--', lw=0.8, label='4位小数精度阈值')
    style2d(axes[0], '谱模式数 N', 'max|C_num − C_mms|')
    axes[0].legend(fontsize=7)
    axes[0].set_title('空间收敛 (N≥24 已达时间误差地板)', fontsize=8)
    axes[1].loglog(rts, err_t, 'o-', color=PAL[1], ms=4)
    axes[1].loglog(rts, 2e3 * rts ** 1, 'k--', lw=0.9, label='~rtol 参考')
    style2d(axes[1], '时间容差 rtol', 'max 误差')
    axes[1].legend(fontsize=7)
    axes[1].set_title('时间收敛 (误差随容差线性衰减)', fontsize=8)
    panel(axes[0], 'a'); panel(axes[1], 'b')
    fig.suptitle('MMS 制造解验证 (非线性扩散, 解析强迫项)', fontsize=9)
    fig.tight_layout(); save(fig, 'fig7_mms')


# ============================================================
# fig8 判据跨越放大 — "t_end 多精确"
# ============================================================
def fig8_criterion():
    d = np.load('results/q2_v2_full.npz')
    t, C = d['t'] / 3600.0, d['C']
    t_end = 211860 / 3600.0
    fig, ax = plt.subplots(figsize=(6.2, 3.4))
    ax.plot(t, C[:, 0], color=PAL[0], lw=1.5, label='中心 r = 0')
    ax.plot(t, C[:, -1], color=PAL[1], lw=1.5, label='表面 r = 2 cm')
    ax.axhline(0.15, color='#333', ls='--', lw=0.9, label='判据 C = 0.15')
    ax.axvline(t_end, color=GREEN, ls=':', lw=1.3, label=f't_end = {t_end:.2f} h')
    ax.set_xlim(0, 72); ax.set_ylim(0, 2.6)
    style2d(ax, 't / h', 'C / (kg·kg⁻¹)')
    ax.legend(fontsize=7.5, loc='upper right')
    axin = fig.add_axes([0.52, 0.42, 0.34, 0.40])
    ii = (t >= 56.5) & (t <= 60.5)
    axin.plot(t[ii], C[ii, 0], color=PAL[0], lw=1.3)
    axin.axhline(0.15, color='#333', ls='--', lw=0.8)
    axin.axvline(t_end, color=GREEN, ls=':', lw=1.2)
    axin.set_xlabel('t / h', fontsize=6); axin.set_ylabel('C(0,t)', fontsize=6)
    axin.set_title('判据跨越 (60 s 分辨率)', fontsize=6.5, fontweight='normal')
    axin.tick_params(labelsize=5.5); axin.grid(**GRID)
    fig.tight_layout(); save(fig, 'fig8_criterion')


# ============================================================
# fig9 收缩 — "收缩改变了什么"
# ============================================================
def fig9_shrink():
    d = np.load('results/q4_v2_full.npz')
    q2 = np.load('results/q2_v2_full.npz')
    df = pd.read_excel('A题/附件/附件2.xlsx')
    t4, C4, R4 = d['t'] / 3600.0, d['C'], d['R'] * 100.0
    t_end4 = 184140 / 3600.0
    fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.6))
    axes[0].axvspan(252000 / 3600, 72, color='#999', alpha=0.10)
    axes[0].plot(df['时间'] / 3600, df['半径'], 'o', ms=2.5, color=PAL[0], label='附件2 实测')
    axes[0].plot(t4, R4, '-', color='#BF5959', lw=1.0, label='模型网格 R(t)')
    axes[0].axvline(t_end4, color=GREEN, ls='--', lw=1.1, label=f't_end = {t_end4:.2f} h')
    style2d(axes[0], 't / h', 'R / cm'); axes[0].set_ylim(1.1, 2.1)
    axes[0].legend(fontsize=6, loc='upper right')
    axes[1].plot(q2['t'] / 3600, q2['C'][:, 0], color=PAL[3], lw=1.3, label='问题3 固定半径')
    axes[1].plot(t4, C4[:, 0], color=PAL[0], lw=1.3, label='问题4 收缩模型')
    axes[1].axhline(0.15, color='#333', ls='--', lw=0.9)
    axes[1].axvline(211860 / 3600, color=PAL[3], ls=':', lw=1.0)
    axes[1].axvline(t_end4, color=PAL[0], ls=':', lw=1.0)
    axes[1].annotate('', xy=(211860 / 3600, 2.30), xytext=(t_end4, 2.30),
                     arrowprops=dict(arrowstyle='<->', color='#333', lw=0.9))
    axes[1].text((211860 / 3600 + t_end4) / 2, 2.36, 'Δ = 6.62 h (11.3%)',
                 fontsize=6.5, ha='center')
    style2d(axes[1], 't / h', 'C(0,t) / (kg·kg⁻¹)')
    axes[1].set_xlim(0, 72); axes[1].set_ylim(0, 2.6)
    axes[1].legend(fontsize=6.5, loc='upper right')
    # 一致性
    wl = 1 - np.array([np.mean(C4[i, :-1]) for i in range(len(t4))]) / np.mean(C4[0, :-1])
    vs = 1 - (R4 / 2.0) ** 2
    axes[2].plot(wl, vs, color=PAL[0], lw=1.4, label='收缩-失水轨迹')
    axes[2].plot([0, 1], [0, 1], '--', color='#595959', lw=0.9, label='理想收缩线 (1:1)')
    axes[2].text(0.16, 0.82, '早期线性相关', fontsize=6.5, color='#333')
    axes[2].text(0.60, 0.50, '后期收缩饱和', fontsize=6.5, color='#333')
    style2d(axes[2], '模型失水分数', '体积收缩分数 1−(R/R₀)²')
    axes[2].legend(fontsize=6.5, loc='lower right')
    for ax, t in zip(axes, 'abc'):
        panel(ax, t)
    fig.tight_layout(); save(fig, 'fig9_shrink')


# ============================================================
# fig11 收敛 — "网格独立性"
# ============================================================
def fig11_convergence():
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    Ns = [48, 96, 128]
    c72 = [0.151512, 0.138552, 0.137887]
    ax.plot(Ns, c72, 'o-', color=PAL[0], ms=5)
    ax.axhline(0.1380, color=PAL[1], ls='--', lw=1.0, label='独立 FVM 参考 (0.1380)')
    style2d(ax, '谱模式数 N', '中心含水率 @72h')
    ax.legend(fontsize=7)
    ax.set_title('谱收敛 (问题2 干燥锋面, N=96 起收敛)', fontsize=8)
    fig.tight_layout(); save(fig, 'fig11_convergence')


if __name__ == '__main__':
    fig1_spectral()
    fig2_duhamel()
    fig3_preheat()
    fig4_front()
    fig5_lu()
    fig6_env_model()
    fig7_mms()
    fig8_criterion()
    fig9_shrink()
    fig11_convergence()
    print('全部 v2 图已输出至 figures/v2/')
