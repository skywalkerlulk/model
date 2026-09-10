# -*- coding: utf-8 -*-
"""
preview_nature_figs.py — Nature 风格图预览 v2 (与 figures_matlab/*.m 同构)
============================================================
v2 升级: 摩尔纹修正 / 0.15 等值线+结束平面 / 多温度 D-C 曲线 / Δt 标注 /
         收缩平台阴影 / ±σ 波动带 / 守恒残差标注 / 环境驱动线
本机无 MATLAB 时生成预览; 正式出图用 figures_matlab/run_all_figs.m
输出: figures/preview/fig1_q1 ~ fig6_env (300 dpi PNG + 矢量 PDF)
"""

import sys, os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
os.chdir(_ROOT)

import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.sans-serif': ['Arial Unicode MS', 'PingFang SC', 'Heiti TC', 'STHeiti'],
    'axes.unicode_minus': False,
    'font.size': 8, 'axes.linewidth': 0.6,
    'axes.titlesize': 8, 'xtick.direction': 'out', 'ytick.direction': 'out',
    'legend.frameon': False, 'figure.facecolor': 'white',
    'savefig.facecolor': 'white', 'pdf.fonttype': 42, 'ps.fonttype': 42,
})

PAL = ['#0C7BDC', '#D64541', '#3CB478', '#8C78C8', '#C8AA3C', '#787878']
GRID = dict(alpha=0.12, color='#595959', linewidth=0.5)
GREEN = '#199B4D'
OUT = 'figures/preview'
os.makedirs(OUT, exist_ok=True)

Q1 = loadmat('results/matlab/q1.mat')
Q2 = loadmat('results/matlab/q2.mat')
Q4 = loadmat('results/matlab/q4.mat')
ENV = loadmat('results/matlab/env.mat')
VAL = loadmat('results/matlab/validation.mat')

T_END3 = 207960 / 3600.0     # 问题3 烘干结束 (h)
T_END4 = 184200 / 3600.0     # 问题4 烘干结束 (h)


def panel(ax, txt):
    t = ax.text2D if hasattr(ax, 'text2D') else ax.text
    t(-0.13, 1.06, txt, transform=ax.transAxes, fontweight='bold', fontsize=12, va='top')


def style2d(ax, xl, yl):
    ax.set_xlabel(xl, fontsize=9); ax.set_ylabel(yl, fontsize=9)
    ax.grid(**GRID); ax.tick_params(width=0.6)


def style3d(ax, xl, yl, zl, elev=26, azim=-38):
    ax.set_xlabel(xl, fontsize=8, labelpad=2)
    ax.set_ylabel(yl, fontsize=8, labelpad=2)
    ax.set_zlabel(zl, fontsize=8, labelpad=2)
    ax.view_init(elev=elev, azim=azim)
    ax.xaxis.set_pane_color((1, 1, 1, 0)); ax.yaxis.set_pane_color((1, 1, 1, 0))
    ax.zaxis.set_pane_color((1, 1, 1, 0))
    ax.tick_params(labelsize=6.5, pad=0.5)


def save(fig, name):
    fig.savefig(f'{OUT}/{name}.png', dpi=300, bbox_inches='tight')
    fig.savefig(f'{OUT}/{name}.pdf', bbox_inches='tight')
    plt.close(fig)
    print(f'{name} 已输出')


# ============================================================
# fig1_q1 — 问题1: 三维网格(带环境驱动线) + 剖面
# ============================================================
def fig1_q1():
    r, t, T1, C1 = Q1['r_cm'].ravel(), Q1['t1'].ravel(), Q1['T1'], Q1['C1']
    tt = Q1['TABLE_TIMES'].ravel()
    t_env = ENV['t'].ravel(); T_env = ENV['T_env'].ravel()
    fig = plt.figure(figsize=(9.8, 6.8))
    ax_a = fig.add_subplot(2, 2, 1, projection='3d')
    ax_b = fig.add_subplot(2, 2, 2, projection='3d')
    ax_c = fig.add_subplot(2, 2, 3)
    ax_d = fig.add_subplot(2, 2, 4)
    st = 5
    t_m = t[::st] / 60.0
    R, TT = np.meshgrid(r, t_m)
    for ax, Z, cb_lab in ((ax_a, T1[::st], 'T / °C'),
                          (ax_b, C1[::st], 'C / (kg·kg⁻¹)')):
        surf = ax.plot_surface(R, TT, Z, cmap='viridis', linewidth=0.08,
                               edgecolors=(0.08, 0.10, 0.16, 0.22),
                               rstride=10, cstride=2, antialiased=True)
        fig.colorbar(surf, ax=ax, shrink=0.62, pad=0.10, label=cb_lab)
    # 面板 a: 环境温度驱动线 (r = 2 cm 处的边界条件)
    T_env_i = np.interp(t[::st], t_env, T_env)
    ax_a.plot([2.0] * len(t_m), t_m, T_env_i, color=PAL[1], lw=1.4,
              label='烘房温度 T_env(t)')
    ax_a.legend(fontsize=6, loc='upper left')
    style3d(ax_a, 'r / cm', 't / min', 'T / °C')
    style3d(ax_b, 'r / cm', 't / min', 'C / (kg·kg⁻¹)')
    panel(ax_a, 'a'); panel(ax_b, 'b')
    for k, ts in enumerate(tt):
        i = int(ts) - 1
        ax_c.plot(r, T1[i], color=PAL[k % 6], lw=1.1, label=f'{int(ts)} s')
        ax_d.plot(r, C1[i], color=PAL[k % 6], lw=1.1, label=f'{int(ts)} s')
    style2d(ax_c, 'r / cm', 'T / °C'); style2d(ax_d, 'r / cm', 'C / (kg·kg⁻¹)')
    ax_c.legend(fontsize=6.5, loc='lower right')
    ax_d.legend(fontsize=6.5, loc='upper right')
    panel(ax_c, 'c'); panel(ax_d, 'd')
    fig.tight_layout()
    save(fig, 'fig1_q1')


# ============================================================
# fig2_q2 — 问题2: 剖面 + 历程(带结束标记) + 多温度 D-C 机理
# ============================================================
def fig2_q2():
    t60, C60, t3, T3, C3, r_out = (Q2['t60'].ravel(), Q2['C60'], Q2['t3'].ravel(),
                                   Q2['T3'], Q2['C3'], Q2['r_out_cm'].ravel())
    fig, axes = plt.subplots(2, 2, figsize=(9.8, 6.8))
    iT = np.array([int(t) - 1 for t in range(1800, 10801, 1800)])
    lab = [f'{t/3600:.1f} h' for t in range(1800, 10801, 1800)]
    for k, i in enumerate(iT):
        axes[0, 0].plot(r_out, T3[i], color=PAL[k % 6], lw=1.1, label=lab[k])
        axes[0, 1].plot(r_out, C3[i], color=PAL[k % 6], lw=1.1, label=lab[k])
    style2d(axes[0, 0], 'r / cm', 'T / °C'); style2d(axes[0, 1], 'r / cm', 'C / (kg·kg⁻¹)')
    axes[0, 0].legend(fontsize=6.5, loc='lower right')
    axes[0, 1].legend(fontsize=6.5, loc='upper right')
    axes[1, 0].plot(t60 / 3600, C60[:, 0], color=PAL[0], lw=1.4, label='中心 r = 0')
    axes[1, 0].plot(t60 / 3600, C60[:, -1], color=PAL[1], lw=1.4, label='表面 r = 2 cm')
    axes[1, 0].axhline(0.15, color='#333', ls='--', lw=0.8, label='判据 0.15')
    axes[1, 0].axvline(T_END3, color=GREEN, ls=':', lw=1.1, label=f't_end = {T_END3:.2f} h')
    style2d(axes[1, 0], 't / h', 'C / (kg·kg⁻¹)')
    axes[1, 0].set_xlim(0, 72)
    axes[1, 0].legend(fontsize=7, loc='upper right')
    # 多温度 D-C 曲线 (Arrhenius 机理)
    Cg = np.linspace(0.05, 2.55, 300)
    for Tc, c in ((28.0, PAL[4]), (40.0, PAL[3]), (50.0, PAL[0])):
        Dg = 2.4e-3 * np.exp(-0.45 / Cg) * np.exp(-3850 / (Tc + 273.15))
        axes[1, 1].semilogy(Cg, Dg, color=c, lw=1.3, label=f'T = {Tc:.0f} °C')
    axes[1, 1].axvline(0.15, color='#333', ls='--', lw=0.8)
    axes[1, 1].axvline(2.55, color='#595959', ls=':', lw=0.8)
    axes[1, 1].text(0.30, 3e-8, '判据 C = 0.15', fontsize=6.5, color='#333', rotation=90)
    axes[1, 1].text(2.28, 3e-8, '湿态冻结点', fontsize=6.5, color='#595959')
    style2d(axes[1, 1], 'C / (kg·kg⁻¹)', 'D / (m²·s⁻¹)')
    axes[1, 1].legend(fontsize=6.5, loc='lower left')
    for ax, t in zip(axes.flat, 'abcd'):
        panel(ax, t)
    fig.tight_layout()
    save(fig, 'fig2_q2')


# ============================================================
# fig3_q3 — 问题3: 三维曲面(0.15等值线+结束平面+等高线投影) + 判据跨越
# ============================================================
def fig3_q3():
    t60, C60, r = Q2['t60'].ravel(), Q2['C60'], Q2['r_cm'].ravel()
    fig = plt.figure(figsize=(9.8, 4.3))
    ax = fig.add_subplot(121, projection='3d')
    st_t, st_r = 4, 2
    r4, t4, Z4 = r[::st_r], t60[::st_t] / 3600.0, C60[::st_t][:, ::st_r]
    R, TT = np.meshgrid(r4, t4)
    surf = ax.plot_surface(R, TT, Z4, cmap='viridis', linewidth=0, rstride=4,
                           cstride=1, antialiased=False)
    ax.contourf(R, TT, Z4, zdir='z', offset=0.0, cmap='viridis', alpha=0.55, levels=12)
    # 0.15 等值线 (落在曲面上)
    ax.contour(R, TT, Z4, levels=[0.15], zdir='z', offset=0.15, colors=PAL[1],
               linewidths=1.2, linestyles='-')
    # 烘干结束时刻半透明平面
    RRp, ZZp = np.meshgrid(r4, np.linspace(0, Z4.max(), 8))
    TTp = np.full_like(ZZp, T_END3)
    ax.plot_surface(RRp, TTp, ZZp, color=GREEN, alpha=0.16, linewidth=0)
    fig.colorbar(surf, ax=ax, shrink=0.62, pad=0.10, label='C / (kg·kg⁻¹)')
    style3d(ax, 'r / cm', 't / h', 'C / (kg·kg⁻¹)', elev=22, azim=-42)
    panel(ax, 'a')
    ax2 = fig.add_subplot(122)
    ax2.plot(t60 / 3600, C60[:, 0], color=PAL[0], lw=1.4, label='中心 r = 0')
    ax2.plot(t60 / 3600, C60[:, -1], color=PAL[1], lw=1.4, label='表面 r = 2 cm')
    ax2.axhline(0.15, color='#333', ls='--', lw=0.9, label='判据 C = 0.15')
    ax2.axvline(T_END3, color=GREEN, ls=':', lw=1.2, label=f't_end = {T_END3:.2f} h')
    style2d(ax2, 't / h', 'C / (kg·kg⁻¹)')
    ax2.set_xlim(0, 72); ax2.set_ylim(0, 2.6)
    ax2.legend(fontsize=7, loc='upper right')
    axin = fig.add_axes([0.668, 0.40, 0.19, 0.42])
    ii = (t60 >= 55 * 3600) & (t60 <= 59 * 3600)
    axin.plot(t60[ii] / 3600, C60[ii, 0], color=PAL[0], lw=1.2)
    axin.axhline(0.15, color='#333', ls='--', lw=0.7)
    axin.axvline(T_END3, color=GREEN, ls=':', lw=1.1)
    axin.set_xlabel('t / h', fontsize=6); axin.set_ylabel('C(0,t)', fontsize=6)
    axin.set_title('判据跨越 (60 s 分辨率)', fontsize=6.5, fontweight='normal')
    axin.tick_params(labelsize=5.5); axin.grid(**GRID)
    panel(ax2, 'b')
    fig.tight_layout()
    save(fig, 'fig3_q3')


# ============================================================
# fig4_q4 — 问题4: 收缩(平台阴影) + 对比(Δt标注) + 一致性
# ============================================================
def fig4_q4():
    t, C, R, t_R, R_cm = (Q4['t'].ravel(), Q4['C'], Q4['R'].ravel(),
                          Q4['t_R'].ravel(), Q4['R_cm'].ravel())
    fig, axes = plt.subplots(1, 3, figsize=(10.6, 3.7))
    # a) 半径收缩 + 平台阴影
    axes[0].axvspan(252000 / 3600, 72, color='#999999', alpha=0.12)
    axes[0].plot(t_R / 3600, R_cm, 'o', ms=2.5, color=PAL[0], label='附件2 实测')
    axes[0].plot(t / 3600, R * 100, '-', color='#BF5959', lw=1.0, label='模型 PCHIP 插值')
    axes[0].axvline(70, color='#333', ls=':', lw=0.9, label='收缩平台 70 h')
    axes[0].axvline(T_END4, color=GREEN, ls='--', lw=1.1, label=f't_end = {T_END4:.2f} h')
    axes[0].text(70.4, 1.16, '收缩平台', fontsize=6, color='#595959')
    style2d(axes[0], 't / h', 'R / cm'); axes[0].set_ylim(1.1, 2.1)
    axes[0].legend(fontsize=6, loc='upper right')
    # b) 问题3 vs 问题4 + Δt 双箭头标注
    axes[1].plot(Q2['t60'].ravel() / 3600, Q2['C60'][:, 0], color=PAL[3], lw=1.3, label='问题3 固定半径')
    axes[1].plot(t / 3600, C[:, 0], color=PAL[0], lw=1.3, label='问题4 收缩模型')
    axes[1].axhline(0.15, color='#333', ls='--', lw=0.9, label='判据 0.15')
    axes[1].axvline(T_END3, color=PAL[3], ls=':', lw=1.0)
    axes[1].axvline(T_END4, color=PAL[0], ls=':', lw=1.0)
    axes[1].annotate('', xy=(T_END3, 2.30), xytext=(T_END4, 2.30),
                     arrowprops=dict(arrowstyle='<->', color='#333', lw=0.9))
    axes[1].text((T_END3 + T_END4) / 2, 2.36, f'Δ = {T_END3 - T_END4:.2f} h (11.4%)',
                 fontsize=6.5, ha='center', color='#333')
    style2d(axes[1], 't / h', 'C(0,t) / (kg·kg⁻¹)')
    axes[1].set_xlim(0, 72); axes[1].set_ylim(0, 2.6)
    axes[1].legend(fontsize=6.5, loc='upper right')
    # c) 一致性
    xi = Q4['xi'].ravel(); dxi = 1 / 400
    xif = (np.arange(400) + 0.5) * dxi
    w = np.zeros(401)
    w[0] = xif[0] ** 2
    w[1:400] = xif[1:] ** 2 - xif[:-1] ** 2
    w[400] = 1 - xif[-1] ** 2
    Mw = C * w
    wl = 1 - Mw.sum(axis=1) / Mw[0].sum()
    vs = 1 - (R / 0.02) ** 2
    axes[2].plot(wl, vs, color=PAL[0], lw=1.4, label='收缩-失水轨迹')
    axes[2].plot([0, 1], [0, 1], '--', color='#595959', lw=0.9, label='理想收缩线 (1:1)')
    axes[2].annotate('早期线性相关\n后期收缩饱和', xy=(0.62, 0.62), xytext=(0.18, 0.82),
                     fontsize=6.5, color='#333',
                     arrowprops=dict(arrowstyle='->', color='#333', lw=0.7))
    style2d(axes[2], '模型失水分数 1 − M_w/M_w0', '体积收缩分数 1 − (R/R₀)²')
    axes[2].legend(fontsize=6.5, loc='lower right')
    for ax, t in zip(axes, 'abc'):
        panel(ax, t)
    fig.tight_layout()
    save(fig, 'fig4_q4')


# ============================================================
# fig5_validation — 验证 (解析解 + 收敛 + 守恒残差)
# ============================================================
def fig5_validation():
    r, tt, Tn, Tr, Ns, eT, eC = (VAL['r_cm'].ravel(), VAL['t_targets'].ravel(), VAL['T_num'],
                                 VAL['T_ref'], VAL['Ns'].ravel(), VAL['err_T'].ravel(),
                                 VAL['err_C'].ravel())
    fig = plt.figure(figsize=(9.8, 4.0))
    ax = fig.add_subplot(121)
    for k, ts in enumerate(tt):
        ax.plot(r, Tn[k], 'o', ms=3, color=PAL[k], mfc='none', mew=0.8,
                label=f'数值 t={int(ts)} s')
        ax.plot(r, Tr[k], '-', color=PAL[k], lw=1.2, label=f'解析 t={int(ts)} s')
    style2d(ax, 'r / cm', 'T / °C'); ax.legend(fontsize=6, loc='lower right')
    axin = fig.add_axes([0.155, 0.42, 0.20, 0.38])
    for k in range(3):
        axin.semilogy(r, np.abs(Tn[k] - Tr[k]), color=PAL[k], lw=0.9)
    axin.set_xlabel('r / cm', fontsize=6); axin.set_ylabel('|误差| / °C', fontsize=6)
    axin.set_title('误差分布 (max 2.5×10⁻⁶)', fontsize=6.5, fontweight='normal')
    axin.tick_params(labelsize=5.5); axin.grid(**GRID)
    panel(ax, 'a')
    ax2 = fig.add_subplot(122)
    ax2.loglog(Ns, eT, 'o-', color=PAL[0], ms=4, label='温度')
    ax2.loglog(Ns, eC, 's-', color=PAL[1], ms=4, label='水分浓度')
    ax2.loglog(Ns, 1e-4 * (400 / Ns) ** 2, 'k--', lw=0.9, label='二阶参考线')
    style2d(ax2, '网格数 N', '相邻网格最大差')
    ax2.legend(fontsize=6.5, loc='lower left')
    ax2.text(0.97, 0.14,
             '守恒残差 (Q1):\n水分 7.6×10⁻⁷\n能量 7.8×10⁻⁸\n时间收敛: dt减半差 < 10⁻⁷',
             transform=ax2.transAxes, fontsize=6.5, va='bottom', ha='right',
             bbox=dict(boxstyle='round,pad=0.35', fc='white', ec='#999999', lw=0.4))
    panel(ax2, 'b')
    fig.tight_layout()
    save(fig, 'fig5_validation')


# ============================================================
# fig6_env — 环境条件: 平台阴影 + ±σ 波动带
# ============================================================
def fig6_env():
    t, Te, Ce = ENV['t'].ravel(), ENV['T_env'].ravel(), ENV['C_env'].ravel()
    Tb = float(np.asarray(ENV['T_bar']).ravel()[0])
    Cb = float(np.asarray(ENV['C_bar']).ravel()[0])
    Ts = float(np.asarray(ENV['T_std']).ravel()[0])
    Cs = float(np.asarray(ENV['C_std']).ravel()[0])
    fig, axes = plt.subplots(1, 2, figsize=(9.8, 3.4))
    for ax, y, yb, ys, ylab, fmt in (
            (axes[0], Te, Tb, Ts, '烘房温度 / °C', '均值 {:.2f} ± {:.2f} °C'),
            (axes[1], Ce, Cb, Cs, '烘房水分浓度 / (kg·kg⁻¹)', '均值 {:.4f} ± {:.4f}')):
        ax.axvspan(7200, 14400, color='#999999', alpha=0.12)
        ax.fill_between([7200 / 3600, 14400 / 3600], yb - ys, yb + ys,
                        color=PAL[1], alpha=0.20, linewidth=0)
        ax.plot(t / 3600, y, '-', color=PAL[0], lw=1.3, label='附件1 实测')
        ax.plot([7200 / 3600, 14400 / 3600], [yb, yb], '--', color=PAL[1], lw=1.0,
                label=fmt.format(yb, ys))
        ax.axvline(0.5, color='#333', ls=':', lw=0.9, label='问题1窗口 (0.5 h)')
        style2d(ax, 't / h', ylab)
        ax.legend(fontsize=6.5, loc='lower right')
    axes[0].text(2.2, Te.min(), '平台期 (t ≥ 2 h)', fontsize=6.5, color='#595959', va='bottom')
    axes[1].text(2.2, Ce.min(), '平台期 (t ≥ 2 h)', fontsize=6.5, color='#595959', va='bottom')
    panel(axes[0], 'a'); panel(axes[1], 'b')
    fig.tight_layout()
    save(fig, 'fig6_env')


if __name__ == '__main__':
    fig1_q1(); fig2_q2(); fig3_q3(); fig4_q4(); fig5_validation(); fig6_env()
    print('全部预览图 v2 已输出至 figures/preview/')
