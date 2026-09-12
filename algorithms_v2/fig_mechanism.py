# -*- coding: utf-8 -*-
"""
fig_mechanism.py — 物理机理总览示意图 (参照获奖论文三面板结构, 独立视觉语言)
======================================================================
布局逻辑 (借鉴): (a) 烘干室中位置 → (b) 单根药材表面传热传质 → (c) 横截面径向模型+控制方程
视觉区分 (不抄袭): 手绘双线框 / 低饱和六色板 / 一阶惯性环境辨识曲线(本文特色) /
                 谱配点节点聚类(本文特色) / β 记号体系
输出: figures/v2/fig0_mechanism.png (Nature 风格)
"""

import sys, os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
sys.path.insert(0, os.path.join(_ROOT, 'algorithms_v2'))
os.chdir(_ROOT)

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import (Rectangle, Circle, FancyArrow, FancyBboxPatch,
                                Wedge, PathPatch)
from matplotlib.path import Path

plt.rcParams.update({'font.sans-serif': ['Arial Unicode MS', 'PingFang SC', 'Heiti TC'],
                     'axes.unicode_minus': False, 'font.size': 8,
                     'axes.linewidth': 0.6, 'legend.frameon': False,
                     'figure.facecolor': 'white'})
PAL = ['#0C7BDC', '#D64541', '#3CB478', '#8C78C8', '#C8AA3C', '#787878']
INK = '#3A3F47'          # 手绘线稿墨色
GRID = dict(alpha=0.12, color='#595959', linewidth=0.5)
OUT = 'figures/v2'
os.makedirs(OUT, exist_ok=True)


def hdbox(ax, xy, w, h, fc='none', ec=INK, lw=1.0, ls='-'):
    """双线铅笔描边矩形 (手绘感)"""
    ax.add_patch(Rectangle(xy, w, h, fill=(fc != 'none'), facecolor=fc,
                           edgecolor=ec, lw=0.35, ls=ls, zorder=1))
    ax.add_patch(Rectangle((xy[0] + 1.2, xy[1] + 1.2), w - 2.4, h - 2.4,
                           fill=False, edgecolor=ec, lw=0.35, ls=ls, zorder=1))


def hdarrow(ax, xy0, xy1, color=INK, lw=1.3, style='-|>', ms=11):
    ax.annotate('', xy=xy1, xytext=xy0,
                arrowprops=dict(arrowstyle=style, color=color, lw=lw,
                                mutation_scale=ms, shrinkA=0, shrinkB=0))


def cyl(ax, cx, cy, w, h, ec=INK, fc='#DCE7F2', shrink=False):
    """圆柱 (胶囊形)"""
    r = h / 2
    body = Rectangle((cx - w / 2, cy - r), w, 2 * r, fill=True, facecolor=fc,
                     edgecolor=ec, lw=0.8, zorder=2)
    ax.add_patch(body)
    for sgn in (-1, 1):
        ax.add_patch(Wedge((cx + sgn * w / 2, cy), r, 90, 270, width=0,
                           facecolor='none', edgecolor=ec, lw=0.8, zorder=3))


def zz(ax, x0, y0, x1, n=5, color=INK, lw=1.0):
    """锯齿线 (加热器符号)"""
    pts = [(x0, y0)]
    dx = (x1 - x0) / (2 * n)
    for i in range(2 * n):
        pts.append((x0 + (i + 1) * dx, y0 + (0.9 if i % 2 == 0 else -0.9) * 3.2))
    ax.add_patch(PathPatch(Path(pts), fill=False, edgecolor=color, lw=lw, zorder=2))


def fan(ax, cx, cy, r=5.5, color=INK):
    """循环风机 (圆+叶片)"""
    ax.add_patch(Circle((cx, cy), r, fill=False, edgecolor=color, lw=0.7, zorder=2))
    for k in range(3):
        th = np.deg2rad(30 + k * 120)
        x0, y0 = cx + r * 0.35 * np.cos(th), cy + r * 0.35 * np.sin(th)
        x1, y1 = cx + r * 0.95 * np.cos(th), cy + r * 0.95 * np.sin(th)
        ax.plot([x0, x1], [y0, y1], color=color, lw=0.8, zorder=2)


def main():
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 4.0),
                             gridspec_kw={'width_ratios': [1.05, 1.0, 1.35]})

    # ============ (a) 烘干室 ============
    ax = axes[0]
    ax.set_xlim(0, 100); ax.set_ylim(0, 78); ax.axis('off')
    hdbox(ax, (4, 6), 92, 64, fc='#F7F9FB')
    ax.text(50, 66.5, '烘 干 室', fontsize=9.5, ha='center', color=INK)
    # 加热器 (左下)
    zz(ax, 12, 14, 27, color=PAL[1])
    ax.text(19.5, 9, '加热器', fontsize=7, ha='center', color=PAL[1])
    hdarrow(ax, (7, 16.5), (10.5, 16.5), color=PAL[1], lw=1.2)
    # 循环风机 (右上)
    fan(ax, 82, 58)
    ax.text(82, 49, '循环风机', fontsize=7, ha='center', color=INK)
    # 新风入口 / 排出口
    hdarrow(ax, (2, 40), (8.5, 40), color=PAL[3], lw=1.3)
    ax.text(4, 44, '新风', fontsize=7, color=PAL[3])
    hdarrow(ax, (92, 36), (97.5, 36), color=PAL[3], lw=1.3)
    ax.text(90, 32, '湿气排出', fontsize=7, color=PAL[3])
    # 药材 (胶囊形) 置于中部
    cyl(ax, 52, 32, 26, 15)
    ax.text(52, 27.8, '圆柱药材', fontsize=8, ha='center')
    # 环境标识 + 辨识曲线 (本文特色)
    ax.text(52, 55.5, '热湿空气:  $T_{env}(t)$, $C_{env}(t)$', fontsize=8.5,
            ha='center', color=PAL[0])
    axin = fig.add_axes([0.045, 0.46, 0.16, 0.24])
    tt = np.linspace(0, 24, 200)
    Tfit = 49.97 - (49.97 - 28) * np.exp(-tt / 1.776)
    axin.plot(tt, Tfit, color=PAL[0], lw=1.2)
    axin.axvspan(2, 24, color='#999', alpha=0.15, linewidth=0)
    axin.text(13, 31, '平台期 49.97±0.17°C', fontsize=5.2, ha='center', color='#555')
    axin.set_xlim(0, 24); axin.set_ylim(26, 52)
    axin.set_xlabel('t / h', fontsize=5.5); axin.set_ylabel('$T_{env}$ / °C', fontsize=5.5)
    axin.tick_params(labelsize=4.5, width=0.5)
    axin.set_title('环境辨识 (一阶惯性)', fontsize=5.8, color=PAL[0])
    for s in axin.spines.values():
        s.set_linewidth(0.5)
    ax.text(50, 3, '(a) 药材在烘干室中的位置', fontsize=8.5, ha='center', color=INK)

    # ============ (b) 单根药材 ============
    ax = axes[1]
    ax.set_xlim(0, 100); ax.set_ylim(0, 78); ax.axis('off')
    cyl(ax, 50, 44, 58, 20)
    # 尺寸标注
    ax.annotate('', xy=(21, 38), xytext=(79, 38),
                arrowprops=dict(arrowstyle='<->', color=INK, lw=0.8))
    ax.text(50, 35.5, 'L = 25 cm（轴向长度保持不变）', fontsize=7.5, ha='center')
    ax.annotate('', xy=(54, 54), xytext=(54, 44.5),
                arrowprops=dict(arrowstyle='<->', color=INK, lw=0.8))
    ax.text(56.5, 48.5, 'R(t)：2.00→1.198 cm', fontsize=7, va='center', color=PAL[1])
    # 收缩箭头
    for ang in (-40, 40):
        x0 = 50 + 24 * np.cos(np.deg2rad(ang))
        y0 = 44 + 24 * np.sin(np.deg2rad(ang))
        ax.annotate('', xy=(50 + 18 * np.cos(np.deg2rad(ang)), 44 + 18 * np.sin(np.deg2rad(ang))),
                    xytext=(x0, y0),
                    arrowprops=dict(arrowstyle='->', color=PAL[1], lw=1.0,
                                    mutation_scale=9))
    # 表面通量
    hdarrow(ax, (79, 58), (86, 64), color=PAL[1], lw=1.4)
    ax.text(86, 66.5, '$q_T = h\\,(T-T_{env})$', fontsize=7, color=PAL[1], ha='center')
    hdarrow(ax, (21, 58), (14, 64), color=PAL[0], lw=1.4)
    ax.text(14, 66.5, '$q_C = \\beta\\,(C-C_{env})$', fontsize=7, color=PAL[0], ha='center')
    ax.text(50, 60, '周围热湿空气 $T_{env}(t)$, $C_{env}(t)$', fontsize=7.5,
            ha='center', color='#555')
    ax.text(50, 3, '(b) 单根药材的表面传热传质', fontsize=8.5, ha='center', color=INK)

    # ============ (c) 横截面径向模型 ============
    ax = axes[2]
    ax.set_xlim(0, 112); ax.set_ylim(0, 78); ax.axis('off')
    # 截面圆
    ax.add_patch(Circle((36, 42), 24, fill=True, facecolor='#DCE7F2',
                        edgecolor=INK, lw=0.8, zorder=1))
    # 谱节点 (Chebyshev 聚类 — 本文特色)
    x_gl = np.cos(np.pi * np.arange(49) / 48)
    r_nodes = 24 * (1 - x_gl) / 2
    ax.plot(36 + r_nodes, np.full(49, 42), 'o', ms=1.6, color=PAL[0], zorder=4)
    ax.plot([36, 60], [42, 42], color=INK, lw=0.7, zorder=3)
    ax.plot(36, 42, 'o', ms=4, color=PAL[3], zorder=5)
    ax.text(36, 36.5, '轴心 $r=0$（对称）', fontsize=7, ha='center', color=PAL[3])
    ax.annotate('', xy=(48, 51), xytext=(38, 44.5),
                arrowprops=dict(arrowstyle='->', color=PAL[0], lw=0.9))
    ax.text(49, 52.5, '径向坐标 $r$', fontsize=7, color=PAL[0])
    ax.annotate('', xy=(61, 40), xytext=(61, 44),
                arrowprops=dict(arrowstyle='->', color=PAL[1], lw=0.9))
    ax.text(63, 38, '表面 $r=R(t)$（Robin 条件）', fontsize=7, color=PAL[1])
    # 控制方程
    ax.text(78, 68, '内部传热:', fontsize=7.5, ha='left', color=INK)
    ax.text(78, 63, r'$\rho c_p\frac{\partial T}{\partial t}=\frac{1}{r}\frac{\partial}{\partial r}\left(r k\frac{\partial T}{\partial r}\right)$',
            fontsize=8, ha='left')
    ax.text(78, 55.5, '内部传质:', fontsize=7.5, ha='left', color=INK)
    ax.text(78, 50.5, r'$\frac{\partial C}{\partial t}=\frac{1}{r}\frac{\partial}{\partial r}\left(r D(C,T)\frac{\partial C}{\partial r}\right)$',
            fontsize=8, ha='left')
    ax.text(78, 42.5, '表面条件:', fontsize=7.5, ha='left', color=INK)
    ax.text(78, 37.5, r'$-k\frac{\partial T}{\partial r}=h(T-T_{env})$', fontsize=7.5, ha='left')
    ax.text(78, 32.5, r'$-D\frac{\partial C}{\partial r}=\beta(C-C_{env})$', fontsize=7.5, ha='left')
    ax.text(56, 3, '(c) 横截面径向模型与谱离散节点', fontsize=8.5, ha='center', color=INK)

    fig.tight_layout()
    fig.savefig(f'{OUT}/fig0_mechanism.png', dpi=300, bbox_inches='tight')
    fig.savefig(f'{OUT}/fig0_mechanism.pdf', bbox_inches='tight')
    plt.close(fig)
    print('fig0_mechanism 已输出')


if __name__ == '__main__':
    main()
