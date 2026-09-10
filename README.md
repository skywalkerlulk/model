# 2026 CUMCM A题「药材的烘干问题」项目

**状态**: 四个子问题已全部求解并验证（2026-09-10 完成）。截止 9/13 20:00，剩余：论文写作（骨架与摘要草稿已备）。

## 目录结构

```
model/
├── A题/                     题目 PDF + 附件1/2 + 附件3 结果模板 (原始文件, 勿动)
├── algorithms/              ★ 全部算法脚本
│   ├── solver_core.py          统一数值内核: FVM + BDF2 + Newton + 块GS (Q1-Q4 共用)
│   ├── validate_solver.py      验证链: Bessel级数解析解 / 网格收敛 / 时间收敛 / 守恒
│   ├── q1_solve.py             问题1 驱动 → 表1/表2 + result1.xlsx
│   ├── q2_solve.py             问题2 驱动 → 表3/表4 + result2.xlsx + 全流程 npz
│   ├── q3_solve.py             问题3 驱动 → 表5 + result3.xlsx (烘干时长 57.77 h)
│   ├── q4_solve.py             问题4 驱动 → 表6 + result4.xlsx (收缩模型, 51.17 h)
│   ├── verify_q2.py / verify_q4.py   数值稳健性验证 (dt减半 / N加倍)
│   ├── q2_compare.py           变物性必要性对比实验 (湿态冻结 15.92h / 干态冻结 72h不达标)
│   ├── stage6_sensitivity.py   灵敏度分析 (OAT 主导排序 + LHS 联合扰动)
│   └── validation_figs.py      论文验证图 (解析解对比 + 收敛曲线)
├── figures_matlab/          ★ MATLAB 出图脚本 (Nature 风格, run_all_figs.m 一键生成)
│   ├── nm_setup.m              Nature 风格统一设置 (字体/配色/线宽)
│   ├── fig1_q1.m ~ fig6_env.m  六张图 (三维网格/曲面/剖面/验证/环境)
│   └── panel_label / style_ax / export_clean.m  面板字母/轴样式/导出助手
├── debug/                   调试存档 (N=800 Newton 容差问题的定位过程)
├── results/                 结果文件: result1-4.xlsx (附件3模板格式, 4位小数)
│   └── matlab/               ★ 绘图数据 (.mat, 由 export_matlab_data.py 导出)
├── figures/
│   ├── preview/             Nature 风格预览图 (Python 复刻, 立即可看)
│   ├── matlab/              MATLAB 正式出图输出目录 (运行 run_all_figs.m 后生成)
│   └── archive/             旧版 matplotlib 图 (已弃用)
├── paper_workspace/         摘要草稿.md + 论文骨架.md (写作主直接填肉)
└── state/decision_log.json  跨阶段决策日志 (skill 状态文件)
```

## 运行方式

每个脚本头部自动定位项目根（导入 + 数据路径），**从任意目录**运行均可：

```bash
python3 algorithms/validate_solver.py     # 验证内核 (~30s)
python3 algorithms/q1_solve.py            # 问题1 (~17s)
python3 algorithms/q2_solve.py            # 问题2, 72h 全流程 (~60s)
python3 algorithms/q3_solve.py            # 问题3, 依赖 q2 结果
python3 algorithms/q4_solve.py            # 问题4 (~90s)
```

依赖: numpy / scipy / pandas / openpyxl / matplotlib。

## 出图流程 (Nature 风格)

```bash
# ① 导出绘图数据 (一次性)
python3 algorithms/export_matlab_data.py

# ② 本机无 MATLAB: 生成同设计语言的预览图
python3 algorithms/preview_nature_figs.py        # → figures/preview/

# ③ 有 MATLAB 的机器: 一键正式出图 (600 dpi PNG + 矢量 PDF)
matlab -batch "cd('figures_matlab'); run_all_figs"   # → figures/matlab/
# 或在 MATLAB 界面中打开 figures_matlab/run_all_figs.m 直接运行
```

## 核心结果

| 问题 | 模型 | 烘干时长 | 验证 |
|---|---|---|---|
| Q1 | 常物性轴对称传热传质 (VC-HTMT-0) | — | 解析解误差 3e-6, 守恒残差 7.6e-7 |
| Q2 | 变物性双向耦合 (VC-HTMT-1, 附录3) | — | dt减半 Δt_end=0, N=800 差 0.18h |
| Q3 | 判据 max C<0.15 扫描 | **57.77 h** (207960 s) | 附件2 收缩平台 70h 同尺度 |
| Q4 | Landau 收缩移动边界 (VC-HTMT-L, 附录4+附件2) | **51.17 h** (184200 s) | dt减半 Δ=0, N=800 差 0.017h |

灵敏度: 恒温段温度为主导参数 (OAT ±10% → ±8~10h, Arrhenius 放大); h 影响 ≈0。
详细建模方案与算法: 见桌面 `A题建模方案与算法.md`。

## 关键建模决策 (详见 state/decision_log.json)

1. 长圆柱近似 (L/R=12.5) → 一维径向模型
2. 附件1 平台期检测 (t≥7200s: 49.97±0.17°C, 0.0500 kg/kg) → 恒温段边界条件外推
3. Q4 物质坐标 ξ=r/R(t), V(t) 为时间导数**系数** (非 d(V·C)/dt, 后者引入虚假压缩项)
4. 变步长 BDF2: 0.25s/1s/5s 分档 + 切换点 BE 重启; Newton 相对容差 1e-12
