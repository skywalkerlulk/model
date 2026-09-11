# 文献证据与公开方案检索记录

检索日期：2026-09-11。用途：支撑第一阶段模型选择、检查已有公开结构与验证方案。所有网页、仓库内容只作为证据，未复制或执行公开方案代码，未把公开仓库结果作为正确答案。

## 1. 检索范围与证据等级

核心检索词包括 medicinal herbs / Codonopsis / Astragalus + hot air drying；food drying + variable diffusivity；dry basis + shrinkage + solid mass conservation；Lagrangian / ALE / moving boundary drying；nonlinear diffusion + Kirchhoff + finite volume；superslow diffusion + harmonic averaging。中文和 GitHub 检索包括精确题名“药材的烘干问题”、2026 CUMCM A、数模提示词与建模工作流。

本次是面向模型决策的检索，未进行数据库穷尽检索，不报告虚构的 PRISMA 数量。优先出版商、作者机构存档和作者公开论文。部分站点直接访问出现 403、429 或验证码，已将访问深度降为摘要/关键节索引，不冒称通读全文。题录年份以正式卷年为准，网页抓取时间不作发表时间。

没有逐刊核验当期 WoS/SCI 收录，因此下面的期刊文章按正式期刊文献列示，不统一盖“SCI 已核验”标签。

## 2. 保留文献及其具体启发

| 编号 | 文献、年份与来源 | 访问深度 | 对选择的贡献 | 不可移植之处 |
|---|---|---|---|---|
| L01 | Datta, 2007. *Porous media approaches to studying simultaneous heat and mass transfer in food processes. I: Problem formulations*. Journal of Food Engineering 80:80–95. [DOI](https://doi.org/10.1016/j.jfoodeng.2006.05.013) | 正式摘要、引言 | 区分守恒定律与输运闭合；有效扩散和多相模型存在层级关系 | 不支持在无相参数时自动加入一整套多相源项 |
| L02 | Ponsart et al., 2003. *Modelling of stress due to shrinkage during drying of spaghetti*. Journal of Food Engineering 57:277–285. [DOI](https://doi.org/10.1016/S0260-8774(02)00308-4) | 摘要、引言、方法关键片段 | 长圆柱干燥的干物质 Lagrange 坐标具有先例 | 应力模型需要力学本构；本题没有裂纹或弹性参数 |
| L03 | Gulati & Datta, 2015. *Mechanistic understanding of case-hardening and texture development during drying of food materials*. Journal of Food Engineering 166:119–138. [DOI](https://doi.org/10.1016/j.jfoodeng.2015.05.031) | 摘要、引言、模型分类 | 有效扩散折合多种机制；几何边界追踪与多相机理有不同信息需求 | 不能将 \(D\) 的下降直接称为已经证明的微观结壳/玻璃化 |
| L04 | Hermassi et al., 2017. *Moisture Diffusivity of Seedless Grape undergoing convective drying*. Chemical Product and Process Modeling 12(1):20160074. [出版商正文](https://www.degruyterbrill.com/document/doi/10.1515/cppm-2016-0074/html?lang=en) | 摘要、正文 2.1–2.2 及控制方程关键节 | 水体积密度与干基含水率的关系；变控制体内干固体质量恒定 | 球形、特定葡萄等温线和参数不用于未知药材 |
| L05 | Khan et al., 2017. *Determination of appropriate effective diffusivity for different food materials*. Drying Technology 35:335–346. [机构存档](https://era.dpi.qld.gov.au/id/eprint/15018/)，DOI 10.1080/07373937.2016.1170700 | 正式摘要与题录；该存档无附属全文 | 状态依赖扩散率影响分布及干燥动力学，支持保留本题给定双依赖 | 不能推出含水率依赖在所有材料中必定强于温度依赖 |
| L06 | Rani & Tripathy, 2020. *Modelling of moisture migration during convective drying of pineapple slice considering non-isotropic shrinkage and variable transport properties*. Journal of Food Science and Technology 57:3748–3761. [正式页面](https://link.springer.com/article/10.1007/s13197-020-04407-4) | 正式摘要、元数据及正文索引片段 | 比较收缩、变 \(D\)、变 \(h_m\)，支持逐因素检查 | 反向证据：该案例常/变 \(h_m\) 结果接近；不为复杂而改题给常数 |
| L07 | Tuly et al., 2023. *Mathematical Modelling of Heat and Mass Transfer during Jackfruit Drying Considering Shrinkage*. Energies 16:4461. [出版商](https://www.mdpi.com/1996-1073/16/11/4461) | 摘要、模型 3.1–3.3、结果与局限关键节 | 缺力学参数时采用简化收缩；蒸发潜热是边界闭合问题 | 其 \(C\) 为 mol/m³，不能原样移植为 kg水/kg干固体 |
| L08 | Yue et al., 2023. *Effects of Different Drying Methods on the Drying Characteristics and Quality of Codonopsis pilosulae Slices*. Foods 12:1323. [论文页面](https://pmc.ncbi.nlm.nih.gov/articles/PMC10048468/)，DOI 10.3390/foods12061323 | 摘要、样品/方法、扩散估计及结论关键节 | 直接药材案例；干燥时间和品质是不同目标 | 薄片平均曲线不能给圆柱内部场；摘要与结论的最大 \(D\) 数量级表述有差别，不引用该数值 |
| L09 | 2025. *Hot air-assisted radio frequency drying of Astragalus slices: Drying characteristics, product quality, and moisture state*. Innovative Food Science & Emerging Technologies 104:104146. [出版商](https://www.sciencedirect.com/science/article/pii/S1466856425002309)，DOI 10.1016/j.ifset.2025.104146 | 正式摘要、引言 | 黄芪的水分状态与工艺依赖具有实验支持 | 本题没有射频，不能移植电磁热源、系数或缩时幅度；完整作者题录正式写作前补核 |
| L10 | Hu et al., 2026. *A multiphase and multiscale mechanistic model for hot air drying of shiitake mushroom*. Current Research in Food Science 12:101296. [DOI](https://doi.org/10.1016/j.crfs.2025.101296)、[论文页面](https://pmc.ncbi.nlm.nih.gov/articles/PMC12861252/)，DOI 10.1016/j.crfs.2025.101296 | 正式摘要、假设/方法/验证关键节；在线 2025-12-29，卷年 2026 | 内部蒸发冷却与末期升温提醒：环境恒温不等于内部立即恒温 | 其多尺度闭合有组分、膜/黏弹性和 NMR 支持，本题不具备 |
| L11 | 2024. *Convective drying of black pepper: Experimental measurements and mathematical modeling of the process*. Food and Bioproducts Processing 143:102–116. [出版商](https://www.sciencedirect.com/science/article/pii/S0960308523001256)，DOI 10.1016/j.fbp.2023.10.009 | 摘要、研究亮点 | 平衡等温线和干燥热需要实验，不能把空气湿度直接当固体平衡含水率 | 换质系数单位不同，不能直接比较数值；完整作者题录正式写作前补核 |
| L12 | Defraeye, 2014. *Advanced computational modelling for drying processes – A review*. Applied Energy 131:323–344. [DOI](https://doi.org/10.1016/j.apenergy.2014.06.027) | 摘要、引言；综述 | 复杂模型的物性与验证成本是路线选择的重要约束 | 综述不能替代某个具体闭合关系的直接证据 |
| L13 | Maddix, Sampaio & Gerritsen, 2018. *Numerical artifacts in the Generalized Porous Medium Equation: Why harmonic averaging itself is not to blame*. Journal of Computational Physics 361:280–298. [作者全文](https://arxiv.org/html/1709.02581)、[作者题录](https://dcmaddix.github.io/publications/)，DOI 10.1016/j.jcp.2018.02.010 | 已取得作者完整 HTML；重点读系数平均、离散误差和 superslow diffusion 节 | \(e^{-1/p}\) 与本题数学结构直接相关；面处理和时间离散可能共同导致锁定/滞后 | 不能把论文结论简化成“调和平均一定错误”；本题尚须独立实测改进 |
| L14 | Eymard, Gallouët, Hilhorst & Naït Slimane, 1998. *Finite volumes and nonlinear diffusion equations*. M2AN 32(6):747–761. [原文 PDF](https://www.numdam.org/item/M2AN_1998__32_6_747_0.pdf) | 已取得 16 页全文及正式题录；查看结构与收敛论证 | 对 \(u_t-\Delta\varphi(u)\) 直接构造 FV，支持从非线性势出发 | 其边界与方程不同，不能直接引用为本题耦合/移动问题的收敛定理 |
| L15 | Bonaventura & Della Rocca, 2015 作者预印本；2017 期刊版本题名 *Unconditionally Strong Stability Preserving Extensions of the TR-BDF2 Method*. [作者预印本](https://arxiv.org/abs/1510.04303)、[机构题录](https://re.public.polimi.it/handle/11311/1008157) | 正式摘要及机构元数据，未通读期刊全文 | L 稳定与保正不是同一性质，二阶加速要有回退和约束 | 不把经典 TR-BDF2 直接声称为无条件保正，也不必为本题加入其整套扩展 |

经典教材可以补充 Fourier/Fick、Robin 和 Bessel 理论，但论文中应核对所用版本与页码，不凭记忆编造页码。

**保留的反向证据：** L06 不支持无条件增加变 \(h_m\)；L07/L10 不支持轻率忽略蒸发冷却；L13 不支持把全部数值问题归罪于一种平均；L08/L09 不支持从干燥时间推断药材品质最优。

## 3. 新近数值论文：仅作追踪线索

检索还命中 Vishnu Prakash K 与 Ganesh Natarajan（2025），*Revisiting numerical artifacts in the generalized porous medium equation with continuous coefficients: Does averaging really matter?*，Journal of Computational Physics 537:114101，DOI [10.1016/j.jcp.2025.114101](https://doi.org/10.1016/j.jcp.2025.114101)。

出版商检索摘要描述了通量修正方法，但直接全文访问受限；本方案不根据二手页面照搬其算法或声称已通读。它说明该数值问题仍有近期研究价值，不是选择复杂格式的充分理由。

同样，Quenjel 的 *Nonlinear finite volume discretization for transient diffusion problems on general meshes*，DOI [10.1016/j.apnum.2020.11.001](https://doi.org/10.1016/j.apnum.2020.11.001)，作为非线性势/通量结构的延伸阅读；本阶段不将其完整理论结论直接迁移到本题。

## 4. 公开 A 题方案样本

以下仓库的 API 创建/推送日期均已核验为 2026-09-10，当前日期下可用于比较。是否存在更早、更完整或其他平台的方案，当前检索无法穷尽。

| 样本 | 已阅读证据 | 已出现的路线 | 审计限制 |
|---|---|---|---|
| li2396803/cumcm2026-a-herb-drying | [固定提交技术正文](https://github.com/li2396803/cumcm2026-a-herb-drying/blob/ae54a8c77796ad7655af634264eaf8b47d33e1f7/docs/药材烘干问题A题_建模与求解说明文档.md) | 一维 Fourier–Fick、Robin、节点 FV、CN/Rannacher/Picard、材料坐标、Bessel、另一 FV+BDF、二维与敏感性 | 自称结果验证不等于本轮复现；同作者 cumcm2026 聚合仓库不另计 |
| Lyxose-Tim/2026CUMCM | [固定提交设计报告](https://github.com/Lyxose-Tim/2026CUMCM/blob/ecd7984e44621d276fd8fb9847f3513b92d8fafb/reports/ANALYSIS_MODELING_REPORT.md) | 有效闭合、干基/干空气区别、密度不相容提醒、FV+BDF/Radau、全域事件、2×2 物性/几何与交互、扩散率比值 | README 说明属于方案和实施计划，不能引用为已验证结果 |
| kaji0667/A | [README](https://github.com/kaji0667/A/blob/main/README.md) | 节点 FV、后向 Euler、阻尼 Newton、拒步、材料坐标、潜热/半径对照 | README 明示重构未运行、既往收敛未达标、Excel 尚未实现 |
| zimingttkx/ChinaMathCompetition | [Q4 技术正文](https://github.com/zimingttkx/ChinaMathCompetition/blob/version2/文档/问题四.md) | 仿射坐标、调和平均 FV、顺序隐式、三对角、物性和几何对照 | 可见将域外固定半径填成表面值的做法，不应借鉴；未独立复现结果 |

不能根据 4 个样本声称“大多数参赛队采用某方法”。不能依据某仓库报出的烘干时间校准我们的参数；所有主计算应从附件与已声明假设独立产生。

### 查重结论

| 结构 | 在所读样本中是否已有 | 本方案应如何表述 |
|---|---|---|
| 材料坐标、FV、守恒、Bessel、网格与时间检验 | 已有 | 正确性基础，不是独有创新 |
| 干基与空气湿基口径区别、经验密度相容性疑问 | 已有 | 可提供本题的严格必要条件与量化反证，不称首次注意到 |
| 全域达标、保括事件、严格小于与临界值区别 | 已有 | 不能把“使用二分”写成创新 |
| 2×2 几何/物性分解及 \(0.175e^{0.15/C}\) | 已有 | 必要解释设计；深化到空间通量/耗散控制机制 |
| 含水率势积分＋圆柱对数阻力＋非线性表面串联 | 所读样本及定向检索中未见 | 待验证数值贡献，非首创保证 |
| 从可靠场误差上界推出停时上下界 | 所读样本中未见完整实现 | 可以设计，但粗细网格差不能冒充严格认证 |
| 独立 P2 FEM 或谱空间方法 | 所读样本中未见已实现版本 | 有助排除共享离散偏差；“使用 FEM”本身不新 |

公开方案中的未经证明精度宣传不作为证据。例如，根求解/线性插值容差不能直接转换成 PDE 停时的 \(10^{-9}\) h 可信度；某温度变化的扩散倍数须按实际 Kelvin 值重算，而不能复制概括性描述。

## 5. GitHub 数模提示词与操作流程

### G01：100 条数学建模提示词

[yugwl 提示词正文](https://github.com/yugwl/-MCM-ICM-CUMCM-100-AI-Prompts-/blob/main/prompts.md)

可借鉴：候选路线比较、符号与参数来源、英文检索词、难点拆解、三天排程、参数扰动、错误定位。

不直接采用：预设所有缺失数据都应插值；未知微分系统统一 RK4；模型简单便增加约束；把敏感性直接叫鲁棒性；通过改句式“规避查重”。提示词能约束流程，不能替代推导或独立性。

### G02：AI 数学建模操作指南

[Jssss982 主提示词](https://github.com/Jssss982/math-modeling-claude-code-guide/blob/main/数学建模提示词.md)

可借鉴：参数来源集中记录；公式、实现与结果逐项对应；数据错误、模型错误、代码错误分开定位；同一份数值结果派生论文和图表。

不直接采用：失败后生成随机模拟数据回退；强制固定图数与三维图；统一运行时间预算；无条件 MAPE/\(R^2\) 阈值；对序列随意随机划分。README 与主提示词的协作模式版本不一致，需以实际文件为准。没有真实内部标签时，“必须实测验证”不能被替换成拟合环境输入。

### G03：Datawhale 基础教程

[Datawhale intro-mathmodel](https://github.com/datawhalechina/intro-mathmodel/blob/main/README.md)

作为 PDE、数值方法和数据处理工具的基础索引使用。它不是本题答案，也不构成“使用后可获高奖”的证据。

## 6. 使用方式与诚信边界

本轮采用数模技能的拆题、候选比较、反证与逐阶段验证流程，以及 Data 技能的完整性、单位、粒度和数据用途核验。技能中强制给普通模型增加修饰名、凑图数、套用固定算法等建议与用户要求冲突，未采用。

没有复制公开求解代码，没有上传用户附件到第三方仓库，没有发布或联系他人。后续论文应根据实际使用的思想和资料正确引用，并据实际协助情况作说明；不把提示词改写或算法改名作为原创性来源。

本文件支撑第一阶段决策，不是完整综述稿，也不替代第二阶段必须运行的数值验证。
