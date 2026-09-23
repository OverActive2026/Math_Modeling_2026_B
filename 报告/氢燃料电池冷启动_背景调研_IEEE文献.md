# 氢燃料电池低温冷启动建模与控制策略研究 —— 问题背景调研（IEEE 文献）

> 调研日期：2026-09-23
> 检索工具：IEEE Xplore Metadata Search API（共 12 组查询，去重后选取 16 篇）
> 说明：本报告仅作背景参考。按赛题备注3要求，正式论文中的模型与公式引用须出自正式发表文献，请勿直接引用 AI 生成内容。

## 1. 赛题结构速览

赛题围绕 PEMFC 低温冷启动的"反应产热 vs 环境散热、反应产水 vs 冻结"动态竞争，层层递进：

| 问题 | 尺度 | 核心任务 |
|---|---|---|
| 问题1 | 一维单电池 | 常温模型 + 冰相（mv/ml/mi 三水形态、相变控制方程、冰占孔修正），用 -20 ℃/-25 ℃ 实验数据校准验证 |
| 问题2 | 5 片电堆 | 计入片间导热与端部对流，优化恒流/线性/阶梯加载 j(t)，约束 qmax=20 C·cm⁻²、jmax=0.5 A·cm⁻² |
| 问题3 | 电堆+辅助加热 | -30 ℃ 下纯预加热 vs 恒功率协同启动，最小化辅助加热总能耗 |
| 问题4 | 电堆+动态控制 | 以 Tk(t)、Vk(t) 反馈的 qk(t) 动态功率控制；三种预冷工况（均匀 / 20 min / 40 min 非均匀温度场） |

## 2. 背景机理要点（与文献互证）

1. **冰堵是冷启动失败的主因之一**：冰点以下，阴极催化层（CL）与气体扩散层（GDL）孔隙被冰占据，氧传输受阻 → 电压快速衰减。失败模式存在"阳极脱水 vs 阴极堵孔"之争（赛题文献[4]，Energies 2020）。
2. **启动成败 = 产热/产水竞争的时序结果**：电流太小产热不足、太大加速产水结冰——这正是问题2 优化 j(t) 的物理动机（赛题文献[5]，JES 2010 的 current ramping 思路）。
3. **电堆尺度的端部效应**：端片散热快、温升慢，最先出现局部结冰与低电压，是启动失败的"短板"（赛题文献[6]，CEJ 2026，明确指出端板单电池热-水失衡并主张通过加载控制缓解）。问题2(2) 要求定位"关键单电池"，答案大概率就是端片。
4. **辅助加热的收益-能耗权衡**：预热或协同加热可扩展启动温度下限，但消耗额外能量 → 问题3/4 的"最小能耗"目标。

## 3. IEEE 文献调研结果（按赛题问题分组）

### 3.1 综述类（写引言/背景用）

- **[R1]** Amamou A. A., Kelouwani S., Boulon L. *A Comprehensive Review of Solutions and Strategies for Cold Start of Automotive PEMFCs*. **IEEE Access**, 2016. DOI: 10.1109/ACCESS.2016.2597058
  系统梳理冰点以下结冰、冻融循环导致的性能衰减，以及各种防冻/快速启动方案。冷启动背景的权威综述，适合作为引言引证。
- **[R2]** Yan H., Dou Y. *Research Progress in Cold Start and Control Strategies of PEMFC*. **IEEE PSGEC**, 2023. DOI: 10.1109/PSGEC58411.2023.10255832
  中文团队综述，聚焦极寒下结冰堵孔、反应气受阻机理与控制策略进展。

### 3.2 冷启动建模（→ 问题1）

- **[R3]** Li Z., Wan Y., Liu J. *PEMFC Cold Start Model Simulation in Low Temperature Environment*. **IEEE AEES**, 2024. DOI: 10.1109/AEES63781.2024.10872517
  低温环境冷启动模型仿真，可作一维瞬态模型建立的直接参照。
- **[R4]** Tatschl R., Ritzberger D., Pötsch C. *Scalable Multi-physics Simulation to Support PEM Fuel Cell System Development*. **IEEE ITEC-India**, 2023. DOI: 10.1109/ITEC-India59098.2023.10471450
  多物理场可扩展仿真框架，说明"机理模型 → 快速降阶模型"的工程化路径。
- **[R5]** Ma Y., Wang X., Tao A. *Onboard Impedance Measurement and Advanced Application of PEMFCs for Cold Start Scenarios*. **IEEE Trans. Industrial Electronics**, 2025. DOI: 10.1109/TIE.2025.3587102
  用 EIS 在线表征冷启动过程内部水含量——可作为模型验证手段的参考（冰/水状态的可观测性）。

### 3.3 自冷启动策略与加载优化（→ 问题2）

- **[R6]** Amamou A., Boulon L., Kelouwani S. *Comparison of self cold start strategies of automotive PEMFC*. **IEEE ICIT**, 2018. DOI: 10.1109/ICIT.2018.8352298
  实验对比恒压 vs 恒流自启动策略，恒压在能耗与系统成本上更优。提示：问题2 除题目指定的三种 j(t) 形式外，可在讨论中提及恒压/最大功率类策略。
- **[R7]** Amamou A., Kandidayeni M., Kelouwani S. *An Online Self Cold Startup Methodology for PEM Fuel Cells in Vehicular Applications*. **IEEE Trans. Vehicular Technology**, 2020. DOI: 10.1109/TVT.2020.3011381
  基于最大功率模式的自适应冷启动（停机排水 + 启动自加热），并做参数实验寻优——"自适应加载"思想与问题2 阶梯/线性优化、问题4 反馈控制直接相关。
- **[R8]** Song R., Wei Z., Sun D. *Adaptive Cold Start of PEM Fuel Cell by Tracking the Maximum Power Point*. **IEEE ITEC+EATS**, 2025. DOI: 10.1109/ITEC63604.2025.11098133
  控制导向模型（水平衡 + 热平衡 + 性能估计三部分耦合），跟踪最大功率点实现快速冷启动；显式引入冰含量（ice fraction）状态。这是与赛题"水-热-电压耦合 + 冰堵约束"框架最接近的 IEEE 文献。
- **[R9]** Gantzer M., Giurgea S., Hissel D. *Statistical investigation of the most influent parameters on the cold start of PEMFC*. **IEEE VPPC**, 2024. DOI: 10.1109/VPPC63154.2024.10755303
  用物理模型 + 统计方法辨识冷启动最敏感的控制与设计参数——可借鉴其"参数敏感性"思路支撑问题2/3 的优化变量选取。

### 3.4 辅助加热与能耗优化（→ 问题3）

- **[R10]** Jin K., Ruan X., Yang M. *Power Management for Fuel-Cell Power System Cold Start*. **IEEE Trans. Power Electronics**, 2009. DOI: 10.1109/TPEL.2009.2020559
  冷启动期间用蓄电池经双向变换器带载、燃料电池自启动完成后滑入/滑出的功率管理方案——辅助能量来源与功率流设计的经典文献。
- **[R11]** Chen F., Zhang H., Li Y. *A Three-Objective Optimization for Fuel Cell-Metal Hydride System Cold Start Using MOPSO Algorithm*. **IEEE ICMA**, 2026. DOI: 10.1109/ICMA69663.2026.11647564
  代理模型（Kriging）+ MOPSO 多目标优化冷启动，说明"高保真模型太贵 → 代理模型 + 群智能优化"是冷启动优化问题的常见解法，可用于问题3/4 的算法选型论证。

### 3.5 动态/智能控制策略（→ 问题4）

- **[R12]** Li W., Dou Y., Li P. *Simulated Annealing Algorithm-based Fuzzy Control Strategy for PEMFC Cold-Start*. **IEEE EI2**, 2025. DOI: 10.1109/EI268505.2025.11425168
  模拟退火优化模糊控制用于冷启动——模糊规则（温度偏低→增强加热等）与赛题问题4 要求的"识别温度偏低/温升不足/冰堵风险/接近成功并分别处置"在形式上高度吻合，可直接作为控制规则设计的参考范式。
- **[R13]** Xu H., Wang L., Guo Z. *Decoupling Multi-Physical Conflicts in PEMFC Cold Start: A Physics-Guided Non-cooperative Game Framework*. **IEEE Trans. Transportation Electrification**, 2026. DOI: 10.1109/TTE.2026.3736413
  物理引导的非合作博弈多目标优化框架，处理"快速活化 vs 空间均匀性"冲突——恰好对应问题4 中"最小能耗 + 电堆最大温差"的多目标权衡，且避免主观赋权。
- **[R14]** Song R., Wei Z., Pan F. *Anodic Cold Start Control of PEM Fuel Cell System With Temperature-Dependent Solenoid Valve Model*. **IEEE Trans. Transportation Electrification**, 2025. DOI: 10.1109/TTE.2025.3526184
  阳极侧（氢气供给）的模型化冷启动控制，说明除电加热外气体供给侧也可参与冷启动控制。

### 3.6 其他相关

- **[R15]** Liu S., Wang Z., Li J. *Thermal Management System and Effect Analysis for Winter Operation of Hydrogen Powered Ships*. **IEEE ICIEA**, 2025. DOI: 10.1109/ICIEA65512.2025.11149251  船舶冬季运行热管理，呼应赛题背景中的"船舶动力、寒区供能"应用场景。
- **[R16]** Guo T., Wang G., Tian Y. *Analysis of Key Technologies for an Advanced Fuel Cell Vehicle*. **ICSGGE**, 2026. DOI: 10.1109/ICSGGE69348.2026.11508988  含 -30 ℃ 低温冷启动测试与整车热管理，可佐证赛题问题3/4 取 -30 ℃ 的工程合理性。

## 4. 对四个问题的启示小结

| 问题 | 可借鉴的文献思路 |
|---|---|
| 1（建模验证） | [R3][R4] 一维多物理建模；[R8] 水/热/电压三耦合 + ice fraction 状态；[R5] 水含量的实验可观测性 |
| 2（加载优化） | [R6][R7][R8] 恒流/恒压/最大功率策略对比；[R9] 参数敏感性定位关键优化变量；端片为最弱环节（赛题文献[6]） |
| 3（辅助加热） | [R10] 辅助能量功率流设计；[R11] 代理模型 + 多目标优化降低计算成本 |
| 4（动态控制） | [R12] 模糊规则 + 退火调参（与题目要求的规则形态一致）；[R13] 博弈论多目标解耦（能耗 vs 温差）；[R7][R8] 自适应/反馈加载 |

## 5. 检索记录

- API：IEEE Xplore Metadata Search API v1（`/api/v1/search/articles`）
- 查询词（12 组）："PEM fuel cell" cold start modeling / fuel cell cold start ice formation / current ramp strategy / stack end cell temperature / auxiliary heating / control strategy / PEMFC cold start optimization / "fuel cell stack" "cold start" / subfreezing startup / preheating energy / experiment voltage temperature / freeze start model validation
- 原始返回已存盘：`.workbuddy/tmp/ieee_raw.json`（第一批），脚本：`.workbuddy/tmp/ieee_search.py`、`ieee_search2.py`
