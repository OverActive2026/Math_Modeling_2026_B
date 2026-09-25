# Q2_Optimized_Step_Loading

这是可单独放在 MATLAB 当前文件夹运行的五片燃料电池阶梯加载包，包含真实物理模型、
未优化四档基准、五档九变量优化算法、已得到的结果及独立复算入口。建议 MATLAB R2024b；
优化需要 Statistics and Machine Learning Toolbox（`lhsdesign`、`fitrgp`）。

在 MATLAB 将当前文件夹切换到本文件夹后：

1. 运行 `run_q2_step_loading_comparison`：重算未优化四档与优化五档，输出对比数值、
   `Q2_step_loading_comparison.png` 和 `.mat`。这是写论文最直接的入口。
2. 运行 `verify_q2_step5_full_search`：对五档最优解做快速网格及更严格容差复核。
3. 如需继续寻优，运行 `report = run_q2_step5_full_search(20)`。首次在新路径执行会重新
   仿真已发表的最优点，而不是盲用旧目录缓存；之后相同源码及设置下可复用本地缓存。
   一次真实模型搜索可能耗时较长。`run_q2_step_current` 可单独运行和修改四档基准。

九维向量的前五项是五档电流密度，后四项是四个切换时刻；**九项全部参与优化**。
电流范围 0--0.5 A/cm²，切换时刻的数值搜索窗分别是 `[1,8]`、
`[8.1,16.9]`、`[17,23.9]`、`[24,38]` 秒。这些是搜索范围，不是题目给定的
物理界限；若最优点碰到边界，应扩窗复核。模型为 `-10 °C` 起始、快速网格
`[2 3 4 4 2]`、`gammaIce=3.5`、`kFreeze=0.4 s^-1`。模拟时间和电脑耗时不是一回事。

当前结果：未优化四档为 **45.8918 s**；优化五档为 **36.0266 s**，严格容差复算
**36.0275 s**。优化五档的电流 `[0.1752,0.2914,0.2427,0.2618,0.4598]`
A/cm²，切换时刻 `[2.5713,12.1507,20.5979,27.5905]` 秒。严格复算的最低
单电池电压 0.314724 V、累计电荷 11.002117 C/cm²、启动时最大冰量 0.192016。
完整结果在 `q2_step5_full_report.mat` 和 `q2_step5_full_verification.mat`。

论文中请区分两种改进来源：**45.8918→36.0266 s 同时改变了曲线结构
（四档变五档）和参数**，不能全归因于算法；同属五档曲线，原固定前缀方案
36.5010 s，九变量联合优化达到 36.0266 s。算法借鉴 FuRBO 的可行性驱动
信赖域与 LogEI 排序，但含工程改动，不能称为论文方法的严格复现。
结果只在当前模型、参数与快速网格上验证；30 秒以下尚未实现，也尚未完成
完整网格和实验数据的最终验证。

算法参考：Ascia 等，*Feasibility-Driven Trust Region Bayesian Optimization*
（https://arxiv.org/abs/2506.14619）；Ament 等，*Unexpected Improvements to
Expected Improvement for Bayesian Optimization*
（https://arxiv.org/abs/2310.20708）。论文写作应明确“受其启发”，而非严格复现。
