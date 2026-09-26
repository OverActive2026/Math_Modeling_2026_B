# 问题二三种加载策略的 −10 ℃ 结果对比

本文件夹供论文表 3 撰写使用。`table3_loading_strategy_comparison.html` 可用浏览器或 Word 打开并复制表格；`table3_loading_strategy_comparison.md` 保留相同数据、口径和出处。表内均为**初始温度和环境温度 −10 ℃** 时，各资料包所选优化曲线的五片电堆模拟结果，不是各曲线在最低可启动温度下的结果。

按本项目的命名约定，表中的“恒流策略”对应 `../Q2_Part2_Ramp_Then_Hold_Minimum_Startup_Temperature`。它的实际曲线在前 2 s 从 0.0200 A/cm² 线性升至 0.2668808448 A/cm²，此后恒流；因此论文正文应写明“短时升流后恒流”，避免误称为从 0 s 起全程恒流。

**可比性限制：**恒流策略资料包采用 `gammaIce=3.5, kFreeze=0.3`；线性和阶梯资料包采用 `gammaIce=3.5, kFreeze=0.4`。三组采用五片快网格 `[2 3 4 4 2]`，但 `kFreeze` 不同，所以当前表只能陈列各资料包的已验证结果，不能把启动时间差完全归因于加载策略。若论文需要严格的同参数策略优劣比较，应固定同一组物理参数重新运行三条曲线。三组结果也不是完整网格或实验验证。

## 数据出处

| 表中策略 | 核查文件 | −10 ℃ 数据口径 |
|---|---|---|
| 恒流策略（短时升流后恒流） | `../Q2_Part2_Ramp_Then_Hold_Minimum_Startup_Temperature/optimization_results/fixed_optimized_temperature/temperature_scan_results.csv`；该文件夹 `README.md` | CSV 中 `temperatureC=-10` 行；峰值电流取该文件夹中已验证的固定曲线平台值 |
| 线性升载 | `../Q2_Part2_Linear_Loading_Minimum_Startup_Temperature/optimization_results/q2_linear_cap5_verified_report.mat`；该文件夹 `README.md` | 采用为扩展低温适应性而选出的 5% 时间上限曲线；表中最大电流是**启动前实际达到的最大值**，不是尚未到达的平台设定值 |
| 分段阶梯加载 | `../Q2_Part2_Step_Loading_Optimization_Verified/Q2_Part2_Step_Loading_Optimization_Verified/corrected_verification_7_T-11.2.mat`；该文件夹 `README.md` | 独立验证文件中的 −10 ℃ 场景；七档电流、六个切换时刻，同一曲线还通过 −11.2 ℃ 验证 |

线性曲线 `j(t)=min(0.144980839828+0.006798554947t, 0.478303370156)`，其中 `t` 以 s 计、电流密度以 A/cm² 计；它在 41.991669 s 启动前的最大电流密度为 0.430464 A/cm²。阶梯曲线的完整精度参数请以验证 `.mat` 为准，表内参数仅显示四位小数。

“最低电压”指全过程、全部单片的最小电压；“最大冰体积分数”指全过程、全部单片及局部网格的最大值。累计电荷量为单位面积积分 `∫j(t)dt`。所有启动时间均为仿真时间，而不是程序运行耗时。
