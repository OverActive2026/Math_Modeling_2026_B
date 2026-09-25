# 恒流固定参考组与文献算法

## 30次真实评价结果

| 指标 | 固定参考组 | 优化组 |
|---|---:|---:|
| 平台电流密度 / (A/cm²) | 0.259937538852887 | 0.266880844760886 |
| 启动时间 / s | 68.2624274374444 | 65.0058211159030 |
| 全程最低单片电压 / V | 0.3324556900 | 0.3240780645 |

优化组启动时间缩短3.2566063215秒，即4.7707156686%。最佳点来自第28次真实评价，属于第13次起的FuRBO＋LogCEI自适应选点阶段；独立重跑得到相同启动结果。30次中18次成功、12次物理失败、0次数值未完成。最终平台电流未达到0.5上限，因为更高电流可能使单片电压提前触限。结果为本轮预算内的最佳实测值，不声称全局最优。

完整报告位于`optimization_results/constant_fixed_reference_kfreeze03/tp8d194731_ec31_4801_8c43_3e9515cb2c92_report.mat`，精确数值见同目录`comparison_30eval.json`，图为`reference_vs_optimized_30eval.png`。

在MATLAB中切换到本文件夹 `Q2_Optimized_Constant_Loading` 后，运行`run_q2_constant_fixed_reference`可单独复现固定参考；运行`result=run_q2_constant_best_current()`可复现优化组并显示对比图；运行`report=run_q2_constant_fixed_reference_optimization()`会重新发起30次搜索。

固定参考脚本为 `run_q2_constant_fixed_reference.m`。它从0.02 A/cm²起，2秒线性升至0.259937538852887 A/cm²，此后保持平台电流。gammaIce=3.5、kFreeze=0.3、初温-10 ℃、粗网格[2 3 4 4 2]、原全时域单次ode15s。独立复核启动时间68.2624274374444秒，最低单片电压0.332455690000669 V。该电流最初由前轮拉丁超立方初采样发现；现在作为固定运行的对照组，但不能称作未经过任何参数选择的原始0.27 A/cm²加载。

运行 `report=run_q2_constant_fixed_reference_optimization()`，进行30次真实模型评价。平台电流全局搜索范围[0.02,0.5] A/cm²，全部候选及模型电流上限均不超过0.5 A/cm²。初始12个样本在[0.245,0.270] A/cm²内采用拉丁超立方采样，其中第一个是固定参考点。完成初始采样后，FuRBO检查点排序与信赖域决定局部搜索范围，解析稳定LogCEI从候选点选取下一次真实评价。方法遵循Word中的FuRBO＋LogCEI组合，非FuRBO原文Thompson采样的完整复现。

参考组与优化组只改变平台电流；物理参数、初始电流、2秒升流、粗网格、积分、精度及约束保持一致。每次优化评价设置300秒实际计算时间保护；超时为数值未完成，不是物理失败。旧轮搜索样本不会用于这轮代理模型。报告保存在`optimization_results/constant_fixed_reference_kfreeze03/`，结束后独立重跑最佳候选。若需续跑，可调用`run_q2_constant_fixed_reference_optimization(40,'上一次完整报告的路径')`，把总预算提高到40次。
