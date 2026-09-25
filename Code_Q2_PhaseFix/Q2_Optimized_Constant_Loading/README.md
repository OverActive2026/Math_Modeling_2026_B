# 问题 2：线性升流后恒流加载优化

此文件夹可单独复制到 GitHub 仓库并在 MATLAB 中运行。加载方式是从 **0.02 A/cm²** 起，在 **2 s** 内线性升至目标电流，然后保持该电流。当前只优化目标平台电流；这不是从零时刻就施加目标值的瞬时恒流实验。

## 已验证结果

| 指标 | 固定参考组 | 优化组 |
|---|---:|---:|
| 平台电流 / A·cm⁻² | 0.259937538852887 | 0.266880844760886 |
| 启动时间 / s | 68.2624274374444 | 65.0058211159030 |
| 最低单片电压 / V | 0.3324556900 | 0.3240780645 |

在相同模型及约束下，本次优化使启动时间缩短 **3.2566063215 s（4.7707156686%）**。优化组是 30 次真实仿真评估中第 28 次找到的最佳成功候选，并已独立重跑验证；不代表全局最优。30 次中有 18 次成功、12 次触发物理失败条件，没有数值未完成。固定参考平台电流最初由前轮拉丁超立方采样发现，随后被固定并重新仿真验证；论文中应如实说明其来源。

## 如何运行

要求 MATLAB 与 Statistics and Machine Learning Toolbox（搜索使用 `lhsdesign` 和 `fitrgp`）。将整个文件夹复制到任意位置，在 MATLAB 中进入该文件夹：

```matlab
cd('你的路径/Q2_Optimized_Constant_Loading')
run_q2_constant_fixed_reference                 % 单独运行固定参考组
result = run_q2_constant_best_current();         % 重跑已验证的最佳电流并画对比图
report = run_q2_constant_fixed_reference_optimization(); % 重新进行 30 次搜索，耗时较长
```

仿真运行时会在命令窗口持续输出模拟时间和实际计算时间。前两个入口使用已经确定的电流，最后一个入口才会进行完整的寻优。重跑搜索时若要延长预算，可调用 `run_q2_constant_fixed_reference_optimization(40, '先前完整报告的路径')`。

## 实验条件与算法

- `gammaIce = 3.5`，`kFreeze = 0.3`，初温 `−10 °C`。
- 电流密度在任意阶段均不得超过 **0.5 A/cm²**；最低单片电压 **0.30 V**；电荷上限 **20 C/cm²**。其他终止条件见模型代码。
- 仅采用粗网格 `[2 3 4 4 2]`；求解使用原始单次 `ode15s` 积分。此结果尚未用完整网格验证。
- 首 12 个平台电流在 `[0.245, 0.270] A/cm²` 内进行拉丁超立方初采样，其中包含固定参考点；后续按文档所述 FuRBO 检查点／信赖域与 LogCEI 自适应选点。候选平台电流全局范围为 `[0.02, 0.5] A/cm²`。这是 Word 方案中的组合实现，不能称作两篇文献算法的逐字完整复现。
- 两组仅改变平台电流，起始电流、2 s 升流、物理参数、网格、积分设置和约束相同。

## 文件

- `run_q2_constant_fixed_reference.m`：固定参考组入口。
- `run_q2_constant_best_current.m`：读取已验证报告并复跑优化组、绘图。
- `run_q2_constant_fixed_reference_optimization.m`：从头进行寻优。
- `q2_*.m`、`pemfc_*.m`、`*_state_ice.m`：搜索、加载与模型依赖，运行时须保留在同一文件夹。
- `optimization_results/q2_constant_fixed_reference_gamma3p5_kfreeze03.mat`：固定参考轨迹。
- `optimization_results/constant_fixed_reference_kfreeze03/`：30 次评估的完整 `.mat` 报告、`comparison_30eval.json` 精确数值和 `reference_vs_optimized_30eval.png` 对比图。
- [详细实验说明](README_恒流固定参考组与文献算法.md)：更多算法和对照组说明。

历史报告的 `.mat` 内有生成时的原始绝对路径元数据；新文件夹中的运行入口使用自身所在目录查找代码和报告。
