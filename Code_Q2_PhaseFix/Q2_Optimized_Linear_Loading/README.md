# 问题2：优化后线性加载（独立 MATLAB 代码）

本文件夹可单独放进 GitHub 仓库，不依赖上一级目录的 `.m` 文件。`run_q2_linear_comparison.m` 在同一快速网格和同一物理参数下，分别重跑未优化与已优化的五片电堆冷启动，并生成对比图；无需历史检查点。

## 给队友的最短运行方法

在 MATLAB 中把“当前文件夹”设为本目录，执行：

```matlab
verification = run_q2_linear_comparison();
```

仓库内已附上一份验证过的 `q2_linear_comparison_result.mat` 和 `q2_linear_comparison.png`；重新运行会在本地更新它们。命令窗口显示两条曲线的**仿真启动时间**、缩短比例以及约束检查；这不是电脑运行耗时。此入口会实算两次五片模型，耗时取决于电脑。依赖检查识别到 MATLAB、System Identification Toolbox、Statistics and Machine Learning Toolbox；优化器明确使用后者的 `fitrgp` 与 `lhsdesign`。

## 策略及当前结果

两种策略均为 `j(t)=min(j_0+r t,j_p)`，单位分别为 A/cm²、A/(cm²·s)、A/cm²。

| 策略 | `j_0` | `r` | `j_p` | 启动时间 |
|---|---:|---:|---:|---:|
| 未优化（原 `run_q2_linear_current.m`） | 0.020000000 | 0.010000000 | 0.400000000 | 48.838844519 s |
| 优化后（真实仿真得到） | 0.117641245 | 0.008079398 | 0.447158213 | 40.726546388 s |

相同设置下缩短 8.112298131 s（16.61034%）。优化方案的最小单片电压为 0.311869712 V，启动时累计电荷为 11.491596076 C/cm²，启动时最大局部冰体积分数为 0.197513767。将求解相对、绝对容差再收紧一倍后，该方案仍启动成功，得到 40.724230793 s。

以上是**五片模型快速网格** `[2 3 4 4 2]`、初始温度 −10 °C、`gammaIce=3.5`、`kFreeze=0.4 s^-1`、`RelTol=5e-5`、`AbsTol=5e-9` 的模型计算结果；尚未完成完整网格验证，不能写成实验实测或全局最优。原脚本单独运行时容差较宽，其末位数字可能与本公平对比略有不同。

## 代码分工

- `run_q2_linear_comparison.m`：无需缓存的论文对比复现入口，运行两次真实模型并画图。
- `run_q2_linear_current.m`：原始未优化线性加载的独立演示脚本；不建议把它的较宽容差输出直接与严格容差输出比较。
- `run_q2_optimized_expanded.m`：重新启动优化搜索，每运行一次会增加 **6 次新的昂贵物理评估**；可能得到不同结果，不保证再次找到上述最优值。首次在新文件夹中运行没有旧检查点，会从头搜索。
- `q2_optimize_furbo_real.m`、`q2_optimize_trbo.m`、`q2_logei.m`：优化算法和真实模型接口。
- `q2_decode_strategy.m`、`q2_current_strategy.m`：归一化参数与加载曲线。
- `pemfc_stack5_simulate.m`、`pemfc_setup_ice.m`、`*_state_ice.m`：五片电堆物理模型。
- `verify_q2_literature_best.m`：如果本目录已有优化器生成的 `q2_furbo_expanded_strict_latest_report.mat`，可对该报告的当前最优值另做无截断验证；否则请使用上面的对比入口。

## 算法口径与论文表述

优化变量为 `(j_0,r,j_p)`，目标是**在真实五片模型全部冷启动约束满足时，最小化启动仿真时间**。失败方案没有有效启动时间，只提供可行性标签；成功方案才用于拟合时间代理模型。代码使用受 [FuRBO：Feasibility-Driven Trust Region Bayesian Optimization](https://arxiv.org/abs/2506.14619) 启发的检查点/信赖域候选搜索，并使用 [Ament 等的 LogEI](https://arxiv.org/abs/2310.20708) 对可行候选排序。最终候选仍由真实五片物理仿真确认。这里是组合式、经简化的实现，**不能称为 FuRBO 原论文的严格复现**。

约束包括 `j≤0.50 A/cm²`、最小单片电压不低于 0.30 V、累计电荷不超过 20 C/cm²、局部冰体积分数低于 0.99，以及模型中的孔隙占用上限 0.98。当前进一步提速主要受电压下限限制。
