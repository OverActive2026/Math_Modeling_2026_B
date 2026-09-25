# 第一版含冰模型说明

本目录是在 `Code_Simplified` 基础上独立建立的含冰版本。所有函数均增加 `_ice` 后缀，避免 MATLAB 路径中同时存在两个版本时误调用。

## 文件结构

- `pemfc_setup_ice.m`：网格、状态编号和初值。
- `water_ice_state_ice.m`：总水守恒、cCL离聚物水合、冰守恒、相分配和水输运。
- `thermal_temperature_state_ice.m`：电压、冰覆盖修正、反应热、相变潜热和温度。
- `gas_transport_state_ice.m`：含冰动态孔隙率下的氢气和氧气输运。
- `pemfc_calculate_ice.m`：读取附件2的电流参数、组装RHS并调用 `ode15s`。
- `run_ice_model.m`：集中调整水合参数和三个冰参数，并调用-20/-25 degC工况。
- `run_ice_model_voltage_losses.m`：调用脚本副本，新增Erev、三项电压损失、PEM及cCL含水量lambda随时间的绘图。

## 修改标记

- `[ICE-1]`：初始冰。吹扫基准下，多孔层 `mi0=0`；PEM仍保留 `lambda0=3` 的结合水。
- `[ICE-2]`：冰守恒。`dmi/dt=Rfreeze-Rmelt`，相变不进入总水源项。
- `[ICE-3]`：孔隙率。`eps_i=mi/rho_i`，气相孔隙率直接扣除冰和液态水体积。
- `[ICE-4]`：有效反应面积。只在阴极催化层用冰饱和度修正交换电流密度。
- `[ICE-5]`：相变潜热。冻结放热、融化吸热。
- `[FIX-1]`：CL/PEM界面吸附和解吸使用对称导通关系，移除0.001单向修正；当前两侧速率因子恢复为30。
- `[FIX-2]`：极限电流使用cCL+cGDL全路径串联扩散阻力。
- `[FIX-3]`：电渗拖曳覆盖aCL/PEM和PEM/cCL两个边界面。
- `[FIX-5]`：欧姆损失加入两个CL离聚物内的质子传导损失；使用附件1给出的CL离聚物含量0.3和CL内已有的质子电流分布，按 `etaOhmCL=int(im/kappaCL)dx`计算。
- `[HYD-1]`：增加一个cCL平均离聚物含水量状态。cCL总水明确拆为`mIon+mv+ml+mi`，反应产水优先水合离聚物，结合水不能冻结。
- `[HYD-2]`：cCL水合上限采用附件所示的低温最大非冻结含水量`lambdaSat(T)`；达到上限后的产水才进入孔隙液水、排水和结冰。
- `[WATER-1]`：水蒸气与液水分开输运。水蒸气使用Fick扩散，液水使用Darcy/Leverett毛细通量，不再把`mv+ml`整体套用水蒸气扩散系数。

PEM没有设置孔隙冰状态，因为当前PEM水采用结合水含量 `lambda` 描述，不宜直接套用GDL/CL的孔隙冰模型。第一版对四个多孔层采用相同相变速率系数，但各层原始孔隙率不同，因此相同冰质量对应的冰饱和度不同。

## 水合与冰参数

```matlab
opt.init.lambdaCCL0 = 1.0; % 吹扫后cCL初始离聚物含水量
opt.tauHyd = 10;           % cCL水合时间常数 [s]
opt.kFreeze = 1;           % 冻结非平衡速率系数 [1/s]
opt.kMelt = 0;             % 当前低温数据不能辨识融化速率
opt.gammaIce = 0;          % 电压约束拟合的下界结果
```

低温最大非冻结含水量采用：`T<223.15 K`时`lambdaSat=4.837`；
`223.15<=T<273.15 K`时
`lambdaSat=(-1.304+0.01479*T-3.594e-5*T^2)^(-1)`。
该上限不是新的拟合参数。

多孔层相变采用分段关系：低于或等于冰点时
`Rfreeze=kFreeze*rhoLiquid*eps0*epsLiquid`，高于冰点时
`Rmelt=kMelt*rhoIce*epsIce`；不再显式乘过冷度、过热度或剩余孔隙因子。

冰孔隙率不是独立拟合参数，而由 `eps_i=mi/920` 直接计算。若后续实验显示冰形貌造成的有效堵塞体积偏离几何体积，再增加单独的堵塞修正系数。

## 活化参数标定

使用-20 degC和-25 degC两个零时刻实验电压联合反算得到：

```matlab
j0Ref = 0.00496453590833; % [A/m^2]，参考温度298.15 K
Ea = 7028.980902;         % [J/mol]
```

上述数值是在`lambdaCCL0=1`并加入`[HYD-1]`后，重新利用-20/-25 degC零时刻电压联合反算的结果。选择零时刻是为了避开冰累积及后续水传输结构对活化参数的干扰。代码中的修改位置使用`[ACT-1]`标记。

## 快速升流阶段拟合诊断

- 6-13 s内，-20/-25 degC实验电压分别下降约0.316/0.318 V；未计CL质子电阻时，模型只下降约0.177/0.173 V。
- 该阶段模型中的cCL冰饱和度基本为0，将界面速率因子由30降到1、或调整 `gammaIce`，对早期压降几乎没有影响；因此不能用冰参数拟合这一段。
- 加入独立cCL水合状态后，`lambdaCCL0=1、tauHyd=10 s`时，两个温度工况6-13 s区间的联合电压RMSE约0.019 V；初始CL质子电阻随产水水合自然衰减，不再使用人为时间开关。
- `tauHyd`按快速升流区间优先选择。若只优化36.6 s全区间会得到更小的1-2 s，但会明显削弱6-13 s的快速压降，因此没有采用。
- 在守恒水量和毛细排水同时存在时，电压拟合把`gammaIce`推到非负约束下界0；`kMelt`因温度未越过冰点不可辨识。`kFreeze=1 1/s`保留为文献量级，而不是宣称已由电压唯一标定。
- 当前模型在20 s以后仍普遍低估电压恢复，这部分不能再用增大结冰惩罚修正；后续应优先检查热源/热容、接触电阻随温度变化或实验电流边界解释。

## 运行方法

直接打开 `run_ice_model.m`，可修改脚本顶部的`lambdaCCL0`、`tauHyd`、`j0Ref`、`Ea`、`kFreeze`、`kMelt`和`gammaIce`后运行。脚本输出保存在`results.T20`和`results.T25`。

```matlab
addpath('Code_ICE');

result20 = pemfc_calculate_ice(-20);
result25 = pemfc_calculate_ice(-25,struct('plot',true));

opt = struct('tauHyd',10,'kFreeze',1,'kMelt',0, ...
             'gammaIce',0,'tEnd',10,'plot',true, ...
             'init',struct('lambdaCCL0',1));
result20Tune = pemfc_calculate_ice(-20,opt);
```

主要冰输出位于：

- `result.iceMassPerArea`
- `result.maxIceSaturation`
- `result.iceSaturationCCL`
- `result.iceAreaFactor`
- `result.latentHeatPerArea`
- `result.lambdaCCL`
- `result.lambdaCCLEquilibrium`
- `result.lambdaSaturationCCL`
- `result.cCLIonomerWaterMassPerArea`
- `result.cCLLiquidWaterMassPerArea`

当 `kFreeze=0`、初始冰为0时，冰动力学被关闭，便于逐项调试新增机理。
