# 第一版含冰模型说明

本目录是在 `Code_Simplified` 基础上独立建立的含冰版本。所有函数均增加 `_ice` 后缀，避免 MATLAB 路径中同时存在两个版本时误调用。

## 文件结构

- `pemfc_setup_ice.m`：网格、状态编号和初值。
- `water_ice_state_ice.m`：总水守恒、冰守恒、三相分配和水输运。
- `thermal_temperature_state_ice.m`：电压、冰覆盖修正、反应热、相变潜热和温度。
- `gas_transport_state_ice.m`：含冰动态孔隙率下的氢气和氧气输运。
- `pemfc_calculate_ice.m`：读取附件2的电流参数、组装RHS并调用 `ode15s`。
- `run_ice_model.m`：集中调整三个标定参数并调用-20/-25 degC工况。
- `run_ice_model_voltage_losses.m`：调用脚本副本，新增Erev、三项电压损失及膜含水量lambda随时间的绘图。

## 修改标记

- `[ICE-1]`：初始冰。吹扫基准下，多孔层 `mi0=0`；PEM仍保留 `lambda0=3` 的结合水。
- `[ICE-2]`：冰守恒。`dmi/dt=Rfreeze-Rmelt`，相变不进入总水源项。
- `[ICE-3]`：孔隙率。`eps_i=mi/rho_i`，气相孔隙率直接扣除冰和液态水体积。
- `[ICE-4]`：有效反应面积。只在阴极催化层用冰饱和度修正交换电流密度。
- `[ICE-5]`：相变潜热。冻结放热、融化吸热。
- `[FIX-1]`：CL/PEM界面吸附和解吸使用对称导通关系，移除0.001单向修正。
- `[FIX-2]`：极限电流使用cCL+cGDL全路径串联扩散阻力。
- `[FIX-3]`：电渗拖曳覆盖aCL/PEM和PEM/cCL两个边界面。

PEM没有设置孔隙冰状态，因为当前PEM水采用结合水含量 `lambda` 描述，不宜直接套用GDL/CL的孔隙冰模型。第一版对四个多孔层采用相同相变速率系数，但各层原始孔隙率不同，因此相同冰质量对应的冰饱和度不同。

## 三个待标定参数

```matlab
opt.kFreeze = 0.01;  % 冻结系数 [1/(K s)]，当前仅为初值
opt.kMelt = 0.01;    % 融化系数 [1/(K s)]，当前仅为初值
opt.gammaIce = 3.5;  % cCL冰覆盖修正指数，当前仅为初值
```

冰孔隙率不是独立拟合参数，而由 `eps_i=mi/920` 直接计算。若后续实验显示冰形貌造成的有效堵塞体积偏离几何体积，再增加单独的堵塞修正系数。

## 活化参数标定

使用-20 degC和-25 degC两个零时刻实验电压联合反算得到：

```matlab
j0Ref = 0.00807579349527; % [A/m^2]，参考温度298.15 K
Ea = 22204.0159133;       % [J/mol]
```

选择零时刻是为了避开冰累积及后续水传输结构对活化参数的干扰。代码中的修改位置使用 `[ACT-1]` 标记。

## 运行方法

直接打开 `run_ice_model.m`，可修改脚本顶部的 `j0Ref`、`Ea`、`kFreeze`、`kMelt` 和 `gammaIce` 后运行。脚本输出保存在 `results.T20` 和 `results.T25`。

```matlab
addpath('Code_ICE');

result20 = pemfc_calculate_ice(-20);
result25 = pemfc_calculate_ice(-25,struct('plot',true));

opt = struct('kFreeze',0.02,'kMelt',0.01, ...
             'gammaIce',4.0,'tEnd',10,'plot',true);
result20Tune = pemfc_calculate_ice(-20,opt);
```

主要冰输出位于：

- `result.iceMassPerArea`
- `result.maxIceSaturation`
- `result.iceSaturationCCL`
- `result.iceAreaFactor`
- `result.latentHeatPerArea`

当 `kFreeze=0`、初始冰为0时，冰动力学被关闭，便于逐项调试新增机理。
