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
- `[HYD-3]`：增加有界的cCL水合有效反应面积。`lambdaCCL`仍独立影响质子欧姆电阻；新系数只表示干态催化表面接入连续质子通道的比例，充分水合后最多恢复到1。
- `[WATER-1]`：水蒸气与液水分开输运。水蒸气使用Fick扩散，液水使用Darcy/Leverett毛细通量，不再把`mv+ml`整体套用水蒸气扩散系数。
- `[WATER-2]`：阴极水蒸气边界由`mv=0`理想水汇改为`Nout=hVaporCathode*(mv_surface-mv_inlet)`有限传质；附件1给出干空气入口，因此`mv_inlet=0`。
- `[ICE-3a]`：正式区分题目式(9)的冰体积分数`epsIce=mi/rhoIce`与孔隙冰饱和度`sIce=epsIce/eps0`。

PEM没有设置孔隙冰状态，因为当前PEM水采用结合水含量 `lambda` 描述，不宜直接套用GDL/CL的孔隙冰模型。第一版对四个多孔层采用相同相变速率系数，但各层原始孔隙率不同，因此相同冰质量对应的冰饱和度不同。

## 水合与冰参数

```matlab
opt.init.lambdaCCL0 = 1.0; % 吹扫后cCL初始离聚物含水量
opt.tauHyd = 10;           % cCL水合时间常数 [s]
opt.hVaporCathode = 0.01;  % cGDL/阴极气道水蒸气传质系数 [m/s]
opt.fHydDry = 0.06;        % 干态cCL可利用反应面积比例
opt.lambdaHydOn = 3.5;     % 反应面积开始恢复的lambda
opt.lambdaHydWet = 8.0;    % 反应面积完全恢复的lambda
opt.nHyd = 3.0;            % 水合面积恢复指数
opt.kFreeze = 1;           % 冻结非平衡速率系数 [1/s]
opt.kMelt = 0;             % 当前低温数据不能辨识融化速率
opt.gammaIce = 0.1;        % 保留冰覆盖弱正耦合
```

低温最大非冻结含水量采用：`T<223.15 K`时`lambdaSat=4.837`；
`223.15<=T<273.15 K`时
`lambdaSat=(-1.304+0.01479*T-3.594e-5*T^2)^(-1)`。
该上限不是新的拟合参数。

多孔层相变采用分段关系：低于或等于冰点时
`Rfreeze=kFreeze*rhoLiquid*eps0*epsLiquid`，高于冰点时
`Rmelt=kMelt*rhoIce*epsIce`；不再显式乘过冷度、过热度或剩余孔隙因子。

冰孔隙率不是独立拟合参数，而由 `eps_i=mi/920` 直接计算。若后续实验显示冰形貌造成的有效堵塞体积偏离几何体积，再增加单独的堵塞修正系数。

`hVaporCathode`不是相变系数，而是综合cGDL表面至气道的未解析传质阻力。附件1未给出气道尺寸、流量和停留时间，所以暂不虚构0D气道状态；待有相应输入后可将本边界换成完整气道守恒。

有限排水灵敏度扫描中，`hVaporCathode=0.1 m/s`时-20 degC仍不结冰；降至`0.03/0.01/0.003 m/s`时，最大冰体积分数分别约为`0.004/0.035/0.050`。当前取`0.01 m/s`作为能保留有限排水且不使水完全封闭的中间值；没有冰量实验前，该参数不应声称已唯一辨识。

## 活化参数与水合反应面积标定

原模型只让`lambdaCCL`改善CL质子电阻，但不改善催化剂的质子可达面积。因此它能拟合快速升流段，却会在14 s以后持续低估电压恢复。新模型采用：

`thetaHyd=clip((lambdaCCL-lambdaHydOn)/(lambdaHydWet-lambdaHydOn),0,1)`

`fHyd=fHydDry+(1-fHydDry)*thetaHyd^nHyd`

`j0Eff=j0Wet(T)*fHyd*fIce`

参数标定遵循以下约束：

1. `lambdaHydOn=3.5`由实验电压最低点附近的`lambdaCCL`确定，保证水合面积恢复不会抬高前期压降。
2. `lambdaHydWet=8`用于表示cCL进入接近充分水合状态。
3. `fHydDry=0.06`表示初始干态只有6%的湿态有效反应面积。
4. `nHyd=3`使反应面积在电压最低点之后加速恢复，而非一启动就立即增大。
5. 为保持两个零时刻电压不变，使`j0RefWet*fHydDry`等于旧模型的干态有效`j0Ref`。

重新标定后：

```matlab
j0Ref = 0.0827422651388;  % [A/m^2]，298.15 K充分水合参考值
Ea = 7028.980902;         % [J/mol]
fHydDry = 0.06;
lambdaHydOn = 3.5;
lambdaHydWet = 8.0;
nHyd = 3.0;
```

干态零时刻有效参考值仍为`0.0827422651388*0.06=0.00496453590833 A/m^2`，因此原来的零时刻标定被完整保留。修改位置使用`[ACT-2]`和`[HYD-3]`标记。

加入有限排水边界并用`MaxStep=0.01 s`完整求解后，两工况结果为：

- -20 degC：电压RMSE为`0.01443 V`，温度RMSE为`0.05432 degC`，末端电压为`0.67657 V`。约`31.6 s`开始结冰，35 s最大冰体积分数为`0.01961`，全时段最大值为`0.03492`。
- -25 degC：电压RMSE为`0.01284 V`，温度RMSE为`0.13406 degC`，末端电压为`0.63958 V`。约`24.6 s`开始结冰，35 s最大冰体积分数为`0.10615`，全时段最大值为`0.12519`。
- 结冰位置均位于cCL，不是cGDL外边界的数值堆积。-20/-25 degC总水守恒误差分别为`-1.40e-6/-8.48e-7 kg/m^2`，其他健康检查均通过。

## 快速升流阶段拟合诊断

- 6-13 s内，-20/-25 degC实验电压分别下降约0.316/0.318 V；未计CL质子电阻时，模型只下降约0.177/0.173 V。
- 该阶段模型中的cCL冰饱和度基本为0，将界面速率因子由30降到1、或调整 `gammaIce`，对早期压降几乎没有影响；因此不能用冰参数拟合这一段。
- 加入独立cCL水合状态后，`lambdaCCL0=1、tauHyd=10 s`时，两个温度工况6-13 s区间的联合电压RMSE约0.019 V；初始CL质子电阻随产水水合自然衰减，不再使用人为时间开关。
- `tauHyd`按快速升流区间优先选择。若只优化36.6 s全区间会得到更小的1-2 s，但会明显削弱6-13 s的快速压降，因此没有采用。
- 加入有限排水后，电压单目标扫描仍把`gammaIce`推向非负下界0，说明它无法只由当前两条电压曲线唯一辨识。附件1给出的参考指数为3.5，但直接使用会将当前-25 degC末段冰影响放大过度；因此这一版暂取`gammaIce=0.1`保留弱正耦合，后续应用冰体积分数实验或冻结形貌数据标定。
- `kMelt`因温度未越过冰点不可辨识；`kFreeze=1 1/s`保留为文献量级，而不是宣称已由电压唯一标定。
- `[HYD-3]`已将原先的14 s后系统性低估大幅消除。当前最大残差转移到10.2 s左右，后续应优先检查升流阶段的局部活化/氧浓度耦合，而不是再增大冰惩罚。

## 运行方法

直接打开 `run_ice_model.m`，可修改脚本顶部的`lambdaCCL0`、`tauHyd`、`hVaporCathode`、`j0Ref`、`Ea`、`kFreeze`、`kMelt`和`gammaIce`后运行。脚本输出保存在`results.T20`和`results.T25`。

```matlab
addpath('Code_ICE');

result20 = pemfc_calculate_ice(-20);
result25 = pemfc_calculate_ice(-25,struct('plot',true));

opt = struct('tauHyd',10,'fHydDry',0.06, ...
             'lambdaHydOn',3.5,'lambdaHydWet',8,'nHyd',3, ...
             'hVaporCathode',0.01, ...
             'kFreeze',1,'kMelt',0,'gammaIce',0.1, ...
             'tEnd',10,'plot',true, ...
             'init',struct('lambdaCCL0',1));
result20Tune = pemfc_calculate_ice(-20,opt);
```

主要冰输出位于：

- `result.iceMassPerArea`
- `result.maxIceVolumeFraction`：题目式(9)，各时刻全域最大冰体积分数
- `result.iceVolumeFractionCCL`：cCL厚度平均冰体积分数
- `result.validationTable`：题目表1/表2所需的0、5、...、35 s电压、温度、相对误差和最大冰体积分数
- `result.waterBalance`：总水守恒的实际变化、理论变化和误差
- `result.maxIceSaturation`
- `result.iceSaturationCCL`
- `result.iceAreaFactor`
- `result.latentHeatPerArea`
- `result.lambdaCCL`
- `result.lambdaCCLEquilibrium`
- `result.lambdaSaturationCCL`
- `result.hydrationDegree`
- `result.hydrationAreaFactor`
- `result.cCLIonomerWaterMassPerArea`
- `result.cCLLiquidWaterMassPerArea`

当 `kFreeze=0`、初始冰为0时，冰动力学被关闭，便于逐项调试新增机理。
