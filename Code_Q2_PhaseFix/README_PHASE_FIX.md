# Q2水-冰相变修正版

本目录从`Code_Q2`完整复制而来，原目录未修改。本版专门处理
“-20 degC几乎无冰”以及五片恒流中“不结冰却启动失败”的机理问题。

## 主要修改

1. `[PHASE-1]` CL/PEM界面两侧统一使用离聚物归一化水活度，
   不再用PEM的lambda与cCL孔隙蒸气活度直接比较。界面速率倍率
   `interfaceRateFactor`显式化，保守默认值为`0.01`。
   `[PHASE-1b]`在此基础上使实际导通率随界面水合程度对称增强：
   `factor=interfaceRateFactor*(1+29*aWet^4)`。它在干态仍接近0.01，
   充分水合时最高约0.3，且正反两个方向使用同一关系，不是单向阀。
2. `[PHASE-2]` 阳极和阴极外边界都改为对干气道的有限水蒸气传质，
   默认`hVaporAnode=hVaporCathode=0.01 m/s`，不再把阳极当作无限快零浓度水汇。
3. `[PHASE-3]` cCL离聚物水合速率受实际孔隙水和反应产水速率限制。
   最多只有`1-directFreezeFraction`的当时产水可用于水合，避免所有产水
   被lambda状态瞬时吸收。
4. `[PHASE-4]` 增加低温产水直接冻结通道。它不增加总水，只在液水/冰之间
   重分配；受`freezeNucleationLambdaFraction`和局部可移动水成核因子约束。
5. `[PHASE-5]` 保留分段冻结/融化式，并给出正的`kMelt=1 1/s`默认值。
   为避免刚性求解中的非物理Newton试探直接中断，增加了温度、反应气体和
   气相孔隙率的本构计算投影；接受解仍由守恒与约束检查判定。
6. `[NUM-PHASE]` 问题2在0 degC附近使用可调的窄平滑相变开关，默认
   `phaseTransitionWidthK=0.2 K`。它只替代冰点邻域内的0/1硬切换，远离
   冰点仍退化为原冻结/融化关系；设为0可严格恢复题目分段式。问题1的
   `pemfc_calculate_ice`默认仍为0，不改变原拟合口径。
7. `[NUM-JAC]` 五片模型默认向`ode15s`提供保守稀疏雅可比结构。单片内部
   使用安全的稠密块，片间只保留相邻导热和端板耦合；这不改变物理方程，
   只减少数值雅可比计算时不必要的整套RHS重复调用。可用
   `useJacobianPattern=false`关闭。
8. 若`ode15s`没有触发任何物理事件却提前返回，结果现在会明确标记
   `solverTerminatedUnexpectedly=true`，终止原因不再误报为“达到最大时间”。

## 当前保守默认值

```matlab
simOpt.kFreeze = 0.05;                    % 1/s
simOpt.kMelt = 1.0;                       % 1/s
simOpt.phaseTransitionWidthK = 0.2;       % K，仅问题2默认启用
simOpt.useJacobianPattern = true;
simOpt.directFreezeFraction = 0.01;       % -
simOpt.freezeNucleationLambdaFraction = 0.5;
simOpt.interfaceRateFactor = 0.01;
simOpt.hVaporAnode = 0.01;                % m/s
simOpt.hVaporCathode = 0.01;              % m/s
```

`kFreeze=1 1/s`与`directFreezeFraction=0.1`在当前体积平均模型中会导致
约25 s时最大冰体积分数超过0.30，并引起过强压降，因此不再作默认值。

## 已做的数值检查

- -20 degC短时积分稳定，冰体积分数不再恒为0。
- -20 degC在正式容差下可稳定跨过10 s；较大冻结系数会显著加重后段刚性。
- 五片堆、-10 degC、0.2 A/cm^2的粗网格检查中，10 s时最低历史电压
  约0.380 V，末端电压约0.589 V，最高单片温度约-5.28 degC，
  最大冰体积分数约0.0122。
- 若界面倍率始终固定为0.01，0.10--0.25 A/cm^2都会在约2.4--2.75 C/cm^2时
  触发0.30 V下限；当时冰体积分数仅约0.015--0.027。
- 加入`[PHASE-1b]`后，0.2 A/cm^2粗网格工况在40 s仍未违反电压约束：
  历史最低电压约0.380 V，40 s电压约0.714 V，PEM平均lambda恢复至约8.04，
  最高温度约-1.48 degC，最大冰体积分数约0.0546。55 s时最热片约
  -0.26 degC，但最冷片仍约-3.27 degC，因此还未满足“五片同时0 degC”。

固定0.01界面倍率下的进一步诊断表明：最危险的第3片极限电流仍约
25.2 A/cm^2，浓差损失近似0，最大液水饱和度约0.063；同时活化损失约
0.488 V、总欧姆损失约0.457 V，其中PEM/aCL/cCL约为0.260/0.169/0.026 V；
此时PEM平均lambda仅约0.96，cCL lambda却约9.44。`[PHASE-1b]`正是针对这种
“cCL很湿、PEM和阳极侧很干”的后期失配。修正后的剩余瓶颈主要是五片间
热分布不均以及接近0 degC时的相变刚性，不应继续用`kFreeze`吸收电压误差。

### 57 s附近停滞的数值诊断

在-10 degC、0.2 A/cm^2、`nCellLayer=[2 3 4 4 2]`、`MaxStep=0.1 s`的
对照试验中：

- 保留稀疏雅可比但令`phaseTransitionWidthK=0`时，最热单片接近0 degC后，
  `ode15s`在`t=58.09995 s`报步长必须小于允许最小值，只返回到58.0 s；
- 令`phaseTransitionWidthK=0.2 K`后，模型连续跑到65 s，耗时约15.4 s；
  末端单片平均温度范围为-2.577至0.206 degC，最大局部冰体积分数约0.0940，
  最低单片末端电压约0.718 V；
- 因此原停滞的直接触发因素是部分节点越过冰点时冻结/融化硬开关反复切换，
  而不是电荷上限、冰上限或单纯的网格规模。稀疏雅可比主要降低每次隐式
  Newton迭代的代价，窄平滑层负责消除58 s附近的最小步长失败。

## 新增诊断量

`pemfc_stack5_simulate`的返回值新增：

- `limitingCurrentDensityAcm2`
- `maximumLiquidSaturationCell`
- `etaActCell`
- `etaOhmCell`
- `etaOhmPEMCell`
- `etaOhmACLCell`
- `etaOhmCCLCell`
- `etaConCell`
- `lambdaPEMMean`

这些量可直接判断下一次低电压是来自液水淹没、CL质子电阻，还是冰堵。
