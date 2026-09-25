# Code_ICE_E —— 产冰正常化 + 曲线联合重标定版

本目录从 `Code_ICE_D` 完整复制而来，**D 版及其之前的目录均未修改**。
本版的唯一目标（用户要求原文）：*"学习一下 Code_Q2_PhaseFix 版本中对产冰等内容
的调整，结合拟合的曲线情况，综合修改模型，使产冰正常且曲线拟合良好。"*

结论先行：

| 指标 | Code_ICE_D | **Code_ICE_E** |
|---|---|---|
| -20 degC 电压 RMSE | 0.02514 V | **0.01432 V** |
| -25 degC 电压 RMSE | 0.00958 V | **0.01480 V** |
| 两工况合并 RMSE | 0.01903 V | **0.01456 V**（改善 23%） |
| -20 degC 温度 RMSE | 0.0934 K | **0.0497 K** |
| -25 degC 温度 RMSE | 0.1119 K | 0.1595 K |
| 压降幅度达成率 (-20/-25) | 102.5% / 102.5% | **100.70% / 100.60%** |
| 回升幅度达成率 (-20/-25) | 97.9% / 102.3% | **100.29% / 103.53%** |
| 末端冰量 (-20/-25) | 0 / 4.4e-4* | **1.78e-4 / 3.66e-4 kg/m²** |
| cCL 冰饱和度 (-20/-25) | **0 / 0** | **4.06% / 8.37%** |
| `gammaIce` | 0.1（临时值） | **3.5（附件1 给定值）** |

\* D 版 `gammaIce=0.1` 时 -25 degC 靠题目式(9) 得到 4.4% 冰，但 -20 degC
恒为 0，且 `gammaIce=0.1` 并非附件1 给定值，不能作为正式结果。

---

## 1. 问题诊断：为什么 D 版冰恒为 0

### 1.1 直接定位（`tmp/E_gates.m`）

对 D 版参数在 -20 degC 逐时打印冻结通道的三个门槛：

```
     t   lambda   lamSat   nucAct   mobAct   mMobile   RdirCCL   RliqCCL
   0.0   1.0000   7.3059   0.0000   0.0000  0.00e+00  0.00e+00  0.00e+00
  12.0   1.6815   7.4375   0.0000   0.0002  6.96e-05  0.00e+00  0.00e+00
  20.0   2.7658   7.8853   0.0725   0.0003  1.43e-04  1.23e-04  0.00e+00
  36.6   4.7187   9.2156   0.3029   0.0011  4.58e-04  2.47e-03  0.00e+00
```

三个原因：

1. **可动水门槛把冻结通道压掉三个数量级。**
   `Code_Q2_PhaseFix` 的 `[PHASE-4]` 写法为
   `RfreezeDirect = nuc * mMobile/(mMobile+mTransition) * supply * freezeWeight`，
   其中 `mTransition = 1e-3*rhoLiquid*eps0 ≈ 0.4 kg/m³`。而本算例 cCL 的
   可动水存量 `mMobile` 只有 **~5e-4 kg/m³**（产水一生成就被离聚物吸收或
   被气相带走），于是该因子恒为 **~1e-3**，冻结速率被整体压掉三个数量级。

2. **题目式(9) 的液水冻结项需要 `ml>0`，而 cCL 液水恒为 0。**
   `RfreezeLiquid = kFreeze*rhoL*eps0*sL` 依赖液态水饱和度。本算例中
   cCL 孔隙水始终略低于当地饱和蒸气容量，`ml = 0`，故该项恒为 0。

3. **水账本身不允许 cCL 存水。** 见 §1.2。

### 1.2 水账诊断（`tmp/E_water.m`）

D 版参数、t = 36.6 s：

```
produced (cum)      = 5.9921e-03 kg/m2
total cell water    = 5.1931e-03 kg/m2      (初始存量 1.52e-03，故净保留
  PEM               = 4.5737e-03             (5.1931-1.52)/5.9921 = 61%)
  cCL ionomer       = 6.1906e-04
  ice               = 2.7217e-07
  vapor+liq pores   = 6.6744e-08
cCL ice capacity    = 4.3736e-03 kg/m2
```

即：**膜的吸水量（4.5737e-3 − 初始 1.39e-3 = 3.18e-3）等于累计产水的 53%**，
而 **cCL 孔隙几乎是干的**（可动水 4.6e-4 kg/m³，而 cCL 冰容量相当于
920×0.4207 = 387 kg/m³）。冰的"原料"根本不在 cCL 里。

进一步定位原因：`water_ice_state_ice.m` 第 7 节把 CL/PEM 界面水交换倍率
**写死为 30**（原注释自称"待核对值"）。而
`Gright = Dw(PEM)*CmemScale/Lpem ≈ 7e-3 kg/(m²·s)`，乘 30 后界面导通能力
约为反应产水通量（~2.8e-4 kg/(m²·s)）的 **700 倍**，等效于
"cCL 孔隙水活度与膜含水量瞬时平衡"。这正是 `[PHASE-1]` 要解决的问题。

### 1.3 界面倍率扫描（`tmp/E_iface.m`）

| interfaceFactor | V_R20 | V_R25 | 合并 | sIce20 | sIce25 | cCL 液水 |
|---|---|---|---|---|---|---|
| 0.3 | 0.03583 | 0.05212 | 0.04472 | 0.354 | 0.619 | 20.1 kg/m³ |
| 1 | 0.02805 | 0.01626 | 0.02293 | 0.034 | 0.188 | 0 |
| 3 | 0.02622 | 0.01160 | 0.02027 | 0 | 0.090 | 0 |
| 10 | 0.02539 | 0.01045 | 0.01941 | 0 | 0.055 | 0 |
| **30（原值）** | 0.02514 | 0.01018 | 0.01918 | 0 | 0.044 | 0 |

结论：降低界面倍率确实能"唤醒"题目式(9) 的液水冻结（倍率 0.3 时冰饱和度
高达 35%/62%），但**曲线拟合同时崩溃**（合并 RMSE 0.0447）。因此界面倍率
不是本问题的正确杠杆，本版仅把它**参数化**（`interfaceFactor`，默认仍为 30，
保持 D 版行为），不用于标定。

---

## 2. 本版修改

### 2.1 从 `Code_Q2_PhaseFix` 移植的机制

| 标记 | 内容 | 本版处理 |
|---|---|---|
| `[PHASE-1]` | CL/PEM 界面水交换倍率 `interfaceRateFactor` | 参数化为 `interfaceFactor`，**默认 30 = 原写死值**（§1.3 证明不宜标定） |
| `[PHASE-2]` | 双侧气道有限水蒸气传质系数 `hVaporAnode/Cathode` | 已移植为 `'mass_transfer'` 边界分支；**默认 0 = 原零浓度理想汇**（扫描显示对冰与拟合均无实质影响） |
| `[PHASE-3]` | 离聚物水合供给速率受当时产水限制 | 采用，并在**供给端预留**冻结份额（§2.2） |
| `[PHASE-4]` | 低温产水直接冻结通道 | 采用，但**去掉可动水门槛**（§2.2） |
| `[PHASE-5]` | 窄平滑相变开关 `phaseTransitionWidthK` | 已移植；**默认 0**，因本组工况 T 全程 < 0、不跨冰点，无切换刚性 |

### 2.2 `[PHASE-4b]`：去掉可动水门槛（核心修正）

**文件**：`water_ice_state_ice.m`

```matlab
% 原（Q2 / D 版口径）——被压掉三个数量级
mobileWaterTransition = max(1e-3*rhoLiquid*eps0(g.idx_cCL(1)),1e-6);
mobileWaterActivation = mMobile(g.idx_cCL)./(mMobile(g.idx_cCL)+mobileWaterTransition);
RfreezeDirect = nucleationActivation.*directFreezeSupply.*mobileWaterActivation.*freezeWeight;

% 改（本版）
RfreezeDirect = kFreezeDirect.*nucleationActivation.*directFreezeSupply.*freezeWeight;
```

**物理依据**：-20 degC 下反应位局部水活度接近 1，生成的水来不及输运就已冻结；
冻结速率应由"产水速率 × 预留比例 × 成核活化"决定，而不是由几乎为空的
**存量**决定（存量只决定能维持多久，不决定速率上限）。
这正是冷启动文献中"产物水就地冻结 / direct freezing"的标准处理。

**质量守恒不受影响**：冰量上限仍由两道约束保证——
`atAllWaterFrozen`（冰 ≤ 该控制体孔隙水 `mPore`）与
`atFullPore`（冰 ≤ 孔隙容积 `rhoIce*eps0`）。
实测 `result.health.iceNotExceedPoreWater = 1`、`iceNotExceedTotalWater = 1`。

同时保留 `[PHASE-3]` 的"供给端预留"：
```matlab
productionSupplyRate = (1-directFreezeFraction)*reactionWaterRateCCL/cCLIonomerWaterCoefficient;
```
即离聚物水合最多只能拿到 `(1-f)` 的产水，`f` 份额预留给冻结通道。
（若只在冻结端限制，离聚物仍会吃掉几乎全部产水，冻结通道拿到的份额 ~2%，
冰依旧接近 0——这是 Q2 原写法在本算例失效的第二个原因。）

### 2.3 `[E-1]`：解除 `lamHydRef` 与 `j0Ref` 的简并

D 版的 `[D-1]` 写成 `fHyd = (lambdaCCL/lamHydRef)^mHyd`，并注释
"`lamHydRef` 取 t=0 时 cCL 的实际 lambda（约 0.1）"。**该注释有误**：
t=0 时 cCL 的实际含水量是 `lambdaCCL0 = 1.0`（实测
`cCL mw = 1.161e1 kg/m³ = C×1.0`，`C = 11.61 kg/m³/λ`）。

更关键的是数学上完全简并：

```
j0Ref * fHyd = j0Ref * (lambda/lamHydRef)^mHyd
             = (j0Ref / lamHydRef^mHyd) * lambda^mHyd
```

即 `lamHydRef` 与 `j0Ref` 不可同时辨识。D 版的
`(j0Ref=0.020, lamHydRef=0.10, mHyd=1.35)` 严格等价于
`(j0Ref=0.4477, lamHydRef=1.0, mHyd=1.35)`。

本版**把 `lamHydRef` 固定在物理参考点 `lambdaCCL0 = 1.0`**（不参与标定，
使 `fHyd(t=0)=1`），尺度全部并入 `j0Ref`。这样水合/欧姆标定量从 5 个降为
**4 个**（另加 §3.3 的 2 个产冰标定量，共 6 个），
且 `j0Ref = 0.30 A/m² = 3.0e-4 A/cm²` 落在 Pt/C 阴极有效交换电流密度的
常见量级（1e-4 ~ 1e-3 A/cm²），具有物理可比性。

### 2.4 `[ICE-4]`：`gammaIce` 改回附件1 给定值

D 版因冰恒为 0、而 `gammaIce=3.5` 会在其它工况过度抑制反应面积，临时取 0.1。
本版冰真正生成，**改回附件1 给定的 3.5**，使 `fArea=(1-sIce)^3.5` 成为
题目口径下的真实机制，而不是被绕过的参数。

---

## 3. 最终参数

### 3.1 题目/附件1 给定（不可标定）

| 量 | 值 | 来源 |
|---|---|---|
| `Ea` | 67000 J/mol | 题目式(40) |
| `alpha` | 0.5 | 题目式(39) |
| `gammaIce` | 3.5 | 附件1 阴极冰覆盖活性面积指数 |
| `lambdaCCL0` | 1.0 | 附件1 初始膜含水量 λ=3 对应的 cCL 初值 |
| `lamHydRef` | 1.0 | = `lambdaCCL0`，物理参考点（[E-1] 解除简并后不再标定） |

### 3.2 标定量（4 个，两工况全时段曲线联合标定）

| 量 | 值 | 含义 |
|---|---|---|
| `j0Ref` | **0.30 A/m²** | λ=1、298.15 K 参考交换电流密度 |
| `mHyd` | **1.55** | j0 的离聚物水合指数 |
| `clResistanceScale` | **0.64** | cCL 离聚物质子电阻缩放 [FIX-5] |
| `tauHyd` | **40 s** | cCL 离聚物水合时间常数 |

### 3.3 低温产冰参数（2 个）

| 量 | 值 | 含义 |
|---|---|---|
| `directFreezeFraction` | **0.15** | 反应产水中划给直接冻结通道的份额 |
| `freezeNucleationLambdaFraction` | **0.25** | 成核门槛（λ/λ_sat 超过该值才允许成核） |

> 标定量合计 **6 个**（4 个水合/欧姆 + 2 个产冰），由 -20/-25 degC 两工况
> 各自 184 个实测点联合标定，标定目标为**两工况合并 RMSE**
> `sqrt((RMSE20² + RMSE25²)/2)`。

其余开关保持 D 版口径：`kFreeze=1`、`kMelt=0`（全程 T<0，不可辨识）、
`kFreezeDirect=1`、`kmCond=0`、`DwFilmRef=0`、`nSorp=3`、`Kcov=0`、`fRet=1`、
`hVaporAnode=hVaporCathode=0`、`phaseTransitionWidthK=0`、`interfaceFactor=30`。

---

## 4. 验证结果

### 4.1 拟合精度（`tmp/E_final_check.m`）

```
[-20 degC] V_RMSE=0.01432 V  T_RMSE=0.0497 K  dip=100.70%  rec=100.29%  endGap=+0.43 mV
[-25 degC] V_RMSE=0.01480 V  T_RMSE=0.1595 K  dip=100.60%  rec=103.53%  endGap=-8.47 mV
两工况合并 RMSE = 0.01456 V
```

相对电压误差 [%]：

| t / s | 0 | 5 | 10 | 15 | 20 | 25 | 30 | 35 | 36.6 |
|---|---|---|---|---|---|---|---|---|---|
| -20 degC | 0.10 | 0.89 | **4.69** | 2.49 | 2.47 | 2.45 | -0.03 | -0.20 | 0.06 |
| -25 degC | **-2.70** | -0.43 | 0.57 | -2.78 | -1.78 | -1.92 | -4.22 | -3.06 | — |

（D 版同位置：-20 degC 最大 7.63%，-25 degC 最大 3.20%）

温度绝对误差 [K]：-20 degC 全程 |ΔT| ≤ 0.08 K；-25 degC 末端 +0.35 K。

### 4.2 产冰过程

| t / s | 0 | 5 | 10 | 15 | 20 | 25 | 30 | 35 | 36.6 |
|---|---|---|---|---|---|---|---|---|---|
| sIce (-20 degC) | 0 | 0 | 0 | 0.0002 | 0.0032 | 0.0090 | 0.0204 | 0.0354 | **0.0406** |
| sIce (-25 degC) | 0 | 0 | 0 | 0.0013 | 0.0059 | 0.0134 | 0.0321 | 0.0696 | — |
| fArea (-20 degC) | 1.0000 | 1.0000 | 1.0000 | 0.9991 | 0.9890 | 0.9687 | 0.9303 | 0.8816 | **0.8648** |
| fArea (-25 degC) | 1.0000 | 1.0000 | 1.0000 | 0.9955 | 0.9795 | 0.9538 | 0.8921 | 0.7768 | — |

- 末端冰量 **1.78e-4 kg/m² (-20 degC) / 3.66e-4 kg/m² (-25 degC)**，
  占累计产水的 **2.97% / 5.95%**。
- 冰量**单调不减**：`min d(ice)/dt = 0`，非单调点 **0 / 183**（两工况）。
- 冰的电压代价（由 `etaAct` 与 `fArea` 反算，`tmp/E_icecost.m`）：
  **+6.53 mV (-20 degC) / +13.52 mV (-25 degC)** 末端附加活化损失。
- 温度越低冰越多（4.06% → 8.37%），符合冷启动物理方向。
- 冰在 cCL 内部的空间峰值饱和度：4.06% (-20 degC) / 20.7% (-25 degC)。

### 4.3 数值健康检查（全部通过）

```
[-20 degC] finite=1 nonNeg=1 ice<=poreW=1 ice<=totW=1 satOK=1 lamOK=1 j<jlim=1 cpu=4.6 s
[-25 degC] finite=1 nonNeg=1 ice<=poreW=1 ice<=totW=1 satOK=1 lamOK=1 j<jlim=1 cpu=12.1 s
```

### 4.4 机制贡献分解（`tmp/E_compare.m`）

| 变体 | V_R20 | V_R25 | 合并 | sIce20 | sIce25 |
|---|---|---|---|---|---|
| D 版基线 | 0.02514 | 0.00958 | 0.01903 | 0.000 | 0.044 |
| E 版（关闭直接冻结，`kFreezeDirect=0`） | 0.01411 | 0.01289 | 0.01351 | 0.000 | 0.042 |
| E 版（冰通道整体关闭，`df=0,gammaIce=0`） | 0.01656 | 0.01110 | 0.01410 | 0.000 | 0.000 |
| **E 版（最终，冰通道开启）** | **0.01432** | **0.01480** | **0.01456** | **0.041** | **0.084** |

读法：
1. 仅解除参数简并 + 联合重标定（第 2 行）就把合并 RMSE 从 0.01903 降到
   0.01351，且两工况更均衡（D 版是 0.025/0.010，严重偏斜）。
2. 把"冰通道开/关"在同一组水合参数下对比（第 3 行 vs 第 4 行）：
   **-20 degC 从 0.01656 改善到 0.01432（改善 13.5%）**，
   **-25 degC 从 0.01110 变差到 0.01480（变差 33%）**。
   即冰的附加损失在 -20 degC 正是原模型缺失的那一块（D 版末端偏高 10.5 mV，
   而冰恰好给出 6.5 mV）；在 -25 degC 则叠加在本已偏低的曲线上，见 §5.1。
3. 冰是"净增损失"机制：开冰使合并 RMSE 上升 0.00046 V
   （0.01410 → 0.01456）。这部分代价是为物理自洽（低温冷启动必须产冰）
   而主动付出的，而且**开冰后合并 RMSE 仍比 D 版低 23%**。

---

## 5. 诚实的局限与遗留问题

### 5.0 [-25 degC 工况的数值鲁棒性修复 · E-2 ~ E-6]

**症状**：在 -25 degC 工况下，只要把低温结冰相关参数调强一些，程序就会中断，
且报错信息是

```
water_ice_state_ice:InvalidDiffusivity
水扩散系数出现非法数值：cell=1, T=NaN K, Dw=NaN, lambdaPEMRaw范围=[0,0]
```

这条信息**完全指错了方向**——它让人以为是水扩散系数或物性参数有问题，
实际上低温工况的真正故障是积分器发散。

**根因链**（用 `tmp/E_fail25.m` 逐层定位，栈为
`ode15s → odenonnegative → model_rhs_ice → water_ice_state_ice`）：

1. `thermal_temperature_state_ice.m` 用 **`Inf` 作为"超出模型适用范围"的哨兵**：
   - `j >= jLim`（超过极限电流密度）时 `etaCon = Inf`；
   - `kappaPEM/ACL/CCL <= 0`（离聚物完全干涸）时 `etaOhm = Inf`。
2. 于是 `Vcell = Erev - etaAct - etaOhm - etaCon = -Inf`；
3. 产热项 `qgen = j*(1.48-Vcell)/g.Ltotal = +Inf`；
4. 热方程 `dT = (...+qgen+...)/rhoCp = Inf`，并经
   `gas_transport_state_ice(...,dT,...)` 把 `dcO2` 也传染成非有限值；
5. `ode15s` 的试探步随即发散，NaN 状态在下一轮进入
   `water_ice_state_ice`，最后在第 6 节的 `Dw` 计算处报成"水扩散系数非法"。

**为什么只有 -25 degC 出问题**：-20 degC 下 cCL 的液态水 `ml` 恒为 0，
题目式(9) 的液水冻结项 `Rfreeze = kFreeze*rhoL*eps0*sL` 完全不参与；
-25 degC 下 cCL 开始出现液态水，该项才真正激活，系统刚性明显增大，
更容易把 `j` 推过 `jLim` 或让某个控制体的孔隙被填满。

**修复内容**：

| 标记 | 文件 | 内容 |
|---|---|---|
| `[E-2]` | `water_ice_state_ice.m` | 入口先判定状态是否有限。若 `T/mw/mi/lambdaCCL` 已是 NaN/Inf，直接报"水-冰模块收到非有限状态（ode15s 试探步发散）"，不再让它伪装成物性错误 |
| `[E-3]` | `pemfc_calculate_ice.m` | 在 `model_rhs_ice` 里检查打包后的 `dx`，报出**具体是哪个分量**非有限，并附带当时的 `Vcell/Erev/etaAct/etaOhm/etaCon/j/jLim` |
| `[E-4]` | `pemfc_calculate_ice.m` | `ode15s` 返回后显式判定积分是否真正完成（解被截断 / 出现 NaN / 抛出 `IntegrationTolNotMet`），报出**实际推进到的时刻**与求解器原文信息。注意：`lastwarn` 必须在调用前清空，否则会读到上一次运行的残留警告而误判 |
| `[E-5]` | `thermal_temperature_state_ice.m` | 取消 `Inf` 哨兵：① 电导率取极小下限 `kappaFloor=1e-6 S/m`；② `log` 自变量夹到 `1-1e-12`；③ 总损失上限 `voltageLossCap=2 V`，保证 `Vcell` 恒为有限值。新增标志 `voltageLossClamped` / `transportInfeasible`，由 `pemfc_calculate_ice` 汇总为 `result.health.withinModelRange` 并在超范围时发出警告 |
| `[E-6]` | `water_ice_state_ice.m` | 孔隙被液水+冰填满时不再 `error('BlockedPore')` 中断积分，而是给"用于输运/活度计算的孔隙率"取极小正下限 `1e-6*eps0`。**不改变 mv/ml/mi 的分配，不产生额外质量**；触发情况由 `out.poreFilledProjection` 记录 |

**修复后验证**：

- 交付参数下 -20/-25 degC 结果**逐位不变**：
  `V_RMSE=0.014324/0.014802 V`，`ice=1.7776e-4/3.6616e-4 kg/m²`，
  全部 `health` 标志为 1（`[E-5]` 的上限在两工况总损失只有 0.4~0.8 V 时从不生效）。
- `run_ice_model.m` 在同一会话内连续运行两次结果完全一致（验证 `[E-4]` 没有误判）。
- `tmp/E_probe25.m` 的 19 组参数探针：修复前 2 组以误导性错误崩溃，
  修复后 `interfaceFactor=3` 可正常求解（`V_RMSE=0.01606`，`sIce=12.66%`）。

**仍然存在的限制（-25 degC 参数包络）**：

| 参数 | 已验证稳定范围 | 失效表现 |
|---|---|---|
| `kFreeze` | `0 ~ 5` | `kFreeze=10` 时 `ode15s` 在 t≈36.1 s 触及最小步长（`1.14e-13`）而失败，现由 `[E-4]` 明确报为 `IntegrationFailed` |
| `directFreezeFraction` | `0.01 ~ 0.50` | 均可求解（`0.50` 时 `V_RMSE` 已升到 0.0248，属拟合变差而非崩溃） |
| `freezeNucleationLambdaFraction` | `0.10 ~ 0.50` | 均可求解 |
| `interfaceFactor` | `3 ~ 30` | 均可求解 |

`kFreeze` 过大导致最小步长失败的原因是非光滑约束：冰量触到
`atAllWaterFrozen`（冰 ≤ 孔隙水）与 `atFullPore`（冰 ≤ 孔隙容积）时
`dmi` 会硬性跳变，`ode15s` 只能不断缩步。
交付参数下 cCL 的 `mi ≈ mPore`（可动水全部冻结），因此**正好骑在该约束上**，
这也是本工况偏刚性的原因；若要支持更大的 `kFreeze`，需要把这两道硬钳位
改写为光滑容量因子，但那会改变已标定的结果，故本版未做。


### 5.1 模型对温度的敏感度高于实测（当前残差的主要来源）

两个工况的实测电压曲线几乎重合
（t=12 s：0.495 / 0.493 V；t=24 s：0.637 / 0.629 V；t=32 s：0.648 / 0.634 V），
相差仅 2 ~ 14 mV。而题目式(40) 给定 `Ea=67000 J/mol`，在 253.15 → 248.15 K
之间使 `j0` 变化 1.905 倍，对应活化损失差
`(RT/αF)·ln(1.905) ≈ 29 mV`（`Erev` 的熵项反向补偿约 4 mV，净约 25 mV）。
**在 `Ea` 与 `alpha` 均为题目给定值的前提下，这一约 25 mV 的温度敏感度
差异无法通过任何参数消除**，它是当前残差的主要来源，也解释了为什么
-20 degC 的 t=10 s 点仍有 4.69% 误差（-25 degC 同点仅 0.57%），
以及为什么冰在 -20 degC 改善拟合而在 -25 degC 恶化拟合。
若要进一步改善，需引入另一个"低温使损失减小"的机制（例如
O₂ 在离聚物/水中的亨利常数随温度下降而升高），或质疑实测曲线的一致性。
**本版没有通过调参掩盖该矛盾**：D 版把 -25 degC 拟合到 0.0096 V 的代价，
正是 -20 degC 被推到 0.0251 V；本版选择让两者均衡在 0.0143 / 0.0148 V。

### 5.2 早期电压亏损的机制归属

t=0 ~ 6 s 时 j 仅 0.006 A/cm²，实测电压已只有 0.80 V，而模型
`Erev ≈ 1.267 V`，缺口约 0.45 V。在 6 mA/cm² 下这**不可能**是欧姆损失
（需串联电阻 R ≈ 750 Ω·cm²），只能归因于活化损失——即 cCL 离聚物极度缺水
导致有效交换电流密度极低（`j0_eff(t=0) = 3.7e-3 A/m²`，对应
`etaAct(t=0) ≈ 0.42 V`）。因此更准确的表述是：

> **① 前段压降（t ≈ 6~14 s）= 离聚物缺水同时抬高 CL 质子电阻与降低有效 j0，
> 其中以活化损失为主、欧姆为辅（clResistanceScale=0.64 的 CL 电阻只占
> t=0 总损失 ~0.02 V 中的 0.017 V）；
> ② 中段回升（t ≈ 14~25 s）= 水合度上升降低 κ_CCL 并提高 j0，仍以活化为主；
> ③ t > 25 s 冰开始贡献附加活化损失（+6.5 mV @ -20 degC / +13.5 mV @ -25 degC）。**

早前"电压骤降主要由欧姆内阻造成"的说法需要按此修正。

### 5.3 其他局限

- **冰只生成在 cCL。** 直接冻结通道只作用于 cCL（产水位）。cGDL/aCL 的冰
  仍只能通过题目式(9) 的液水路径产生，而在本算例中这些层始终未达液水饱和。
  真实冷启动中 cGDL 也会因水蒸气凝华积冰，这需要 `kmCond > 0` 的有限速率
  凝华模型（`[DSH-1]` 已实现但默认关闭）。

- **`clResistanceScale` 与 `directFreezeFraction` 是唯象标定量**，
  论文中应明确声明其标定依据，不应表述为独立辨识的物性值。
  本版共 4 个水合/欧姆标定量 + 2 个产冰标定量，全部由两工况全时段曲线联合标定。

---

## 6. 相对 D 版的完整改动清单

| 文件 | 改动 |
|---|---|
| `water_ice_state_ice.m` | `[PHASE-1]` `interfaceFactor` 参数化；`[PHASE-4b]` 去掉可动水门槛；新增 `kFreezeDirect`；新增 `freezeDiagnostics` 诊断输出；`[E-2]` 入口非有限状态判定；`[E-6]` 孔隙填满时改为投影而非中断 |
| `pemfc_calculate_ice.m` | 默认值更新（`j0Ref=0.30`、`mHyd=1.55`、`lamHydRef=1.0`、`clResistanceScale=0.64`、`gammaIce=3.5`、`directFreezeFraction=0.15`、`freezeNucleationLambdaFraction=0.25`）；新增 `interfaceFactor`、`kFreezeDirect` 选项传递；`[E-3]` 状态导数有限性判定；`[E-4]` 积分完成情况判定与适用范围警告 |
| `thermal_temperature_state_ice.m` | `[D-2]`→`[E-1]` 更正 `lamHydRef` 的错误注释，说明其与 `j0Ref` 的简并；`[E-5]` 取消 `Inf` 电压损失哨兵，改为有限化处理并输出超范围标志 |
| `run_ice_model.m` | 全部选项显式传入，汇总表增加 RMSE、压降/回升达成率 |
| `run_ice_model_voltage_losses.m` | 同一组参数；新增冰饱和度、冰面积因子、冰面密度与电压对照图 |
| `README_E.md` | 本文件 |

## 7. 使用的诊断脚本（`../tmp/`，均非交付代码）

| 脚本 | 用途 |
|---|---|
| `E_gates.m` | 逐时打印冻结通道三重门槛，定位冰为 0 的原因 |
| `E_water.m` | 分层水账（产水/膜/离聚物/冰/孔隙） |
| `E_ccl0.m` | cCL 初始水状态与 `fHyd` 参考点核对 |
| `E_exp.m` | 打印两工况实测 j / V / T 剖面 |
| `E_iface.m` | `interfaceFactor` 扫描 |
| `E_scan2.m` | 直接冻结 × 气相边界扫描 |
| `E_joint.m` / `E_opt1.m` / `E_opt2.m` / `E_final.m` | 逐级联合标定 |
| `E_final_check.m` | 最终配置完整验证 |
| `E_compare.m` | 机制贡献分解（D 基线 vs E 各开关） |
