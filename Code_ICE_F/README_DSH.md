# Code_ICE_DSH —— Code_ICE 的机理补全版（含验证结论）

本目录由 `Code_ICE` 复制而来，只做**机理层面**的修改，不改动标定思路。
所有新增改动均带 `[DSH-*]` 标记，新增开关取值均可使模型**严格退化回 Code_ICE**。

## 一、诊断结论（本次验证得出，是全部改动的前提）

用 MATLAB 实跑 `Code_ICE`（`j0Ref=0.0049645, Ea=7028.98, tauHyd=10,
kFreeze=1, kMelt=0, gammaIce=0.1`）并做水收支审计，得到三条硬结论：

1. **多孔层几乎没有水。** 36.6 s 累计产水 5.99e-3 kg/m²，
   而多孔层孔隙水只有 7.1e-8 kg/m²（约万分之一）；
   这些水全部进了 PEM（4.59e-3 kg/m²）。
2. **冰在任何时刻都为零**（−20 ℃ 全程 `s_ice = 0`；−25 ℃ 末值 0.05）。
   原因是 **cCL 水活度全程只升到 0.64~0.79，从未达到 1**，
   因此按式(25)(26) 永远不生成液态水，冰也就没有原料。
3. **后段电压回升不可能来自温度。** 用实验电压反解"有效 j0·fArea"，
   t=20→36.6 s 需要它上升 ×6.8，而温度项只能给 ×1.1；
   要全靠温度解释需要 Ea ≈ 1.94e5 J/mol（题目给定值的 3 倍）。
   整个电压里只剩 `η_act` 一项能提供这个变化，而
   `η_act` 里只有 `fArea=(1-s)^γ` 是随时间可变的 → **必须靠冰的生成/消融**。

**根本原因**：参考模型式(16)(19) 给的水蒸气输运能力过强
（实测边界净流失 ≈ 产水量的 1.3~1.6 倍），cCL 无法过饱和。
因此原 `Code_ICE` 的"先跌后升"是**由 cCL 离聚物水合（`HYD-1`）人造出来的**：
`water_ice_state_ice.m` 原来在 `j>0` 时把水合目标直接设为
`lambdaSaturationCCL=4.837`，与局部水活度（仅 0.007~0.03）无关，
等于凭空向离聚物注入结合水，从而造出 η_ohm,CCL 由 0.135 V 降到 0.051 V
的那段"回升"。

## 二、本目录的修改

| 标记 | 位置 | 内容 |
|---|---|---|
| `[DSH-1]` | `water_ice_state_ice.m` | 有限速率冷凝/蒸发 `Rcond=kmCond*(mv−A·eps_g)`；`kmCond=0` 退化为瞬时平衡 |
| `[DSH-1a]` | `water_ice_state_ice.m` | 多孔层水蒸气有效扩散系数折减因子 `fRet`；`fRet=1` 退化 |
| `[DSH-2]` | `water_ice_state_ice.m` | 冻结速率加入过冷度因子 `(Tf−T)/dTref`，升温后冻结自然减弱 |
| `[DSH-3]` | `water_ice_state_ice.m` | 冰的升华/凝华 `Rsub`（T<0 ℃ 时冰量唯一可下降的通道），使冰堵可逆 |
| `[DSH-4]` | `water_ice_state_ice.m` | 水合目标改由局部水活度决定，删除 `j>0 → λ=λ_sat` 开关 |
| `[DSH-5]` | `thermal_temperature_state_ice.m` | 潜热按过程拆开：`Lf(Rfreeze−Rmelt)`、`−Ls·Rsub`、`Lv·Rcond` |
| `[DSH-6]` | `water_ice_state_ice.m` | **关键修复**：水合供给率由"反应产水速率"限制，并取消"孔隙水存量"钳位（详见下节） |
| `[DSH-DIAG]` | `water_ice_state_ice.m` | 输出边界水通量 `NwLeftBoundary/NwRightBoundary/NliqBoundary*`、`waterProductionRate` |

### `[DSH-6]` 是最关键的一处（原模型电压崩为负的根因）

原代码用**孔隙水存量**限制离聚物水合：

```matlab
lambdaAvailable = min((mw(g.idx_cCL)-miState(...))/cCLIonomerWaterCoefficient);
lambdaUpper     = min(lambdaScale,max(lambdaAvailable,0));      % = mw-mi 的存量
lambdaCCLTarget = min([lambdaCCLEquilibrium,lambdaUpper,lambdaSaturationCCL]);
existingPoreSupplyRate = max(lambdaUpper-lambdaCCLStateRaw,0)/tauHyd;
```

问题在于 cCL 的可迁移水**存量**长期只有约 1e-4 kg/m²（反应水很快被水蒸气带出电池），
于是 `lambdaUpper ≈ 0.086`，把水合目标压到 0.086 以下，
`dlambdaCCL` 被限制在 ~2e-5 1/s —— **36.6 s 只能让 λ 上升 1e-4**。
后果链条：

```
λ_CCL 锁死在初值 1.0
  → κ_CCL = 0.3^1.5 × (0.5139×1 − 0.326) × exp(...) ≈ 0.018 S/m（干态）
  → η_ohm,CCL 达 0.7 V 以上
  → Vcell < 0（电压崩为负）
```

而质量上完全不需要这样限制：把 λ 从 1 提到 λ_sat=4.837 只需要
`11.61 × 3.837 × 11.3e-6 = 5.0e-4 kg/m²`，仅占 36.6 s 累计产水（5.99e-3 kg/m²）的 **8%**。
正确做法是**用工况的产水速率**限制供给：

```matlab
productionSupplyRate = reactionWaterRateCCL/cCLIonomerWaterCoefficient;  % ≈0.711 1/s
poreSupplyRate = max(mPore(g.idx_cCL(1))-mv(g.idx_cCL(1)),0)/ ...
    (cCLIonomerWaterCoefficient*tauHyd);
dlambdaCCL = min(hydrationKineticRate, poreSupplyRate+productionSupplyRate);
lambdaCCLTarget = min(lambdaCCLEquilibrium,lambdaSaturationCCL);   % 不再用存量钳位
```

同时 `[DSH-6]` 把 cCL 的反应水源项按优先级拆分：
反应水先满足离聚物水合，剩余部分才进入孔隙相（可排水/可冻结），
即 `Sw(cCL) = Sw(cCL) − dmIonomer`，而总水方程 `mw = mIon + mv + ml + mi` 保持守恒。

**效果（−20 ℃，fRet=1，MATLAB 实跑）**：

| 量 | 修复前 | 修复后 |
|---|---|---|
| λ_CCL 末值 | 1.0001 | **8.087** |
| κ_CCL 末值 [S/m] | 0.0181 | **0.3251** |
| Vcell 末值 [V] | −0.3503 | **0.5385** |
| V_RMSE [V] | 0.7280 | **0.0754** |
| T_RMSE [K] | 3.308 | **0.312** |

−25 ℃：V_RMSE 0.0676、T_RMSE 0.388、λ_CCL 7.123、V 末值 0.5371。
此时 λ_CCL 与原 `Code_ICE` 的 8.225/7.162 基本一致，
**说明原模型的离聚物水合结果本身是对的，坏掉的只是它的"供给约束"逻辑**；
也说明我先前对 `HYD-1` 的"凭空造水"判断需要修正——真正的问题是
`lambdaUpper` 这个存量钳位，而不是水合目标采用饱和值。


新增求解选项（`pemfc_calculate_ice.m`）：
```matlab
simOpt.kmCond = 0;   % 冷凝传质系数 [1/s]，0=瞬时平衡
simOpt.dTref  = 20;  % 冻结过冷度参考值 [K]
simOpt.nRet   = 0;   % 液水对气相孔隙的阻滞指数
simOpt.rRet   = 0;   % 液水滞留对水蒸气排出的阻滞强度
simOpt.fRet   = 1;   % 多孔层水蒸气扩散直接折减因子，1=原模型
```

## 三、验证结果（−20 ℃，MATLAB 实跑）

| `fRet` | V_RMSE | max s_ice(cCL) | max 水活度 | λ_CCL 末值 | V 末值 |
|---|---|---|---|---|---|
| 1.000（原模型） | 0.7280 | 0.0000 | 0.653 | 1.000 | −0.350 |
| 0.300 | 0.7251 | 0.1375 | 0.939 | 1.000 | −0.345 |
| 0.200 | 0.7261 | 0.2116 | 0.958 | 1.000 | −0.345 |
| 0.100 | 0.7278 | 0.2895 | 0.974 | 1.000 | −0.346 |
| 0.050 | 0.7301 | 0.3190 | 0.976 | 1.000 | −0.348 |

**结论（重要）**：`fRet` 确实能把水留进 cCL 并生成冰（s_ice 达到 0.14~0.32），
但**电压反而崩到负值**。原因是水账里出现了硬冲突：

```
在 j=0.1 A/cm²、-20 ℃ 下
  累计产水            = 9.1e-4 kg/m²
  cCL 孔隙体积(面密度) = 4.76e-6 m³/m² → 最多容纳冰 4.4e-3 kg/m²
  若 s_ice ≈ 0.3 则冰 = 1.3e-3 kg/m²   <-- 已经超过累计产水 9.1e-4
  同时离聚物 λ=4.84 需要 0.28 kg/m³ × 11.3e-6 = 3.2e-3 kg/m²
```

**产水总量不足以同时满足"结冰"和"离聚物水合"。** 冰一旦生成就把可迁移水吃光，
`mPore` 降到 1 kg/m³ 量级，`lambdaUpper≈0.086`，λ_CCL 被锁在 1.0，
κ_CCL 停在 0.018 S/m（干态），η_ohm,CCL 变成 0.7 V 以上 → 电压为负。

**下一步（尚未完成）**，三选一或组合：

1. **降低冰对水的竞争**：给出 cCL 的"优先级"——反应水先满足离聚物水合到
   `lambdaSaturationCCL`，剩余才进入孔隙相（冻结/排水）。这符合离聚物
   亲水、孔隙疏水的物理图像，也是本问题最合理的做法。
2. **提高离聚物容量**：附件1 的低温最大非冻结含水量关系给的是 ~4.84，
   但若按"离聚物体积分数 0.3 + 溶胀"取值会更大；需要回到附件1 核对
   `lambdaSaturationCCL` 的适用对象（膜 vs CL 离聚物）。
3. **把 `lambdaCCL` 从状态量改为代数量**：λ 由局部水活度瞬时决定
   （`lambdaCCLEquilibrium`），只把"何时达到该活度"交给水输运方程，
   避免 λ 与 mw 互相抢水造成刚性。

在此之前，`fRet` 不要启用（保持 1.0 = 原模型），
其余 `[DSH-2][DSH-3][DSH-5]` 是无条件正确的物理补全，可以保留。

## 四、诊断脚本

位于上一级 `tmp/` 目录：

- `audit_water.m`：水收支审计（产物水去向）
- `check_pore.m`：逐层孔隙水/水活度剖面
- `run_diag2.m`：电压四项分解 + 模型/实验 Tafel 斜率对比
- `run_diag3.m`：由实验电压反解"有效 j0·fArea"
- `dsh_run.m`：本目录模型的标准运行 + 诊断表
- `dsh_sweep4.m`：`fRet` 扫描
- `dsh_lam.m`：λ_CCL 状态接线核查
