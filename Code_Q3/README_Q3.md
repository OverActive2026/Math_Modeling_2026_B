# 问题3：外部电热丝辅助冷启动模型

本目录由 `Code_Q2_PhaseFix` 完整复制而来，保留问题2相变修正、五片热耦合、
独立端板节点和全部约束。问题3新增文件不会修改原 `Code_Q2_PhaseFix`。

## 1. 已加入的模型

题目给出每片电热丝恒定面功率密度

```text
0 <= q_k <= 1 W/cm^2,  k=1,...,5
```

现有一维状态只显式离散MEA五层，双极板热容已均匀折算进每片的等效体积热容，
没有独立双极板温度节点。为保持这一降阶口径，电热丝面功率先换算为
`q''_k=10^4 q_k [W/m^2]`，再定义均匀等效体积热源

```text
q_aux,k(x,t) = q''_k(t)/L_MEA
```

因此严格满足

```text
integral_0^L q_aux,k dx = q''_k
```

源项直接进入原能量方程：

```text
rhoCp_eff*dT/dt = div(k_eff*grad(T)) + q_gen + q_phase + q_aux
```

这里的“均匀”不是说真实电热丝位于MEA内，而是把双极板内热源投影到当前已经
均匀化了双极板热容的单片温度状态。若后续需要研究电热丝与MEA之间的瞬态温差，
应再增加显式双极板/电热丝温度节点。

## 2. 两种题设模式

- `preheat`：纯预加热，`j(t)=0`；五片平均温度同时达到0 ℃即完成预热，随后可加载电流。
  若电热丝先于该时刻关闭，即使之后因片间热再分配达到0 ℃，也判为不满足纯预热定义。
- `coheat`：从 `t=0` 同步加载电流和辅助热；电流固定为
  `j(t)=min(0.005t,0.3) A/cm^2`。

两种模式的电热丝都在 `0<=t<t_h` 内保持给定恒功率，之后关闭。模型仍检查问题2
的20 C/cm²电荷、0.5 A/cm²电流、0.30 V最低电压、冰体积分数和孔隙占用约束。

## 3. 能耗

每片能耗和总能耗按题目式(11)解析计算：

```text
E_k = A*q_k*min(t_h,t_stop)
E_aux = sum(E_k)
```

代码使用 `A=25 cm^2`，所以输入 `q_k [W/cm^2]` 后直接得到瓦特和焦耳。解析积分
不依赖ODE输出步长，也不会因加热关断端点的离散采样产生误差。

## 4. 主要文件

- `pemfc_stack5_simulate_q3.m`：问题3五片电堆核心，接入热源并输出能耗。
- `q3_auxiliary_heating.m`：五片恒功率及定时关断模型。
- `run_q3_model.m`：纯预加热与恒功率协同启动的入口脚本；其中参数只是候选值，尚未优化。
- `test_q3_auxiliary_heating.m`：功率开关、面/体热源守恒、能耗和零热源退化检查。
- `pemfc_stack5_simulate.m`：保留的问题2核心，用于零热源对照。

## 5. 调用示例

```matlab
opt = struct('initialTemperatureC',-30, ...
             'ambientTemperatureC',-30, ...
             'plot',true,'verbose',true);

[result,model] = pemfc_stack5_simulate_q3( ...
    'coheat',[0.55 0.35 0.30 0.35 0.55],80,opt);
```

关键输出包括：

- `result.totalAuxiliaryEnergyJ`
- `result.auxiliaryEnergyJCell`
- `result.actualHeatingTimeS`
- `result.stopTimeS`
- `result.minimumCellVoltageV`
- `result.maximumIceVolumeFraction`
- `result.maximumCellTemperatureSpreadC`
- `result.success` 与 `result.stopReason`

## 6. 下一阶段

当前完成的是问题3(1)的热源接入和问题3(2)的可评估接口。下一步可把
`[q1,...,q5,t_h]` 接到有界优化器：先利用五片结构对称性约束
`q1=q5, q2=q4` 降成四个变量，再分别优化纯预热和协同启动；最终候选解需用
完整网格 `[6 13 19 19 6]` 与正式容差复算后填写表4。

## 7. 当前数值连通性检查

在粗网格 `[2 3 4 4 2]`、初温和环境均为-30 ℃的代表工况中：

- 纯预热 `q=[1,1,1,1,1] W/cm^2` 于约33.50 s达到目标，辅助能耗约4187.09 J；
- 协同加热 `q=[0.6,0.5,0.45,0.5,0.6] W/cm^2` 于约62.57 s达到目标，
  辅助能耗约4145.58 J，历史最低单片电压约0.659 V，最大冰体积分数约0.0714。

这些数字只证明热源、事件和能耗链路可运行，不是优化结果，也不能直接填写表4。
