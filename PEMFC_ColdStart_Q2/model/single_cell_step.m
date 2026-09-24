function newState = single_cell_step(oldState,j,q_stack,p)
%SINGLE_CELL_STEP 问题1单片物理模型的待接入接口。
% j [A/cm^2]；q_stack 目前为电堆传来的净热功率 [W]。
% TODO: 接入真实能量方程时，核实 V_cell 后用 q_stack/V_cell [W/m^3]。
if ~isscalar(j) || ~isfinite(j) || j < 0 || j > p.j_max
    error('Q2:InvalidCurrent','电流密度超出允许范围 [A/cm^2]。');
end
if ~isscalar(q_stack) || ~isfinite(q_stack)
    error('Q2:InvalidHeatSource','电堆净热功率必须为有限标量 [W]。');
end
if ~isfield(p,'test') || ~p.test.enabled
    error('Q2:MissingSingleCellPhysics', ...
        '问题1真实物理方程尚未接入，正式计算已停止。');
end

% TEST PLACEHOLDER ONLY：仅测试框架；不预测温度、结冰或电压。
% 其余状态保持原样；电压取已知下限作为可识别的测试哨兵值。
newState = oldState;
newState.V = p.test.V_placeholder;

% TODO 1: 气体守恒
% TODO 2: 氢气/氧气浓度计算
% TODO 3: 水蒸气守恒
% TODO 4: 液态水守恒
% TODO 5: 冰相守恒
% TODO 6: 冻结和融化速率
% TODO 7: PEM膜水冻结
% TODO 8: 冰占孔修正
% TODO 9: 有效扩散系数修正
% TODO 10: 催化层有效反应面积/ECA修正
% TODO 11: 活化极化
% TODO 12: 欧姆极化
% TODO 13: 浓差极化
% TODO 14: 单片电压
% TODO 15: 电化学产热
% TODO 16: 相变潜热
% TODO 17: 能量守恒方程
end
