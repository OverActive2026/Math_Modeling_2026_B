function state = initial_state(p)
%INITIAL_STATE 五片单电池各自独立的初始状态。
n = p.Ncell;
state.T = ones(n,1)*p.T0; % 平均温度 [K]
state.V = NaN(n,1); % TODO: 问题1初始电压 [V]
state.ice_max = zeros(n,1); % 最大局部冰体积分数 [-]
state.lambda_mem = NaN(n,1); % TODO: 问题1膜含水量 [-]
state.m_v = NaN(n,1); % TODO: 气态水初值 [kg/m^3]
state.m_l = NaN(n,1); % TODO: 液态水初值 [kg/m^3]
state.m_i = NaN(n,1); % TODO: 冰质量初值 [kg/m^3]
end
