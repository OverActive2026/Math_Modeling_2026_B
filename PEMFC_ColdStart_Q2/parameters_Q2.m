function p = parameters_Q2()
%PARAMETERS_Q2 问题2参数入口；未确认物理参数保留 NaN。
p.Ncell = 5;
p.A_cm2 = 25; % cm^2
p.A = 25e-4; % m^2
p.h = 40; % W/(m^2 K)
p.j_max = 0.5; % A/cm^2
p.q_max = 20; % C/cm^2
p.V_min_limit = 0.30; % V
p.ice_limit = 0.99; % -
p.T_success = 273.15; % K
p.T0 = 263.15; % K
p.Tamb = p.T0; % K
p.dt = 0.01; % s
p.t_max = 300; % s
p.k_eff_stack = NaN; % TODO: W/(m K)
p.delta_stack = NaN; % TODO: m
p.V_cell = NaN; % TODO: m^3
% TEST PLACEHOLDER ONLY: 仅验证软件框架，非物理参数。
p.test.enabled = false; % 正式参数默认禁止占位计算
p.test.k_eff_stack = 1; % W/(m K), 测试值
p.test.delta_stack = 1; % m, 测试值
p.test.V_placeholder = p.V_min_limit; % V, 测试哨兵
end
