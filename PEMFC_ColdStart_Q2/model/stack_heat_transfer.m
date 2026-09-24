function Qstack = stack_heat_transfer(T,p)
%STACK_HEAT_TRANSFER 五片电池之间的净导热及两端对流散热 [W]。
% 正号表示该片从邻片或环境得到热量。
T = T(:);
if p.Ncell ~= 5 || numel(T) ~= 5 || any(~isfinite(T))
    error('Q2:InvalidTemperature','T 必须是五片的有限温度列向量 [K]。');
end
if ~isfinite(p.k_eff_stack) || p.k_eff_stack <= 0 || ...
        ~isfinite(p.delta_stack) || p.delta_stack <= 0
    error('Q2:MissingStackThermalParameters', ...
        'k_eff_stack 或 delta_stack 未确定；正式计算不能使用 NaN。');
end
G = p.k_eff_stack*p.A/p.delta_stack; % 相邻单片导热导纳 [W/K]
Qstack = zeros(5,1);
Qstack(1) = G*(T(2)-T(1))-p.h*p.A*(T(1)-p.Tamb);
Qstack(2) = G*(T(1)-T(2))+G*(T(3)-T(2));
Qstack(3) = G*(T(2)-T(3))+G*(T(4)-T(3));
Qstack(4) = G*(T(3)-T(4))+G*(T(5)-T(4));
Qstack(5) = G*(T(4)-T(5))-p.h*p.A*(T(5)-p.Tamb);
end
