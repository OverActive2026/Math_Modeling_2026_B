%--------------------------------------------------------------%
% @function: 使用Buck公式计算水的饱和蒸气压，即题目中式(28)
% @author:   PJ, GPT
% @date:     20260923
% @input:    T->温度
% @output:   psat->饱和蒸气压
%--------------------------------------------------------------%

function psat = water_psat(T)

    % 温度单位为开尔文，需大于零
    if any(T(:) <= 0)
        error('water_psat:InvalidTemperature', 'Temperature must be given in Kelvin and greater than zero.');
    end
    
    % 判断水的温度
    Tc = T - 273.15;
    psat = zeros(size(T));
    idxWarm = Tc >= 0;
    idxCold = ~idxWarm;
    
    % 温度不低于冰点时使用液态水分支
    Tw = Tc(idxWarm);
    psat(idxWarm) = 611.21.*exp((18.678 - Tw./234.5).*Tw./(257.14 + Tw));

    % 温度低于冰点时使用冰面饱和蒸气压分支
    Tcold = Tc(idxCold);
    psat(idxCold) = 611.15.*exp((23.036 - Tcold./333.7).*Tcold./(279.82 + Tcold));

end
