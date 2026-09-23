%-------------------------------------------------------------------%
% @function: 使用Springer经验公式计算质子交换膜电导率，即题目中式(42)
% @author:   PJ, GPT
% @date:     20260923
% @input:    Lambda->膜含水量
%            T->温度
%            p->参数列表
% @output:   kappa->膜电导率
%-------------------------------------------------------------------%

function kappa = membrane_conductivity(lambda, T, p)
    
    % 温度单位为开尔文，需大于零
    if any(T(:) <= 0)
        error('Temperature must be greater than zero.');
    end

    kappa = (p.water.kappa_a0.*lambda + p.water.kappa_a1).*exp(p.water.kappa_E.*(1/303.15 - 1./T));

    % 膜电导率应大于零，非正值通常说明Lambda超出经验公式的适用范围
    if any(kappa(:) <= 0)
        warning('membrane_conductivity:NonPositiveConductivity', 'Non-positive PEM conductivity detected. Check membrane water content lambda.');
    end

end
