%--------------------------------------------------------------%
% @function: 计算质子交换膜中吸附态水的扩散系数，即题目中式(24)
% @author:   PJ, GPT
% @date:     20260923
% @input:    Lambda->膜含水量
%            T->温度
%            p->参数列表
% @output:   Dw->扩散系数
%--------------------------------------------------------------%

function Dw = membrane_water_diffusivity(lambda, T, p)
    
    % 温度单位为开尔文，需大于零
    if any(T(:) <= 0)
        error('Temperature must be greater than zero.');
    end
    
    polyLambda = p.water.Dmem_c(1) + p.water.Dmem_c(2).*lambda + p.water.Dmem_c(3).*lambda.^2 + p.water.Dmem_c(4).*lambda.^3;
    
    Dw = p.water.Dmem_prefactor.*exp(p.water.Dmem_Eterm.*(1/303.15 - 1./T) ).*polyLambda;
    
    % 扩散系数需大于零
    if any(Dw(:) < 0)
        warning('membrane_water_diffusivity:NegativeValue', 'Membrane-water diffusivity became negative. Check lambda range.');
    end

end
