%---------------------------------------------------------------------------------%
% @function: 根据温度、压力和气相孔隙率计算多孔介质中的有效气体扩散系数，即题目中式(16)
% @author:   PJ, GPT
% @date:     20260923
% @input:    Dref->参考状态下的气体扩散系数
%            T->局部温度
%            pGas->局部气体绝对压力
%            eps_g->局部气相孔隙率
%            p->参数列表
% @output:   Deff->多孔介质中的有效气体扩散系数
%---------------------------------------------------------------------------------%

function Deff = calc_gas_diffusivity(Dref, T, pGas, eps_g, p)

    % 温度单位为开尔文，需大于零
    if any(T(:) <= 0)
        error('Temperature must be greater than zero.');
    end
    
    if any(pGas(:) <= 0)
        error('Gas pressure must be greater than zero.');
    end
    
    if any(eps_g(:) < 0)
        error('Negative gas porosity detected.');
    end

    % 使用温度修正、压力修正和Bruggeman孔隙率修正
    Deff = Dref .* (T ./ p.const.Tref).^p.trans.Texp .* (p.const.p0 ./ pGas) .* eps_g.^p.trans.epsExp;

end
