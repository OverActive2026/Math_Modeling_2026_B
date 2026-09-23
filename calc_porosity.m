%------------------------------------------------------------------------%
% @function: 由液态水和冰的质量浓度计算多孔介质中各相体积分数
% @author:   PJ, GPT
% @date:     20260923
% @input:    eps0->干燥状态下的孔隙率
%            ml->液态水质量浓度
%            mi->冰质量浓度
%            p->参数列表
% @output:   eps_g->剩余气相孔隙率
%            eps_l->液态水体积分数
%            eps_i->冰体积分数
%            info->孔隙状态诊断信息
%------------------------------------------------------------------------%

function [eps_g, eps_l, eps_i, info] = calc_porosity(eps0, ml, mi, p)

    % 液态水和冰的质量浓度不能为负
    if any(ml(:) < 0)
        error('calc_porosity:NegativeLiquidWater', 'Liquid-water concentration cannot be negative.');
    end

    if any(mi(:) < 0)
        error('calc_porosity:NegativeIce', 'Ice concentration cannot be negative.');
    end

    rhoL = p.ice.rho_liquid;
    rhoI = p.ice.rho_ice;

    if ~isfinite(rhoL) || rhoL <= 0
        error('Invalid liquid-water density.');
    end

    if ~isfinite(rhoI) || rhoI <= 0
        error('Invalid ice density.');
    end

    eps_l = ml ./ rhoL;
    eps_i = mi ./ rhoI;

    % 气相孔隙率等于干孔隙率减去液态水和冰占据的体积
    eps_g = eps0 - eps_l - eps_i;

    % 不对负孔隙率进行截断，保留原始结果以便识别泛水或非法状态
    info.isFlooded = eps_g <= 0;
    info.isValid   = eps_g >= 0;

    info.minGasPorosity = min(eps_g(:));
end
