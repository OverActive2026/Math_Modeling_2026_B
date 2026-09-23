%--------------------------------------------------------------%
% @function: 在已知冰含量下计算多孔介质中水蒸气与液态水的平衡分配
% @author:   PJ, GPT
% @date:     20260923
% @input:    mw->总水质量浓度
%            mi->冰质量浓度
%            T->温度
%            eps0->干燥状态下的孔隙率
%            p->参数列表
% @output:   mv->水蒸气质量浓度
%            ml->液态水质量浓度
%            eps_g->剩余气相孔隙率
%            info->水相平衡诊断信息
%--------------------------------------------------------------%

function [mv, ml, eps_g, info] = water_phase_equilibrium(mw, mi, T, eps0, p)

    % 总水由水蒸气、液态水和冰组成，各部分质量浓度不能为负
    if any(mw(:) < 0) || any(mi(:) < 0)
        error('Water concentrations cannot be negative.');
    end

    if any(mi(:) > mw(:))
        error('Ice concentration cannot exceed total water.');
    end

    rhoL = p.ice.rho_liquid;
    rhoI = p.ice.rho_ice;

    Mw = p.const.Mw;
    R  = p.const.R;

    % 扣除冰后的可流动水由水蒸气和液态水组成
    mMobile = mw - mi;

    % 计算冰占据体积后的剩余孔隙率
    epsAfterIce = eps0 - mi./rhoI;

    if any(epsAfterIce(:) < 0)
        error('Ice volume exceeds available pore volume.');
    end

    % 计算当前温度下的饱和蒸气质量浓度系数
    psat = water_psat(T);
    A = Mw .* psat ./ (R .* T);

    % 假设没有液态水时，剩余孔隙能容纳的饱和水蒸气
    mvCapacity = A .* epsAfterIce;

    % 预分配水蒸气和液态水向量
    mv = zeros(size(mw));
    ml = zeros(size(mw));
    dryOrUnsaturated = mMobile <= mvCapacity;

    % 未饱和网格中的可流动水全部作为水蒸气
    mv(dryOrUnsaturated) = mMobile(dryOrUnsaturated);
    ml(dryOrUnsaturated) = 0;

    % 饱和网格中同时存在水蒸气和液态水
    idx = ~dryOrUnsaturated;
    den = 1 - A(idx)./rhoL;

    if any(den <= 0)
        error('Invalid vapor-liquid equilibrium denominator.');
    end

    ml(idx) = (mMobile(idx) - A(idx).*epsAfterIce(idx)) ./ den;
    mv(idx) = mMobile(idx) - ml(idx);

    % 根据液态水和冰更新各相体积分数
    [eps_g,eps_l,eps_i,poreInfo] = calc_porosity(eps0,ml,mi,p);

    % 输出中间量，用于检查相平衡和孔隙状态
    info.psat = psat;
    info.A = A;

    info.mMobile = mMobile;
    info.mvCapacity = mvCapacity;

    info.eps_l = eps_l;
    info.eps_i = eps_i;

    info.isSaturated = idx;
    info.isFlooded = poreInfo.isFlooded;

end
