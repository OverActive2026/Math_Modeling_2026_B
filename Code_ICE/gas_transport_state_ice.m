%-------------------------------------------------------------------------------
% @function: 计算含冰孔隙中的氢气/氧气扩散、反应消耗及状态导数
% @author:   PJ, GPT
% @date:     20260923
% @input:    T,cH2,cO2,mw->当前状态
%            dmw,dT->总水和温度状态导数
%            water,thermal->水冰模块和热模块的中间结果
%            j,g,s->电流密度、网格和状态编号
% @output:   dcH2,dcO2->气体浓度状态导数 [mol/(m^3 s)]
%            out->扩散通量、源项、动态孔隙率和守恒检查
%
% [ICE-3] 冰与液态水共同占据孔隙，并在气体存储项中计入deps_g/dt。
%-------------------------------------------------------------------------------

function [dcH2,dcO2,out] = gas_transport_state_ice( ...
    T,cH2,cO2,mw,dmw,dT,water,thermal,j,g,~)

    %% 1. 输入处理

    T = T(:);
    cH2 = max(cH2(:),0);
    cO2 = max(cO2(:),0);
    mw = max(mw(:),0);
    dmw = dmw(:);
    dT = dT(:);
    R = 8.314;    % [J/(mol K)]
    F = 96485;    % [C/mol]
    Mw = 0.018;   % [kg/mol]
    rhoLiquid = 990;
    rhoIce = 920;


    %% 2. 有效气体扩散系数
    % [ICE-3] water.eps_g已经同时扣除了液水体积和冰体积。

    idxH2 = g.idx_H2;
    idxO2 = g.idx_O2;
    epsH2 = water.eps_g(idxH2);
    epsO2 = water.eps_g(idxO2);
    DH2 = 1.10e-4.*(T(idxH2)/298.15).^1.75.*epsH2.^1.5;
    DO2 = thermal.DO2;


    %% 3. 气体有限体积扩散

    cathodeMoles = 0.233/31.998e-3+0.767/28.014e-3;
    yO2 = (0.233/31.998e-3)/cathodeMoles;
    cH2Boundary = 101325/(R*T(g.idx_aGDL(1)));
    cO2Boundary = yO2*101325/(R*T(g.idx_cGDL(end)));

    [diffH2,JH2,diffInfoH2] = regional_diffusion( ...
        cH2,DH2,g,idxH2,'dirichlet',cH2Boundary,'noflux',0);
    [diffO2,JO2,diffInfoO2] = regional_diffusion( ...
        cO2,DO2,g,idxO2,'noflux',0,'dirichlet',cO2Boundary);


    %% 4. 电化学反应源项
    % 在给定电流模式下，物质消耗由Faraday定律唯一决定。
    % 冰面积修正已进入活化损失，不再对这些源项重复修正。

    SH2 = zeros(g.N,1);
    SO2 = zeros(g.N,1);
    SH2(g.idx_aCL) = -j/(2*F*g.layerThickness(2));
    SO2(g.idx_cCL) = -j/(4*F*g.layerThickness(4));
    gasRhsH2 = diffH2+SH2(idxH2);
    gasRhsO2 = diffO2+SO2(idxO2);


    %% 5. 动态气相孔隙率
    % 气体存储项：d(eps_g*c)/dt=-dJ/dx+S。
    %
    % [ICE-3] 未饱和时eps_g=eps0-mi/rhoIce，故：
    %   deps_g/dt=-dmi/(rhoIce*dt)
    % 饱和时联立气液平衡可得：
    %   eps_g=[rhoL*eps0-mw+mi*(1-rhoL/rhoI)]/(rhoL-A)
    % 再分别对mw、mi和A(T)求偏导。

    depsGdt = -water.dmiFull/rhoIce;
    saturated = water.ml > 0;

    if any(saturated)
        eps0Layer = [0.8,0.3916,0,0.4207,0.8];
        eps0 = eps0Layer(g.layerId(:)).';
        eps0 = eps0(:);
        deltaT = 1e-3;
        psat = buck_psat(T);
        psatPlus = buck_psat(T+deltaT);
        psatMinus = buck_psat(T-deltaT);
        A = Mw.*psat./(R.*T);
        Aplus = Mw.*psatPlus./(R.*(T+deltaT));
        Aminus = Mw.*psatMinus./(R.*(T-deltaT));
        dAdT = (Aplus-Aminus)/(2*deltaT);

        denominator = rhoLiquid-A;
        numerator = rhoLiquid.*eps0-mw + ...
            water.miFull.*(1-rhoLiquid/rhoIce);
        depsDmw = -1./denominator;
        depsDmi = (1-rhoLiquid/rhoIce)./denominator;
        depsDA = numerator./denominator.^2;
        depsGdt(saturated) = ...
            depsDmw(saturated).*dmw(saturated) + ...
            depsDmi(saturated).*water.dmiFull(saturated) + ...
            depsDA(saturated).*dAdT(saturated).*dT(saturated);
    end

    dcH2 = (gasRhsH2-cH2.*depsGdt(idxH2))./epsH2;
    dcO2 = (gasRhsO2-cO2.*depsGdt(idxO2))./epsO2;


    %% 6. 反应守恒检查

    source.SH2 = SH2;
    source.SO2 = SO2;
    source.Sw = water.Sw;

    reaction.SH2_aCL = -j/(2*F*g.layerThickness(2));
    reaction.SO2_cCL = -j/(4*F*g.layerThickness(4));
    reaction.Sw_cCL = Mw*j/(2*F*g.layerThickness(4));
    reaction.H2Consumption_Area = -sum(SH2.*g.dx);
    reaction.O2Consumption_Area = -sum(SO2.*g.dx);
    reaction.WaterGeneration_Area = sum(water.Sw.*g.dx);
    reaction.H2Expected = j/(2*F);
    reaction.O2Expected = j/(4*F);
    reaction.WaterExpected = Mw*j/(2*F);
    reaction.errorH2 = reaction.H2Consumption_Area-reaction.H2Expected;
    reaction.errorO2 = reaction.O2Consumption_Area-reaction.O2Expected;
    reaction.errorWater = reaction.WaterGeneration_Area-reaction.WaterExpected;
    reaction.consistent = max(abs([reaction.errorH2,reaction.errorO2, ...
        reaction.errorWater])) < max(1e-14,1e-10*max(abs(j/F),1));


    %% 7. 输出诊断信息

    out.DH2 = DH2;
    out.DO2 = DO2;
    out.eps_g = water.eps_g;
    out.depsGdt = depsGdt;
    out.JH2 = JH2;
    out.JO2 = JO2;
    out.diffH2 = diffH2;
    out.diffO2 = diffO2;
    out.diffusionInfo.H2 = diffInfoH2;
    out.diffusionInfo.O2 = diffInfoO2;
    out.source = source;
    out.reaction = reaction;

end


%% ========================================================================
% 局部函数：连续气体区域上的一维有限体积扩散
% ========================================================================

function [divDiff,Jface,info] = regional_diffusion( ...
    phi,D,g,idx,leftType,leftValue,rightType,rightValue)

    phi = phi(:);
    D = D(:);
    idx = idx(:);
    N = numel(idx);
    Jface = zeros(N+1,1);

    for f = 2:N
        xFace = g.xf(idx(1)+f-1);
        dL = xFace-g.x(idx(f-1));
        dR = g.x(idx(f))-xFace;
        if D(f-1) == 0 || D(f) == 0
            Jface(f) = 0;
        else
            Jface(f) = -(phi(f)-phi(f-1))/(dL/D(f-1)+dR/D(f));
        end
    end

    if strcmp(leftType,'dirichlet') && D(1) ~= 0
        d = g.x(idx(1))-g.xf(idx(1));
        Jface(1) = -D(1)*(phi(1)-leftValue)/d;
    elseif strcmp(leftType,'neumann')
        Jface(1) = leftValue;
    end

    if strcmp(rightType,'dirichlet') && D(end) ~= 0
        d = g.xf(idx(end)+1)-g.x(idx(end));
        Jface(end) = -D(end)*(rightValue-phi(end))/d;
    elseif strcmp(rightType,'neumann')
        Jface(end) = rightValue;
    end

    divDiff = (Jface(1:end-1)-Jface(2:end))./g.dx(idx);
    info.Jleft = Jface(1);
    info.Jright = Jface(end);
    info.integralDiv = sum(divDiff.*g.dx(idx));
    info.boundaryNet = Jface(1)-Jface(end);
    info.balanceError = info.integralDiv-info.boundaryNet;
    info.conserved = abs(info.balanceError) < ...
        max(1e-12,1e-10*max(abs(info.boundaryNet),1));
end


%% ========================================================================
% 局部函数：Buck饱和蒸气压，用于动态孔隙率的温度导数
% ========================================================================

function psat = buck_psat(T)

    Tc = T-273.15;
    psat = zeros(size(T));
    warm = Tc >= 0;
    psat(warm) = 611.21.*exp((18.678-Tc(warm)./234.5).* ...
        Tc(warm)./(257.14+Tc(warm)));
    psat(~warm) = 611.15.*exp((23.036-Tc(~warm)./333.7).* ...
        Tc(~warm)./(279.82+Tc(~warm)));
end
