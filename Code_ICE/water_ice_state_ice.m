%-------------------------------------------------------------------------------
% @function: 计算总水守恒、冰守恒以及水在气/液/冰三相间的分配
% @author:   PJ, GPT
% @date:     20260923
% @input:    T->温度状态 [K]
%            mw->总水质量浓度，mw=mv+ml+mi [kg/m^3]
%            mi->四个多孔层中的冰质量浓度 [kg/m^3]
%            j->电流密度 [A/m^2]
%            kFreeze,kMelt->冻结/融化速率系数 [1/(K s)]
%            g,s->网格和状态编号
% @output:   dmw->总水状态导数 [kg/(m^3 s)]
%            dmi->冰状态导数 [kg/(m^3 s)]
%            out->水相、冰相、孔隙率、通量和守恒检查信息
%
% [ICE-2] 冰守恒：dmi/dt=Rfreeze-Rmelt；相变不改变总水mw。
% [ICE-3] 冰体积分数eps_i=mi/rho_i，并从原始孔隙率中直接扣除。
%-------------------------------------------------------------------------------

function [dmw,dmi,out] = water_ice_state_ice( ...
    T,mw,mi,j,kFreeze,kMelt,g,~)

    %% 1. 输入处理

    T = T(:);
    mw = max(mw(:),0);  % 仅处理求解器Newton迭代产生的微小负试探值
    miState = max(mi(:),0);

    R = 8.314;          % [J/(mol K)]
    F = 96485;          % [C/mol]
    Mw = 0.018;         % [kg/mol]
    rhoLiquid = 990;    % [kg/m^3]
    rhoIce = 920;       % [kg/m^3]
    freezingPoint = 273.15; % [K]

    idxPorous = g.idx_porous;
    eps0Layer = [0.8,0.3916,0,0.4207,0.8];
    eps0 = eps0Layer(g.layerId(:)).';
    eps0 = eps0(:);

    % 数值迭代时把冰限制在“现有总水”和“干孔隙容量”以内。
    % 正常积分结果由守恒方程保证满足该约束，projection仅用于试探状态。
    iceCapacity = rhoIce*eps0(idxPorous);
    miEval = min(miState,mw(idxPorous));
    miEval = min(miEval,iceCapacity);
    iceProjectionApplied = any(abs(miEval-miState) > ...
        10*eps(max(max(abs(miState)),1)));

    miFull = zeros(g.N,1);
    miFull(idxPorous) = miEval;


    %% 2. 膜含水量
    % 式(21)：lambda=EW*mw/(rho_pem*Mw)

    lambdaPEM = mw(g.idx_PEM)/(2150*Mw);
    dxPEM = g.dx(g.idx_PEM);
    lambdaMean = sum(lambdaPEM.*dxPEM)/sum(dxPEM);


    %% 3. 多孔区域中的气态水、液态水和冰
    % [ICE-2] 总水定义保持不变：mw=mv+ml+mi。
    % 先扣除冰得到可迁移水mMobile=mw-mi，再作气/液平衡。
    %
    % [ICE-3] 冰占据孔隙：epsAfterIce=eps0-mi/rhoIce。
    % 饱和时同时满足：
    %   mv=A*(epsAfterIce-ml/rhoLiquid)
    %   mMobile=mv+ml
    % 解得：ml=(mMobile-A*epsAfterIce)/(1-A/rhoLiquid)。

    Tc = T-273.15;
    psat = zeros(g.N,1);
    warm = Tc >= 0;
    psat(warm) = 611.21.*exp((18.678-Tc(warm)./234.5).* ...
        Tc(warm)./(257.14+Tc(warm)));
    psat(~warm) = 611.15.*exp((23.036-Tc(~warm)./333.7).* ...
        Tc(~warm)./(279.82+Tc(~warm)));
    A = Mw.*psat./(R.*T);

    eps_i = zeros(g.N,1);
    eps_i(idxPorous) = miEval/rhoIce;
    sIce = zeros(g.N,1);
    sIce(idxPorous) = eps_i(idxPorous)./eps0(idxPorous);

    mMobile = mw;
    mMobile(idxPorous) = max(mw(idxPorous)-miEval,0);
    mPore = zeros(g.N,1);
    mPore(idxPorous) = mw(idxPorous);
    epsAfterIce = eps0(idxPorous)-eps_i(idxPorous);

    mv = zeros(g.N,1);
    ml = zeros(g.N,1);
    eps_g = zeros(g.N,1);
    mvCapacity = A(idxPorous).*epsAfterIce;
    unsaturated = mMobile(idxPorous) <= mvCapacity;
    saturated = ~unsaturated;

    mvPorous = zeros(numel(idxPorous),1);
    mlPorous = zeros(numel(idxPorous),1);
    mvPorous(unsaturated) = mMobile(idxPorous(unsaturated));

    den = 1-A(idxPorous(saturated))/rhoLiquid;
    if any(den <= 0)
        error('water_ice_state_ice:InvalidPhaseEquilibrium', ...
            '气液平衡公式分母小于等于0。');
    end
    mlPorous(saturated) = (mMobile(idxPorous(saturated))- ...
        A(idxPorous(saturated)).*epsAfterIce(saturated))./den;
    mvPorous(saturated) = mMobile(idxPorous(saturated))-mlPorous(saturated);

    mv(idxPorous) = mvPorous;
    ml(idxPorous) = mlPorous;
    eps_g(idxPorous) = epsAfterIce-mlPorous/rhoLiquid;
    if any(eps_g(idxPorous) <= 0)
        error('water_ice_state_ice:BlockedPore', ...
            '液态水和冰已占满至少一个多孔控制体。');
    end


    %% 4. 冻结和融化动力学
    % [ICE-2] 第一版采用线性过冷/过热驱动力：
    %   Rfreeze=kFreeze*ml*(Tf-T)_+*(1-si)
    %   Rmelt  =kMelt*mi*(T-Tf)_+
    % 其中(1-si)使孔隙接近全冰时冻结速率自然衰减。

    remainingPoreFactor = max(1-sIce(idxPorous),0);
    Rfreeze = kFreeze.*mlPorous.*max(freezingPoint-T(idxPorous),0).* ...
        remainingPoreFactor;
    Rmelt = kMelt.*miEval.*max(T(idxPorous)-freezingPoint,0);
    dmi = Rfreeze-Rmelt;

    % 在状态边界上禁止继续越界；正常状态中这两项不改变动力学。
    atZeroIce = miState <= 0 & dmi < 0;
    atAllWaterFrozen = miState >= mw(idxPorous) & dmi > 0;
    atFullPore = miState >= iceCapacity & dmi > 0;
    dmi(atZeroIce | atAllWaterFrozen | atFullPore) = 0;

    dmiFull = zeros(g.N,1);
    dmiFull(idxPorous) = dmi;
    RfreezeFull = zeros(g.N,1);
    RmeltFull = zeros(g.N,1);
    RfreezeFull(idxPorous) = Rfreeze;
    RmeltFull(idxPorous) = Rmelt;


    %% 5. 分区域水扩散系数
    % 多孔层只有气态水和液态水能够迁移，因此扩散驱动量使用mw-mi；
    % 冰质量只通过独立的冰守恒方程变化，不随水扩散通量移动。

    idxA = [g.idx_aGDL;g.idx_aCL];
    idxC = [g.idx_cCL;g.idx_cGDL];
    Dw = zeros(g.N,1);
    Dw(idxA) = 8.69e-5.*(T(idxA)/298.15).^1.75.* ...
        eps_g(idxA).^1.5;
    Dw(idxC) = 2.48e-5.*(T(idxC)/298.15).^1.75.* ...
        eps_g(idxC).^1.5;

    lambdaPoly = 2.563-0.33.*lambdaPEM+0.0264.*lambdaPEM.^2- ...
        0.000671.*lambdaPEM.^3;
    Dw(g.idx_PEM) = 1e-10.*exp(2416.* ...
        (1/303.15-1./T(g.idx_PEM))).*lambdaPoly;
    if any(Dw < 0) || any(~isfinite(Dw))
        error('water_ice_state_ice:InvalidDiffusivity', ...
            '水扩散系数出现非法数值。');
    end

    transportWater = mw;
    transportWater(idxPorous) = mMobile(idxPorous);
    [JdiffA,diffInfoA] = regional_diffusion( ...
        transportWater(idxA),Dw(idxA),g,idxA,'dirichlet',0,'noflux',0);
    [JdiffM,diffInfoM] = regional_diffusion( ...
        transportWater(g.idx_PEM),Dw(g.idx_PEM),g,g.idx_PEM, ...
        'noflux',0,'noflux',0);
    [JdiffC,diffInfoC] = regional_diffusion( ...
        transportWater(idxC),Dw(idxC),g,idxC,'noflux',0,'dirichlet',0);

    Ndiff = zeros(g.N+1,1);
    Ndiff(idxA(1):idxA(end)+1) = JdiffA;
    Ndiff(g.idx_PEM(1):g.idx_PEM(end)+1) = JdiffM;
    Ndiff(idxC(1):idxC(end)+1) = JdiffC;


    %% 6. CL/PEM有限速率界面通量

    mSat = zeros(g.N,1);
    mSat(idxPorous) = eps_g(idxPorous).*A(idxPorous);
    waterActivity = zeros(g.N,1);
    validSat = mSat > 0;
    waterActivity(validSat) = min(max(mv(validSat)./mSat(validSat),0),1);
    waterActivity(ml > 0) = 1;

    lambdaScale = 22;
    thetaACL = waterActivity(g.idx_aCL(end));
    thetaCCL = waterActivity(g.idx_cCL(1));
    thetaMemLeft = min(max(lambdaPEM(1)/lambdaScale,0),1);
    thetaMemRight = min(max(lambdaPEM(end)/lambdaScale,0),1);

    CmemScale = lambdaScale*2150*Mw;
    Gleft = Dw(g.idx_PEM(1))*CmemScale/g.layerThickness(3);
    Gright = Dw(g.idx_PEM(end))*CmemScale/g.layerThickness(3);

    driveLeft = thetaACL-thetaMemLeft;
    driveRight = thetaMemRight-thetaCCL;

    % [FIX-1] 去除原先吸附为1、解吸为0.001的单向阀式修正。
    % 两个方向统一使用同一个有限速率界面导通关系，避免水只容易进入
    % PEM而几乎不能返回CL。后续若有界面动力学实验，再单独标定系数。
    factorLeft = 1.0;
    factorRight = 1.0;

    Ndiff(g.face_aCL_PEM) = Gleft*factorLeft*driveLeft;
    Ndiff(g.face_PEM_cCL) = Gright*factorRight*driveRight;


    %% 7. 质子电流和电渗拖曳水通量

    imCell = zeros(g.N,1);
    xiA = (g.x(g.idx_aCL)-g.layerBoundary(2))/g.layerThickness(2);
    xiC = (g.x(g.idx_cCL)-g.layerBoundary(4))/g.layerThickness(4);
    imCell(g.idx_aCL) = j.*min(max(xiA,0),1);
    imCell(g.idx_PEM) = j;
    imCell(g.idx_cCL) = j.*(1-min(max(xiC,0),1));

    imFace = zeros(g.N+1,1);
    imFace(g.face_aGDL_aCL:g.face_aCL_PEM) = ...
        linspace(0,j,g.face_aCL_PEM-g.face_aGDL_aCL+1).';
    imFace(g.face_aCL_PEM:g.face_PEM_cCL) = j;
    imFace(g.face_PEM_cCL:g.face_cCL_cGDL) = ...
        linspace(j,0,g.face_cCL_cGDL-g.face_PEM_cCL+1).';

    lambdaCell = zeros(g.N,1);
    lambdaCell(g.idx_aCL) = lambdaPEM(1);
    lambdaCell(g.idx_PEM) = lambdaPEM;
    lambdaCell(g.idx_cCL) = lambdaPEM(end);
    lambdaFace = zeros(g.N+1,1);
    lambdaFace(1) = lambdaCell(1);
    lambdaFace(end) = lambdaCell(end);
    lambdaFace(2:g.N) = 0.5*(lambdaCell(1:end-1)+lambdaCell(2:end));
    ndFace = (2.5/22).*lambdaFace;

    Neod = zeros(g.N+1,1);
    % [FIX-3] PEM内电渗拖曳必须覆盖两个CL|PEM边界面。
    % 原范围从face_aCL_PEM+1开始，遗漏了左侧aCL|PEM界面。
    activeEodFaces = g.face_aCL_PEM:g.face_PEM_cCL;
    Neod(activeEodFaces) = ndFace(activeEodFaces).*Mw.* ...
        imFace(activeEodFaces)/F;


    %% 8. 总水守恒方程
    % [ICE-2] 相变只在mv/ml/mi之间重新分配，不是总水的源项：
    %   dmw/dt=-div(Nw)+Sw

    Nw = Ndiff+Neod;
    Sw = zeros(g.N,1);
    Sw(g.idx_cCL) = Mw*j/(2*F*g.layerThickness(4));
    dmw = (Nw(1:end-1)-Nw(2:end))./g.dx + Sw;


    %% 9. 守恒检查和输出

    diffInfo = flux_balance(Ndiff,g.dx);
    diffInfo.anodeRegion = diffInfoA;
    diffInfo.membraneRegion = diffInfoM;
    diffInfo.cathodeRegion = diffInfoC;

    isFace = j-imFace;
    isCell = j-imCell;
    dimdx = zeros(g.N,1);
    dimdx(g.idx_aCL) = j/g.layerThickness(2);
    dimdx(g.idx_cCL) = -j/g.layerThickness(4);

    currentInfo.iv_a = j/g.layerThickness(2);
    currentInfo.iv_c = j/g.layerThickness(4);
    currentInfo.dimdx = dimdx;
    currentInfo.isFace = isFace;
    currentInfo.imFace = imFace;
    currentInfo.currentError = max(abs(isCell+imCell-j));
    currentInfo.faceCurrentError = max(abs(isFace+imFace-j));
    currentInfo.currentConserved = ...
        max(currentInfo.currentError,currentInfo.faceCurrentError) < ...
        max(1e-10,1e-12*max(abs(j),1));

    out.Nw = Nw;
    out.Ndiff = Ndiff;
    out.Neod = Neod;
    out.Dw = Dw;
    out.mv = mv;
    out.ml = ml;
    out.miFull = miFull;
    out.mMobile = mMobile;
    out.mPore = mPore;
    out.eps_g = eps_g;
    out.eps_i = eps_i;
    out.sIce = sIce;
    out.lambdaPEM = lambdaPEM;
    out.lambdaMean = lambdaMean;
    out.nd = (2.5/22)*lambdaMean;
    out.ndFace = ndFace;
    out.lambdaFace = lambdaFace;
    out.imCell = imCell;
    out.imFace = imFace;
    out.diffusionInfo = diffInfo;
    out.waterActivity = waterActivity;
    out.interfaceActivity.aCL = thetaACL;
    out.interfaceActivity.PEMLeft = thetaMemLeft;
    out.interfaceActivity.PEMRight = thetaMemRight;
    out.interfaceActivity.cCL = thetaCCL;
    out.interfaceConductance.aCL_PEM = Gleft;
    out.interfaceConductance.PEM_cCL = Gright;
    out.interfaceRateFactor.aCL_PEM = factorLeft;
    out.interfaceRateFactor.PEM_cCL = factorRight;
    out.interfaceFlux.aCL_PEM_diff = Ndiff(g.face_aCL_PEM);
    out.interfaceFlux.PEM_cCL_diff = Ndiff(g.face_PEM_cCL);
    out.interfaceFlux.aCL_PEM_total = Nw(g.face_aCL_PEM);
    out.interfaceFlux.PEM_cCL_total = Nw(g.face_PEM_cCL);
    out.interfaceModel = ...
        'conservative finite-rate CL/PEM exchange with regional diffusion + local EOD';
    out.currentInfo = currentInfo;
    out.Sw = Sw;
    out.dmiFull = dmiFull;
    out.Rfreeze = RfreezeFull;
    out.Rmelt = RmeltFull;
    out.kFreeze = kFreeze;
    out.kMelt = kMelt;
    out.rhoIce = rhoIce;
    out.freezingPoint = freezingPoint;
    out.iceProjectionApplied = iceProjectionApplied;
    out.phase.psat = psat;
    out.phase.A = A;
    out.phase.isSaturated = ml > 0;
    out.phase.iceModelEnabled = true;

    out.reaction.Sw_cCL = Mw*j/(2*F*g.layerThickness(4));
    out.reaction.WaterGeneration_Area = sum(Sw.*g.dx);
    out.reaction.WaterExpected = Mw*j/(2*F);
    out.reaction.errorWater = ...
        out.reaction.WaterGeneration_Area-out.reaction.WaterExpected;

end


%% ========================================================================
% 局部函数：连续区域上的一维有限体积扩散
% ========================================================================

function [Jface,info] = regional_diffusion(phi,D,g,idx, ...
    leftType,leftValue,rightType,rightValue)

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

    info = flux_balance(Jface,g.dx(idx));
end


%% ========================================================================
% 局部函数：检查面通量是否满足全局守恒
% ========================================================================

function info = flux_balance(Jface,dx)

    divJ = (Jface(1:end-1)-Jface(2:end))./dx;
    info.Jleft = Jface(1);
    info.Jright = Jface(end);
    info.integralDiv = sum(divJ.*dx);
    info.boundaryNet = Jface(1)-Jface(end);
    info.balanceError = info.integralDiv-info.boundaryNet;
    info.conserved = abs(info.balanceError) <= ...
        max(1e-12,1e-10*max(abs(info.boundaryNet),1));
end
