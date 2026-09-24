%-------------------------------------------------------------------------------
% @function: 计算总水守恒、cCL水合、冰守恒以及水在各相间的分配
% @author:   PJ, GPT
% @date:     20260923
% @input:    T->温度状态 [K]
%            mw->总水质量浓度；cCL中mw=mIon+mv+ml+mi [kg/m^3]
%            mi->四个多孔层中的冰质量浓度 [kg/m^3]
%            lambdaCCLState->cCL平均离聚物含水量 [-]
%            j->电流密度 [A/m^2]
%            kFreeze,kMelt->冻结/融化非平衡速率系数 [1/s]
%            tauHyd->cCL离聚物水合时间常数 [s]
%            hVaporCathode->cGDL/阴极气道水蒸气传质系数 [m/s]
%            g,s->网格和状态编号
% @output:   dmw->总水状态导数 [kg/(m^3 s)]
%            dmi->冰状态导数 [kg/(m^3 s)]
%            dlambdaCCL->cCL平均离聚物含水量导数 [1/s]
%            out->水相、冰相、孔隙率、通量和守恒检查信息
%
% [ICE-2] 冰守恒：dmi/dt=Rfreeze-Rmelt；相变不改变总水mw。
% [ICE-3] 冰体积分数eps_i=mi/rho_i，并从原始孔隙率中直接扣除。
% [HYD-1] 反应水先参与cCL离聚物水合；只有扣除结合水后的孔隙水
%   才能成为水蒸气、液水和冰。水合由局部水活度驱动，不使用时间开关。
% [HYD-2] 低温下cCL离聚物采用温度相关最大非冻结含水量lambdaSat(T)；
%   达到该上限后的产水才进入孔隙排水和结冰过程。
% [WATER-1] 水蒸气按Fick定律扩散，液态水按毛细压力梯度排出；不再把
%   mv+ml整体套用水蒸气扩散系数。
% [WATER-2] 阴极外边界由“水蒸气浓度恒为0”改为有限传质：
%   Nout=hVaporCathode*(mv_surface-mv_inlet)。附件1给出干空气入口，
%   因此mv_inlet=0；但有限h不再把整个阴极外侧当成无限快水汇。
%-------------------------------------------------------------------------------

function [dmw,dmi,dlambdaCCL,out] = water_ice_state_ice( ...
    T,mw,mi,lambdaCCLState,j,kFreeze,kMelt,tauHyd, ...
    hVaporCathode,g,s)

    %% 1. 输入处理

    T = T(:);
    mw = max(mw(:),0);  % 仅处理求解器Newton迭代产生的微小负试探值
    miState = max(mi(:),0);
    lambdaCCLStateRaw = lambdaCCLState;
    lambdaCCLState = max(lambdaCCLState,0);

    R = 8.314;          % [J/(mol K)]
    F = 96485;          % [C/mol]
    Mw = 0.018;         % [kg/mol]
    rhoLiquid = 990;    % [kg/m^3]
    rhoIce = 920;       % [kg/m^3]
    freezingPoint = 273.15; % [K]
    lambdaScale = 22;

    idxPorous = g.idx_porous;
    eps0Layer = [0.8,0.3916,0,0.4207,0.8];
    eps0 = eps0Layer(g.layerId(:)).';
    eps0 = eps0(:);

    % [HYD-1] cCL结合水不占孔隙，也不能冻结。附件1的CL离聚物含量
    % 为0.3，故单位lambda对应的结合水质量浓度如下。
    cCLIonomerWaterCoefficient = 0.3*2150*Mw/1.0;
    lambdaAvailable = min((mw(g.idx_cCL)-miState(s.local.mi_cCL))/ ...
        cCLIonomerWaterCoefficient);
    lambdaUpper = min(lambdaScale,max(lambdaAvailable,0));
    lambdaCCLEval = min(lambdaCCLState,lambdaUpper);

    mIonomer = zeros(g.N,1);
    mIonomer(g.idx_cCL) = ...
        cCLIonomerWaterCoefficient*lambdaCCLEval;
    mPore = zeros(g.N,1);
    mPore(idxPorous) = max(mw(idxPorous)-mIonomer(idxPorous),0);

    % 数值迭代时把冰限制在“现有孔隙水”和“干孔隙容量”以内。
    % 正常积分结果由守恒方程保证满足该约束，projection仅用于试探状态。
    iceCapacity = rhoIce*eps0(idxPorous);
    miEval = min(miState,mPore(idxPorous));
    miEval = min(miEval,iceCapacity);
    iceProjectionApplied = any(abs(miEval-miState) > ...
        10*eps(max(max(abs(miState)),1)));

    miFull = zeros(g.N,1);
    miFull(idxPorous) = miEval;


    %% 2. 膜含水量
    % 式(21)：lambda=EW*mw/(rho_pem*Mw)

    lambdaPEMRaw = mw(g.idx_PEM)/(2150*Mw);
    % 本构关系的适用范围取0<=lambda<=22；总水状态本身仍保持守恒，
    % 这里仅限制Newton试探值进入扩散、电导和拖曳公式时的取值。
    lambdaPEM = min(max(lambdaPEMRaw,0),lambdaScale);
    dxPEM = g.dx(g.idx_PEM);
    lambdaMean = sum(lambdaPEM.*dxPEM)/sum(dxPEM);


    %% 3. 多孔区域中的气态水、液态水和冰
    % [ICE-2][HYD-1] 非cCL区域mw=mv+ml+mi；cCL区域还要先扣除
    % 离聚物结合水mIonomer，再以mMobile=mPore-mi作气/液平衡。
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

    mMobile = zeros(g.N,1);
    mMobile(idxPorous) = max(mPore(idxPorous)-miEval,0);
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


    %% 4. cCL离聚物水合动力学

    mSat = zeros(g.N,1);
    mSat(idxPorous) = eps_g(idxPorous).*A(idxPorous);
    waterActivity = zeros(g.N,1);
    validSat = mSat > 0;
    waterActivity(validSat) = min(max(mv(validSat)./mSat(validSat),0),1);
    waterActivity(ml > 0) = 1;

    thetaCCLPore = sum(waterActivity(g.idx_cCL).*g.dx(g.idx_cCL))/ ...
        g.layerThickness(4);
    thetaMemRight = min(max(lambdaPEM(end)/lambdaScale,0),1);
    hydrationActivity = max(thetaCCLPore,thetaMemRight);
    lambdaCCLEquilibrium = lambdaScale*hydrationActivity;

    % [HYD-2] 采用附件给出的Nafion低温最大非冻结膜水含量关系。
    TCCLMean = sum(T(g.idx_cCL).*g.dx(g.idx_cCL))/g.layerThickness(4);
    if TCCLMean < 223.15
        lambdaSaturationCCL = 4.837;
    elseif TCCLMean < freezingPoint
        lambdaSaturationCCL = 1/(-1.304+0.01479*TCCLMean- ...
            3.594e-5*TCCLMean^2);
    else
        lambdaSaturationCCL = lambdaScale;
    end
    lambdaSaturationCCL = min(max(lambdaSaturationCCL,0),lambdaScale);

    % [HYD-1] 反应生成水在cCL内原位产生，应先有机会进入离聚物，
    % 不能先被极快的气相扩散清空后再用瞬时孔隙水反推水合量。这里用
    % 反应产水率限制水合速率，保证结合水增长不会超过实际水源；当已有
    % 孔隙水时，也允许其按tauHyd被吸收。剩余产水才进入孔隙排水/冻结。
    reactionWaterRateCCL = Mw*j/(2*F*g.layerThickness(4));
    if j > 0
        lambdaCCLTarget = lambdaSaturationCCL;
    else
        lambdaCCLTarget = min([lambdaCCLEquilibrium,lambdaUpper, ...
            lambdaSaturationCCL]);
    end
    hydrationKineticRate = ...
        (lambdaCCLTarget-lambdaCCLStateRaw)/tauHyd;
    existingPoreSupplyRate = max(lambdaUpper-lambdaCCLStateRaw,0)/tauHyd;
    productionSupplyRate = reactionWaterRateCCL/ ...
        cCLIonomerWaterCoefficient;
    if hydrationKineticRate > 0
        dlambdaCCL = min(hydrationKineticRate, ...
            existingPoreSupplyRate+productionSupplyRate);
    else
        dlambdaCCL = hydrationKineticRate;
    end
    if lambdaCCLStateRaw <= 0 && dlambdaCCL < 0
        dlambdaCCL = 0;
    elseif lambdaCCLStateRaw >= lambdaUpper && dlambdaCCL > 0
        dlambdaCCL = 0;
    end


    %% 5. 冻结和融化动力学
    % [ICE-2] 按文献中的分段相变关系：
    %   Rfreeze = kFreeze*rhoLiquid*eps0*epsLiquid, T<=Tf
    %   Rmelt   = kMelt*rhoIce*epsIce,          T>Tf
    %
    % eps0*epsLiquid表示液态水占控制体总体积的比例；epsIce同样
    % 使用控制体总体积为基准。由于ml和mi均为单位控制体体积内的
    % 质量浓度，故rhoLiquid*eps0*epsLiquid=ml，rhoIce*epsIce=mi。
    % 该关系不再显式乘过冷度、过热度或剩余孔隙因子。

    liquidSaturation = mlPorous./(rhoLiquid*eps0(idxPorous));
    liquidBulkVolumeFraction = eps0(idxPorous).*liquidSaturation;
    iceBulkVolumeFraction = eps_i(idxPorous);

    freezeRegion = T(idxPorous) <= freezingPoint;
    meltRegion = ~freezeRegion;
    Rfreeze = zeros(numel(idxPorous),1);
    Rmelt = zeros(numel(idxPorous),1);
    Rfreeze(freezeRegion) = kFreeze*rhoLiquid.* ...
        liquidBulkVolumeFraction(freezeRegion);
    Rmelt(meltRegion) = kMelt*rhoIce.* ...
        iceBulkVolumeFraction(meltRegion);
    dmi = Rfreeze-Rmelt;

    % 在状态边界上禁止继续越界；正常状态中这两项不改变动力学。
    atZeroIce = miState <= 0 & dmi < 0;
    atAllWaterFrozen = miState >= mPore(idxPorous) & dmi > 0;
    atFullPore = miState >= iceCapacity & dmi > 0;
    dmi(atZeroIce | atAllWaterFrozen | atFullPore) = 0;

    dmiFull = zeros(g.N,1);
    dmiFull(idxPorous) = dmi;
    RfreezeFull = zeros(g.N,1);
    RmeltFull = zeros(g.N,1);
    RfreezeFull(idxPorous) = Rfreeze;
    RmeltFull(idxPorous) = Rmelt;


    %% 6. 水蒸气扩散系数
    % [WATER-1] 多孔层扩散驱动量只使用mv。液水由下一节的毛细通量
    % 单独输运，冰与离聚物结合水均不随水蒸气扩散移动。

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
    invalidDw = real(Dw) < 0 | ~isfinite(Dw) | abs(imag(Dw)) > 0;
    if any(invalidDw)
        badCell = find(invalidDw,1);
        error('water_ice_state_ice:InvalidDiffusivity', ...
            ['水扩散系数出现非法数值：cell=%d, T=%.6g K, ', ...
            'Dw=%g, lambdaPEMRaw范围=[%.6g,%.6g]。'], ...
            badCell,T(badCell),real(Dw(badCell)), ...
            min(lambdaPEMRaw),max(lambdaPEMRaw));
    end

    transportWater = mw;
    transportWater(idxPorous) = mv(idxPorous);
    [JdiffA,diffInfoA] = regional_diffusion( ...
        transportWater(idxA),Dw(idxA),g,idxA,'dirichlet',0,'noflux',0);
    [JdiffM,diffInfoM] = regional_diffusion( ...
        transportWater(g.idx_PEM),Dw(g.idx_PEM),g,g.idx_PEM, ...
        'noflux',0,'noflux',0);
    % [WATER-2] 干空气是气道外部浓度为0，不等于cGDL边界
    % 自身浓度始终为0。边界通量由有限传质速率决定。
    cathodeInletVaporConcentration = 0;
    [JdiffC,diffInfoC] = regional_diffusion( ...
        transportWater(idxC),Dw(idxC),g,idxC,'noflux',0, ...
        'mass_transfer',[hVaporCathode,cathodeInletVaporConcentration]);

    Ndiff = zeros(g.N+1,1);
    Ndiff(idxA(1):idxA(end)+1) = JdiffA;
    Ndiff(g.idx_PEM(1):g.idx_PEM(end)+1) = JdiffM;
    Ndiff(idxC(1):idxC(end)+1) = JdiffC;


    %% 7. CL/PEM有限速率界面通量
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
    % 两个方向统一使用同一个有限速率界面导通关系。
    % 当前30为待核对值，不作为新的独立标定参数。
    factorLeft = 30;
    factorRight = 30;

    Ndiff(g.face_aCL_PEM) = Gleft*factorLeft*driveLeft;
    Ndiff(g.face_PEM_cCL) = Gright*factorRight*driveRight;


    %% 8. 液态水毛细排水

    % [WATER-1] Darcy/Leverett关系：
    %   Nl=-rhoL*Keff*krl/muL*grad(pc)
    %   pc=-sigma*cos(theta)*sqrt(eps0/K0)*J(sL)
    % 其中K0和接触角直接采用附件1中GDL/CL的数值。
    liquidSaturationCapillary = zeros(g.N,1);
    liquidSaturationCapillary(idxPorous) = min(max(mlPorous./ ...
        (rhoLiquid*max(epsAfterIce,1e-12)),0),1);

    intrinsicPermeabilityLayer = [6.2e-12,6.2e-13,0,6.2e-13,6.2e-12];
    contactAngleLayer = [110,100,0,100,110];
    K0 = intrinsicPermeabilityLayer(g.layerId(:)).';
    K0 = K0(:);
    contactAngle = contactAngleLayer(g.layerId(:)).';
    contactAngle = contactAngle(:);

    Keff = zeros(g.N,1);
    relativePermeability = zeros(g.N,1);
    capillaryPressure = zeros(g.N,1);
    liquidMobility = zeros(g.N,1);
    porosityRatio = max((eps0(idxPorous)-eps_i(idxPorous))./ ...
        eps0(idxPorous),0);
    Keff(idxPorous) = K0(idxPorous).*porosityRatio.^3;
    relativePermeability(idxPorous) = ...
        liquidSaturationCapillary(idxPorous).^3;
    leverett = 1.417*liquidSaturationCapillary(idxPorous)- ...
        2.12*liquidSaturationCapillary(idxPorous).^2+ ...
        1.263*liquidSaturationCapillary(idxPorous).^3;
    capillaryPressure(idxPorous) = -0.075.*cosd(contactAngle(idxPorous)).* ...
        sqrt(eps0(idxPorous)./K0(idxPorous)).*leverett;
    liquidViscosity = 2.414e-5.*10.^(247.8./(T-140));
    liquidMobility(idxPorous) = rhoLiquid.*Keff(idxPorous).* ...
        relativePermeability(idxPorous)./liquidViscosity(idxPorous);

    [JliqA,capillaryInfoA] = regional_capillary_flux( ...
        capillaryPressure(idxA),liquidMobility(idxA),g,idxA, ...
        'dirichlet',0,'noflux',0);
    [JliqC,capillaryInfoC] = regional_capillary_flux( ...
        capillaryPressure(idxC),liquidMobility(idxC),g,idxC, ...
        'noflux',0,'dirichlet',0);
    Nliq = zeros(g.N+1,1);
    Nliq(idxA(1):idxA(end)+1) = JliqA;
    Nliq(idxC(1):idxC(end)+1) = JliqC;


    %% 9. 质子电流和电渗拖曳水通量

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
    lambdaCell(g.idx_cCL) = lambdaCCLEval;
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


    %% 10. 总水守恒方程
    % [ICE-2] 相变只在mv/ml/mi之间重新分配，不是总水的源项：
    %   dmw/dt=-div(Nw)+Sw

    Nw = Ndiff+Nliq+Neod;
    Sw = zeros(g.N,1);
    Sw(g.idx_cCL) = Mw*j/(2*F*g.layerThickness(4));
    dmw = (Nw(1:end-1)-Nw(2:end))./g.dx + Sw;

    dmIonomer = zeros(g.N,1);
    dmIonomer(g.idx_cCL) = ...
        cCLIonomerWaterCoefficient*dlambdaCCL;
    dmPore = dmw-dmIonomer;


    %% 11. 守恒检查和输出

    diffInfo = flux_balance(Ndiff,g.dx);
    diffInfo.anodeRegion = diffInfoA;
    diffInfo.membraneRegion = diffInfoM;
    diffInfo.cathodeRegion = diffInfoC;
    capillaryInfo = flux_balance(Nliq,g.dx);
    capillaryInfo.anodeRegion = capillaryInfoA;
    capillaryInfo.cathodeRegion = capillaryInfoC;

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
    out.Nliq = Nliq;
    out.Neod = Neod;
    out.Dw = Dw;
    out.mv = mv;
    out.ml = ml;
    out.miFull = miFull;
    out.mMobile = mMobile;
    out.mPore = mPore;
    out.mIonomer = mIonomer;
    out.dmIonomer = dmIonomer;
    out.dmPore = dmPore;
    out.eps_g = eps_g;
    out.eps_i = eps_i;
    % [ICE-3a] 题目式(9)的冰体积分数，以控制体总体积为基准。
    out.iceVolumeFraction = eps_i;
    % sIce为孔隙内冰饱和度，只用于堵孔及有效面积修正。
    out.sIce = sIce;
    out.lambdaPEM = lambdaPEM;
    out.lambdaPEMRaw = lambdaPEMRaw;
    out.lambdaMean = lambdaMean;
    out.lambdaCCL = lambdaCCLEval;
    out.lambdaCCLState = lambdaCCLStateRaw;
    out.lambdaCCLEquilibrium = lambdaCCLEquilibrium;
    out.lambdaCCLTarget = lambdaCCLTarget;
    out.lambdaSaturationCCL = lambdaSaturationCCL;
    out.dlambdaCCL = dlambdaCCL;
    out.tauHyd = tauHyd;
    out.cCLIonomerWaterCoefficient = cCLIonomerWaterCoefficient;
    out.hydrationActivity = hydrationActivity;
    out.hydrationKineticRate = hydrationKineticRate;
    out.hydrationProductionSupplyRate = productionSupplyRate;
    out.hydrationExistingPoreSupplyRate = existingPoreSupplyRate;
    out.nd = (2.5/22)*lambdaMean;
    out.ndFace = ndFace;
    out.lambdaFace = lambdaFace;
    out.imCell = imCell;
    out.imFace = imFace;
    out.diffusionInfo = diffInfo;
    out.capillaryInfo = capillaryInfo;
    out.capillaryPressure = capillaryPressure;
    out.liquidMobility = liquidMobility;
    out.liquidSaturationCapillary = liquidSaturationCapillary;
    out.relativePermeability = relativePermeability;
    out.effectiveLiquidPermeability = Keff;
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
        'finite-rate CL/PEM exchange + vapor diffusion + capillary liquid drainage + local EOD';
    out.vaporBoundary.model = ...
        'finite mass transfer from cGDL to dry cathode channel';
    out.vaporBoundary.massTransferCoefficient = hVaporCathode;
    out.vaporBoundary.externalConcentration = ...
        cathodeInletVaporConcentration;
    out.vaporBoundary.outwardFlux = JdiffC(end);
    out.currentInfo = currentInfo;
    out.Sw = Sw;
    out.dmiFull = dmiFull;
    out.Rfreeze = RfreezeFull;
    out.Rmelt = RmeltFull;
    out.liquidSaturation = liquidSaturation;
    out.liquidBulkVolumeFraction = liquidBulkVolumeFraction;
    out.iceBulkVolumeFraction = iceBulkVolumeFraction;
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
% 局部函数：液态水毛细通量
% 干湿前沿若使用扩散型调和平均，干侧零迁移率会把液水永久锁死。
% Darcy通量按压力梯度定方向，并使用上游（流出侧）迁移率。
% ========================================================================

function [Jface,info] = regional_capillary_flux(pc,mobility,g,idx, ...
    leftType,leftValue,rightType,rightValue)

    pc = pc(:);
    mobility = mobility(:);
    idx = idx(:);
    N = numel(idx);
    Jface = zeros(N+1,1);

    for f = 2:N
        gradient = (pc(f)-pc(f-1))/(g.x(idx(f))-g.x(idx(f-1)));
        if gradient <= 0
            mobilityFace = mobility(f-1);
        else
            mobilityFace = mobility(f);
        end
        Jface(f) = -mobilityFace*gradient;
    end

    if strcmp(leftType,'dirichlet')
        gradient = (pc(1)-leftValue)/(g.x(idx(1))-g.xf(idx(1)));
        if gradient > 0
            mobilityFace = mobility(1);
        else
            mobilityFace = 0;
        end
        Jface(1) = -mobilityFace*gradient;
    elseif strcmp(leftType,'neumann') || strcmp(leftType,'noflux')
        Jface(1) = leftValue;
    end

    if strcmp(rightType,'dirichlet')
        gradient = (rightValue-pc(end))/(g.xf(idx(end)+1)-g.x(idx(end)));
        if gradient < 0
            mobilityFace = mobility(end);
        else
            mobilityFace = 0;
        end
        Jface(end) = -mobilityFace*gradient;
    elseif strcmp(rightType,'neumann') || strcmp(rightType,'noflux')
        Jface(end) = rightValue;
    end

    info = flux_balance(Jface,g.dx(idx));
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
    elseif strcmp(leftType,'neumann') || strcmp(leftType,'noflux')
        Jface(1) = leftValue;
    end

    if strcmp(rightType,'dirichlet') && D(end) ~= 0
        d = g.xf(idx(end)+1)-g.x(idx(end));
        Jface(end) = -D(end)*(rightValue-phi(end))/d;
    elseif strcmp(rightType,'mass_transfer')
        % 右边界正方向为流出计算区域。rightValue=[h,phiExternal]。
        h = rightValue(1);
        phiExternal = rightValue(2);
        Jface(end) = h*(phi(end)-phiExternal);
    elseif strcmp(rightType,'neumann') || strcmp(rightType,'noflux')
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
