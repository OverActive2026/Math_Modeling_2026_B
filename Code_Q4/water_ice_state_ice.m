%-------------------------------------------------------------------------------
% @function: 计算总水守恒、cCL水合、冰守恒以及水在各相间的分配
% @author:   PJ, GPT
% @date:     20260923
% @input:    T->温度状态 [K]
%            mw->总水质量浓度；cCL中mw=mIon+mv+ml+mi [kg/m^3]
%            mi->四个多孔层中的冰质量浓度 [kg/m^3]
%            lambdaCCLState->cCL层平均离聚物含水量 [-]
%            j->电流密度 [A/m^2]
%            kFreeze,kMelt->冻结/融化非平衡速率系数 [1/s]
%            phaseTransitionWidthK->0 degC相变开关平滑半宽 [K]
%            directFreezeFraction->剩余反应产水的低温直接冻结比例
%            freezeNucleationLambdaFraction->直接冻结开始的lambda/lambdaSat阈值
%            tauHyd->cCL离聚物水合时间常数 [s]
%            hVaporCathode,hVaporAnode->两侧气道水蒸气传质系数 [m/s]
%            interfaceRateFactor->CL/PEM有限速率界面倍率 [-]
%            g,s->网格和状态编号
% @output:   dmw->总水状态导数 [kg/(m^3 s)]
%            dmi->冰状态导数 [kg/(m^3 s)]
%            dlambdaCCL->cCL层平均离聚物含水量导数 [1/s]
%            out->水相、冰相、孔隙率、通量和守恒检查信息
%
% [ICE-2] 冰守恒：dmi/dt=Rfreeze-Rmelt；相变不改变总水mw。
% [ICE-3] 冰体积分数eps_i=mi/rho_i，并从原始孔隙率中直接扣除。
% [PHASE-1] CL/PEM两侧统一使用离聚物归一化水活度作为驱动力。
% [PHASE-2] 阳极与阴极外边界都采用有限传质，不再使用零浓度理想水汇。
% [PHASE-3] cCL水合速率不得超过现有孔隙水与可分配反应产水，
%   并显式为低温冻结保留directFreezeFraction份额。
% [PHASE-4] 低温下未被离聚物吸收的部分反应产水可直接进入冰相，
%   表示催化层内未解析的局部成核/短寿命超冷液水。
% [HYD-2] 低温下cCL离聚物采用温度相关最大非冻结含水量lambdaSat(T)；
%   达到该上限后的产水才进入孔隙排水和结冰过程。
% [WATER-1] 水蒸气按Fick定律扩散，液态水按毛细压力梯度排出；不再把
%   mv+ml整体套用水蒸气扩散系数。
% [WATER-2] 阴极外边界由“水蒸气浓度恒为0”改为有限传质：
%   Nout=hVaporCathode*(mv_surface-mv_inlet)。附件1给出干空气入口，
%   因此mv_inlet=0；但有限h不再把整个阴极外侧当成无限快水汇。
%-------------------------------------------------------------------------------

function [dmw,dmi,dlambdaCCL,out] = water_ice_state_ice( ...
    T,mw,mi,lambdaCCLState,j,kFreeze,kMelt,phaseTransitionWidthK, ...
    directFreezeFraction, ...
    freezeNucleationLambdaFraction,tauHyd,hVaporCathode,hVaporAnode, ...
    interfaceRateFactor,g,s)

    %% 1. 输入处理

    T = T(:);
    mw = max(mw(:),0);  % 仅处理求解器Newton迭代产生的微小负试探值
    miState = max(mi(:),0);
    lambdaCCLStateRaw = lambdaCCLState(:);
    if numel(lambdaCCLStateRaw) ~= 1
        error('water_ice_state_ice:InvalidCCLLambdaSize', ...
            'lambdaCCLState必须是cCL层平均标量。');
    end
    lambdaCCLState = max(lambdaCCLStateRaw,0);

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
    % 稳定版保留一个cCL层平均lambda状态。上限取整层各
    % 控制体可用总水量的最小值，避免任一网格出现负孔隙水。
    lambdaAvailable = min((mw(g.idx_cCL)- ...
        miState(s.local.mi_cCL))/cCLIonomerWaterCoefficient);
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
    % 粗步长Newton试探可能暂时越过全堵边界。本构关系中只将
    % 气相孔隙率投影到极小正数，使求解器可拒绝该试探步；
    % 接受解是否真正堵塞仍由最终健康检查/事件判定。
    poreProjectionApplied = any(eps_g(idxPorous) <= 0);
    eps_g(idxPorous) = max(eps_g(idxPorous),1e-12);


    %% 4. cCL离聚物水合动力学

    mSat = zeros(g.N,1);
    mSat(idxPorous) = eps_g(idxPorous).*A(idxPorous);
    waterActivity = zeros(g.N,1);
    validSat = mSat > 0;
    waterActivity(validSat) = min(max(mv(validSat)./mSat(validSat),0),1);
    waterActivity(ml > 0) = 1;

    % 水合动力学保留经验证稳定的层平均写法；瞬时水活度
    % 只在停止反应时决定平衡脱水，反应时的实际吸水量由下方
    % 的水源速率限制。
    dxCCL = g.dx(g.idx_cCL);
    thetaCCLPoreCell = waterActivity(g.idx_cCL);
    thetaCCLPore = sum(thetaCCLPoreCell.*dxCCL)/sum(dxCCL);
    thetaMemRightHydration = min(max(lambdaPEM(end)/lambdaScale,0),1);
    hydrationActivity = max(thetaCCLPore,thetaMemRightHydration);
    lambdaCCLEquilibrium = lambdaScale*hydrationActivity;

    TCCLMean = sum(T(g.idx_cCL).*dxCCL)/sum(dxCCL);
    lambdaSaturationCCL = nonfreezing_lambda_saturation( ...
        TCCLMean,lambdaScale,freezingPoint);
    lambdaSaturationCCLCell = nonfreezing_lambda_saturation( ...
        T(g.idx_cCL),lambdaScale,freezingPoint);
    lambdaSaturationPEM = nonfreezing_lambda_saturation( ...
        T(g.idx_PEM),lambdaScale,freezingPoint);
    % 反应时局部生成水使cCL离聚物朝低温非冻结上限水合；
    % 停止反应后才回到由实际孔隙/PEM水活度给出的平衡值。
    % 这里的target只是动力学驱动，并不代表可以无条件吸水；
    % 真实水合速率仍由下方的孔隙水+反应产水上限约束。
    if j > 0
        lambdaCCLTarget = lambdaSaturationCCL;
    else
        lambdaCCLTarget = min(lambdaCCLEquilibrium,lambdaSaturationCCL);
    end

    reactionWaterRateCCL = Mw*j/(2*F*g.layerThickness(4));
    hydrationKineticRate = ...
        (lambdaCCLTarget-lambdaCCLStateRaw)/tauHyd;
    existingPoreSupplyRate = ...
        max(lambdaUpper-lambdaCCLStateRaw,0)/tauHyd;
    % [PHASE-3][PHASE-4] 最多仅把(1-directFreezeFraction)的当时
    % 反应产水用于离聚物水合，为低温局部冻结保留可辨识
    % 的水源。已有孔隙水仍可按tauHyd被吸收。
    productionSupplyRate = (1-directFreezeFraction)* ...
        reactionWaterRateCCL/cCLIonomerWaterCoefficient;
    dlambdaCCL = hydrationKineticRate;
    if hydrationKineticRate > 0
        dlambdaCCL = min(hydrationKineticRate, ...
            existingPoreSupplyRate+productionSupplyRate);
    end
    if lambdaCCLStateRaw <= 0 && dlambdaCCL < 0
        dlambdaCCL = 0;
    end
    % [NUM-HYD] 不再在lambdaCCL=lambdaUpper处把正水合速率突然置零。
    % 上方productionSupplyRate已经保证离聚物吸收量不超过可分配的
    % 当时反应产水，因此反应产水与结合水增长可以在同一微分时刻发生。
    % 原先“先产水、下一时刻再吸收”的硬开关会使Newton迭代在约束面
    % 两侧反复跳变，是第一题开局极慢的主要数值原因。

    dmIonomer = zeros(g.N,1);
    dmIonomer(g.idx_cCL) = ...
        cCLIonomerWaterCoefficient*dlambdaCCL;
    % 未被离聚物当即吸收的反应产水才能进入孔隙/冰相。
    freezableProductionRateCCL = max(reactionWaterRateCCL- ...
        max(cCLIonomerWaterCoefficient*dlambdaCCL,0),0);
    lambdaNonfreezingFraction = min(max(lambdaCCLEval/ ...
        lambdaSaturationCCL,0),1);
    nucleationActivation = min(max((lambdaNonfreezingFraction- ...
        freezeNucleationLambdaFraction)/ ...
        max(1-freezeNucleationLambdaFraction,1e-6),0),1);

    lambdaCCL = lambdaCCLEval;
    lambdaCCLEquilibriumCell = lambdaCCLEquilibrium* ...
        ones(numel(g.idx_cCL),1);
    lambdaCCLTargetCell = lambdaCCLTarget*ones(numel(g.idx_cCL),1);
    hydrationActivityCell = hydrationActivity*ones(numel(g.idx_cCL),1);


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

    % [NUM-PHASE] width=0严格保留题目分段关系；width>0时只在冰点
    % 邻域内用光滑权重替代硬开关，避免求解器在冻结/融化两套方程间
    % 反复切换并把内部步长压到极小。远离冰点时两者数值等价。
    if phaseTransitionWidthK > 0
        freezeWeight = 0.5*(1-tanh((T(idxPorous)-freezingPoint)./ ...
            phaseTransitionWidthK));
    else
        freezeWeight = double(T(idxPorous) <= freezingPoint);
    end
    meltWeight = 1-freezeWeight;
    RfreezeLiquid = kFreeze*rhoLiquid.*liquidBulkVolumeFraction.* ...
        freezeWeight;
    RfreezeDirect = zeros(numel(idxPorous),1);
    Rmelt = kMelt*rhoIce.*iceBulkVolumeFraction.*meltWeight;
    % [PHASE-4] 仅在低于c冰点的cCL内，将剩余反应产水的
    % 可调比例直接分配到冰相。这不增加总水，只改变相分配。
    directFreezeSupply = min(directFreezeFraction* ...
        reactionWaterRateCCL,freezableProductionRateCCL);
    % 直接冻结也需要极少量局部可移动水作为成核/薄膜库。
    % 用0.1%初始孔隙液水容量作平滑过渡尺度，避免在
    % mi=mPore边界上“先产水/先冻结”的高频切换。
    mobileWaterTransition = max(1e-3*rhoLiquid* ...
        eps0(g.idx_cCL),1e-6);
    mobileWaterActivation = mMobile(g.idx_cCL)./ ...
        (mMobile(g.idx_cCL)+mobileWaterTransition);
    RfreezeDirect(s.local.mi_cCL) = nucleationActivation.* ...
        directFreezeSupply.*mobileWaterActivation.* ...
        freezeWeight(s.local.mi_cCL);
    Rfreeze = RfreezeLiquid+RfreezeDirect;
    dmi = Rfreeze-Rmelt;

    % 状态边界上禁止继续越界。当全部孔隙水已结冰时，
    % 新产生的水先进入总水状态，随后时刻再发生相变。
    atZeroIce = miState <= 0 & dmi < 0;
    atAllWaterFrozen = miState >= mPore(idxPorous) & dmi > 0;
    atFullPore = miState >= iceCapacity & dmi > 0;
    dmi(atZeroIce | atAllWaterFrozen | atFullPore) = 0;

    % 受状态边界限幅后的净相变率用于潜热计算。
    Rfreeze = max(dmi+Rmelt,0);
    dmiFull = zeros(g.N,1);
    dmiFull(idxPorous) = dmi;
    RfreezeFull = zeros(g.N,1);
    RmeltFull = zeros(g.N,1);
    RfreezeFull(idxPorous) = Rfreeze;
    RmeltFull(idxPorous) = Rmelt;
    RfreezeLiquidFull = zeros(g.N,1);
    RfreezeDirectFull = zeros(g.N,1);
    freezeWeightFull = zeros(g.N,1);
    RfreezeLiquidFull(idxPorous) = RfreezeLiquid;
    RfreezeDirectFull(idxPorous) = RfreezeDirect;
    freezeWeightFull(idxPorous) = freezeWeight;


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
    % [PHASE-2] 两侧干气都是外部浓度为0的有限传质边界。
    anodeInletVaporConcentration = 0;
    [JdiffA,diffInfoA] = regional_diffusion( ...
        transportWater(idxA),Dw(idxA),g,idxA, ...
        'mass_transfer',[hVaporAnode,anodeInletVaporConcentration], ...
        'noflux',0);
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
    thetaCCLIonomer = min(max(lambdaCCLEval(1)/ ...
        lambdaSaturationCCLCell(1),0),1);
    thetaMemLeft = min(max(lambdaPEM(1)/lambdaSaturationPEM(1),0),1);
    thetaMemRight = min(max(lambdaPEM(end)/lambdaSaturationPEM(end),0),1);

    CmemScale = lambdaScale*2150*Mw;
    Gleft = Dw(g.idx_PEM(1))*CmemScale/g.layerThickness(3);
    Gright = Dw(g.idx_PEM(end))*CmemScale/g.layerThickness(3);

    driveLeft = thetaACL-thetaMemLeft;
    % [PHASE-1] 右界面改为PEM离聚物与cCL离聚物之间的
    % 归一化水活度差，不再把cCL孔隙蒸气活度直接与PEM比较。
    driveRight = thetaMemRight-thetaCCLIonomer;

    % [PHASE-1b] 界面导通率对通量方向保持完全对称，但随界面
    % 两侧较湿一侧的水活度增大。这表示水合后离聚物的界面
    % 水传递能力上升，而不是“吸收/解吸”取不同倍率的单向阀。
    % 0.01*(1+29*a^4)使干态界面接近基准倍率，充分水合时
    % 平滑恢复到约0.3倍有效导通，在后期回扩散与刚性之间折中。
    wetnessLeft = max(thetaACL,thetaMemLeft);
    wetnessRight = max(thetaMemRight,thetaCCLIonomer);
    factorLeft = interfaceRateFactor*(1+29*wetnessLeft^4);
    factorRight = interfaceRateFactor*(1+29*wetnessRight^4);

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
    % [PORE-STOP] 水+冰对原始孔隙的总占用率。该量与冰填充率
    % 分开：即使冰未占满孔隙，剩余孔隙也可能被液水填满。
    poreOccupancyFraction = zeros(g.N,1);
    poreOccupancyFraction(idxPorous) = (eps_i(idxPorous)+ ...
        ml(idxPorous)/rhoLiquid)./eps0(idxPorous);
    out.poreOccupancyFraction = poreOccupancyFraction;
    out.lambdaPEM = lambdaPEM;
    out.lambdaPEMRaw = lambdaPEMRaw;
    out.lambdaMean = lambdaMean;
    out.lambdaCCL = lambdaCCL;
    out.lambdaCCLCell = lambdaCCLEval*ones(numel(g.idx_cCL),1);
    out.lambdaCCLState = lambdaCCLStateRaw;
    out.lambdaCCLEquilibrium = lambdaCCLEquilibrium;
    out.lambdaCCLEquilibriumCell = lambdaCCLEquilibriumCell;
    out.lambdaCCLTarget = lambdaCCLTarget;
    out.lambdaCCLTargetCell = lambdaCCLTargetCell;
    out.lambdaSaturationCCL = lambdaSaturationCCL;
    out.lambdaSaturationCCLCell = lambdaSaturationCCLCell;
    out.dlambdaCCL = dlambdaCCL;
    out.tauHyd = tauHyd;
    out.cCLIonomerWaterCoefficient = cCLIonomerWaterCoefficient;
    out.hydrationActivity = hydrationActivity;
    out.hydrationActivityCell = hydrationActivityCell;
    out.hydrationKineticRate = hydrationKineticRate;
    out.hydrationProductionSupplyRate = productionSupplyRate;
    out.hydrationExistingPoreSupplyRate = existingPoreSupplyRate;
    out.freezableProductionRateCCL = freezableProductionRateCCL;
    out.lambdaNonfreezingFraction = lambdaNonfreezingFraction;
    out.nucleationActivation = nucleationActivation;
    out.directFreezeSupply = directFreezeSupply;
    out.mobileWaterActivation = mobileWaterActivation;
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
    out.interfaceActivity.cCLIonomer = thetaCCLIonomer;
    out.interfaceActivity.cCLPore = waterActivity(g.idx_cCL(1));
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
        'finite mass transfer to dry anode and cathode channels';
    out.vaporBoundary.massTransferCoefficient = hVaporCathode;
    out.vaporBoundary.anodeMassTransferCoefficient = hVaporAnode;
    out.vaporBoundary.cathodeMassTransferCoefficient = hVaporCathode;
    out.vaporBoundary.externalConcentration = ...
        cathodeInletVaporConcentration;
    out.vaporBoundary.outwardFlux = JdiffC(end);
    out.vaporBoundary.cathodeOutwardFlux = JdiffC(end);
    out.vaporBoundary.anodeOutwardFlux = -JdiffA(1);
    out.currentInfo = currentInfo;
    out.Sw = Sw;
    out.dmiFull = dmiFull;
    out.Rfreeze = RfreezeFull;
    out.RfreezeLiquidRaw = RfreezeLiquidFull;
    out.RfreezeDirectRaw = RfreezeDirectFull;
    out.Rmelt = RmeltFull;
    out.liquidSaturation = liquidSaturation;
    out.liquidBulkVolumeFraction = liquidBulkVolumeFraction;
    out.iceBulkVolumeFraction = iceBulkVolumeFraction;
    out.kFreeze = kFreeze;
    out.kMelt = kMelt;
    out.phaseTransitionWidthK = phaseTransitionWidthK;
    out.freezeWeight = freezeWeightFull;
    out.directFreezeFraction = directFreezeFraction;
    out.freezeNucleationLambdaFraction = ...
        freezeNucleationLambdaFraction;
    out.rhoIce = rhoIce;
    out.freezingPoint = freezingPoint;
    out.iceProjectionApplied = iceProjectionApplied;
    out.poreProjectionApplied = poreProjectionApplied;
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
% 局部函数：Nafion低温最大非冻结含水量
% ========================================================================

function lambdaSat = nonfreezing_lambda_saturation(T,lambdaScale,Tf)

    T = T(:);
    lambdaSat = lambdaScale*ones(size(T));
    veryCold = T < 223.15;
    subzero = T >= 223.15 & T < Tf;
    lambdaSat(veryCold) = 4.837;
    lambdaSat(subzero) = 1./(-1.304+0.01479.*T(subzero)- ...
        3.594e-5.*T(subzero).^2);
    lambdaSat = min(max(lambdaSat,1e-6),lambdaScale);

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
    elseif strcmp(leftType,'mass_transfer')
        % 左边界沿+x方向为正，向左流出区域时通量为负。
        h = leftValue(1);
        phiExternal = leftValue(2);
        Jface(1) = h*(phiExternal-phi(1));
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
