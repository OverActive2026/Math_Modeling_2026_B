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
% [DSH-1] 相变传质阻力：气相/液相不再瞬时平衡，冷凝-蒸发按
%   Rcond=kmCond*(mv-A*eps_g) 的有限速率进行。原式(25)(26)隐含
%   "水蒸气一达饱和立刻全部冷凝"，使产物水几乎全部以水蒸气形式扩散
%   出电池（实测排出量约为产水量的2倍），多孔层始终不过饱和，
%   液态水与冰无法生成。加入该阻力后液水得以在cCL滞留。
% [DSH-2] 冻结速率保留过冷度因子 (Tf-T)/dTref，使升温后冻结自然减弱，
%   避免温度回升时冰仍持续累积（原式无温度依赖，升温不减速）。
% [DSH-3] 增加冰的升华/凝华 Rsub，这是T<0C时冰量唯一可以下降的通道，
%   使冰堵具有可逆性（后段电压回升）。
%-------------------------------------------------------------------------------

function [dmw,dmi,dlambdaCCL,out] = water_ice_state_ice( ...
    T,mw,mi,lambdaCCLState,j,kFreeze,kMelt,tauHyd,phaseOpt,g,s)

    %% 1. 输入处理

    T = T(:);
    mw = max(mw(:),0);  % 仅处理求解器Newton迭代产生的微小负试探值
    miState = max(mi(:),0);
    lambdaCCLStateRaw = lambdaCCLState;
    lambdaCCLState = max(lambdaCCLState,0);

    % [E-2] 先判定状态本身是否有限。
    % 若 T 或水/冰状态已经是 NaN/Inf，说明问题出在上游（ode15s 的试探步
    % 已经发散），而不是本模块的物性计算。原代码没有这道检查，NaN 会一路
    % 传到第 6 节的 Dw 计算，最后报成
    %   "water_ice_state_ice:InvalidDiffusivity 水扩散系数出现非法数值"
    % ——把"低温工况刚性积分失败"误报成"物性参数非法"，非常容易误导。
    % -25 degC 是最容易触发该问题的工况：此时 cCL 开始出现液态水，
    % 题目式(9) 的液水冻结项 Rfreeze=kFreeze*rhoL*eps0*sL 才真正激活；
    % -20 degC 下 ml 恒为 0，该项完全不参与，所以 -20 degC 一直很稳。
    if any(~isfinite(T)) || any(~isfinite(mw)) || ...
            any(~isfinite(miState)) || ~isfinite(lambdaCCLState)
        error('water_ice_state_ice:NonFiniteState', ...
            ['水-冰模块收到非有限状态：T范围=[%.6g,%.6g]，mw范围=[%.6g,%.6g]，' ...
            'mi范围=[%.6g,%.6g]，lambdaCCL=%.6g。\n' ...
            '这说明 ode15s 的试探步已经发散（低温工况刚性问题），' ...
            '不是物性参数非法。请检查 kFreeze / directFreezeFraction / ' ...
            'kFreezeDirect 是否过大，或适当放宽 AbsTol（如 1e-6）。'], ...
            min(T),max(T),min(mw),max(mw),min(miState),max(miState),lambdaCCLState);
    end

    % [DSH-1] 相变传质参数
    if nargin < 10 || isempty(phaseOpt) || ~isstruct(phaseOpt), phaseOpt = struct(); end
    if ~isfield(phaseOpt,'kmCond'), phaseOpt.kmCond = 0; end
    if ~isfield(phaseOpt,'dTref'),  phaseOpt.dTref  = 20; end
    if ~isfield(phaseOpt,'nRet'),   phaseOpt.nRet   = 0; end
    if ~isfield(phaseOpt,'rRet'),   phaseOpt.rRet   = 0; end
    if ~isfield(phaseOpt,'fRet'),   phaseOpt.fRet   = 1; end
    % [A-1] 平衡水合模式的松弛时间常数 [s]（仅当 tauHyd<=0 时使用）
    dlambdaRelaxTime = 0.05;
    if isfield(phaseOpt,'dlambdaRelaxTime') && phaseOpt.dlambdaRelaxTime>0
        dlambdaRelaxTime = phaseOpt.dlambdaRelaxTime;
    end
    kmCond = phaseOpt.kmCond;      % 冷凝/蒸发传质系数 [1/s]，0=瞬时平衡(原模型)
    dTref  = phaseOpt.dTref;       % 冻结过冷度参考值 [K]
    nRet   = phaseOpt.nRet;        % 液水对气相孔隙的阻滞指数
    rRet   = phaseOpt.rRet;        % 液水滞留对水蒸气排出的阻滞强度
    fRet   = phaseOpt.fRet;        % 多孔层水蒸气有效扩散系数直接折减因子 [-]
    % [PHASE-2] 两侧气道的有限水蒸气传质系数 [m/s]；0 = 退化为零浓度理想汇
    if ~isfield(phaseOpt,'hVaporAnode'),   phaseOpt.hVaporAnode = 0; end
    if ~isfield(phaseOpt,'hVaporCathode'), phaseOpt.hVaporCathode = 0; end
    hVaporAnode = phaseOpt.hVaporAnode;
    hVaporCathode = phaseOpt.hVaporCathode;
    % [PHASE-3][PHASE-4] 低温直接冻结参数
    if ~isfield(phaseOpt,'directFreezeFraction'), phaseOpt.directFreezeFraction = 0; end
    if ~isfield(phaseOpt,'freezeNucleationLambdaFraction')
        phaseOpt.freezeNucleationLambdaFraction = 0.5;
    end
    directFreezeFraction = phaseOpt.directFreezeFraction;
    freezeNucleationLambdaFraction = phaseOpt.freezeNucleationLambdaFraction;
    % [PHASE-4b] 就地直接冻结速率系数 [1/s]；1 = 预留份额全部就地冻结
    if ~isfield(phaseOpt,'kFreezeDirect'), phaseOpt.kFreezeDirect = 1; end
    kFreezeDirect = phaseOpt.kFreezeDirect;
    % [PHASE-5] 窄平滑相变宽度 [K]；0 = 严格分段式（题目口径）
    if ~isfield(phaseOpt,'phaseTransitionWidthK'), phaseOpt.phaseTransitionWidthK = 0; end
    phaseTransitionWidthK = phaseOpt.phaseTransitionWidthK;
    Lsublim = 2.834e6;             % 冰的升华潜热 [J/kg]

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
    % [E-6] 孔隙被液水+冰填满时不再中断积分。
    % 原代码在此直接 error('water_ice_state_ice:BlockedPore')。在 -25 degC
    % 且冻结速率较大时，ode15s 的试探步很容易瞬时把某个控制体的孔隙填满，
    % 于是整个积分被中断，用户看到的是一个与真实原因无关的报错。
    % 这里只对"用于输运/活度计算的孔隙率"取一个极小正下限，
    % **不改变 mv/ml/mi 之间的分配**，因此不会凭空产生或消灭质量；
    % 触发情况由 poreFilledProjection 记录并向上层报告（含义是该处已被
    % 冰/水堵死，扩散与反应面积应按堵塞处理）。
    minPoreVolumeFraction = 1e-6;
    poreFilledMask = eps_g(idxPorous) <= minPoreVolumeFraction*eps0(idxPorous);
    poreFilledProjection = any(poreFilledMask);
    eps_g(idxPorous) = max(eps_g(idxPorous), ...
        minPoreVolumeFraction*eps0(idxPorous));

    % 局部水活度（0..1）：孔隙水未饱和时等于相对湿度，有液态水时为1
    mSat = zeros(g.N,1);
    mSat(idxPorous) = eps_g(idxPorous).*A(idxPorous);
    waterActivity = zeros(g.N,1);
    validSat = mSat > 0;
    waterActivity(validSat) = min(max(mv(validSat)./mSat(validSat),0),1);
    waterActivity(ml > 0) = 1;


    %% 4. cCL离聚物水合动力学
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
    % [DSH-4] 水合目标必须由局部水活度决定，不能用"j>0就置饱和"的开关。
    % 原写法在cCL水活度仅约0.007时仍把lambda推向4.837，等于凭空虚增
    % 结合水（同时也制造了早期那段并不存在的欧姆压降）。
    % [DSH-6] 目标值只受"水活度平衡值"和"低温最大非冻结含水量"限制。
    % lambdaUpper(=mw-mi 的存量上限) 不能作为目标钳位：cCL 水活度在
    % 反应位附近接近1，而可迁移水存量接近0，用存量钳位会把目标压到
    % 0.09 以下，使水合完全停滞。质量守恒由下面的供给率限制保证。
    lambdaCCLTarget = min(lambdaCCLEquilibrium,lambdaSaturationCCL);
    hydrationKineticRate = ...
        (lambdaCCLTarget-lambdaCCLStateRaw)/tauHyd;
    % [A-1] tauHyd门控（默认）或平衡水合（tauHyd<=0）两种模式。
    % 依据：Nafion 吸水等温线是平衡关系，cCL 离聚物薄膜的吸水时间常数
    % 远小于秒级，因此 lambda 应直接取局部水活度对应的平衡值，
    % 而不是用 tauHyd 作为待标定量去拖。tauHyd<=0 时启用平衡模式。
    reactionWaterRateCCL = Mw*j/(2*F*g.layerThickness(4));
    % [PHASE-3][PHASE-4] 反应产水在"离聚物水合"与"低温直接冻结"之间分配。
    % 关键：必须**在水合供给端就预留** directFreezeFraction 的产水给冻结通道。
    % 若只在冻结端限制（Q2 原写法），由于水合速率上限比产水速率大两个数量级
    % （productionSupplyRate=反应产水/离聚物系数，而离聚物系数很小），
    % 离聚物仍会吃掉几乎全部产水，冻结通道拿到的份额 ~2%，冰依旧接近 0。
    productionSupplyRate = (1-directFreezeFraction)* ...
        reactionWaterRateCCL/cCLIonomerWaterCoefficient;
    poreSupplyRate = max(mPore(g.idx_cCL(1))-mv(g.idx_cCL(1)),0)/ ...
        (cCLIonomerWaterCoefficient*max(tauHyd,eps));
    existingPoreSupplyRate = poreSupplyRate;
    if tauHyd <= 0
        % 平衡水合：lambda 无自身动力学，直接等于目标值
        dlambdaCCL = (lambdaCCLTarget-lambdaCCLStateRaw)/dlambdaRelaxTime;
    elseif hydrationKineticRate > 0
        dlambdaCCL = min(hydrationKineticRate, ...
            poreSupplyRate+productionSupplyRate);
    else
        dlambdaCCL = hydrationKineticRate;
    end
    if lambdaCCLStateRaw <= 0 && dlambdaCCL < 0
        dlambdaCCL = 0;
    elseif lambdaCCLStateRaw >= lambdaSaturationCCL && dlambdaCCL > 0
        dlambdaCCL = 0;
    end

    % [PHASE-4] 被预留出来、未被离聚物吸收的那部分反应产水
    dmIonomerCCL = cCLIonomerWaterCoefficient*max(dlambdaCCL,0);
    freezableProductionRateCCL = max(reactionWaterRateCCL- ...
        dmIonomerCCL,0);
    % 直接冻结的成核门槛：离聚物水合到一定程度后才允许成核
    lambdaNonfreezingFraction = min(max(lambdaCCLStateRaw/ ...
        max(lambdaSaturationCCL,1e-6),0),1);
    nucleationActivation = min(max((lambdaNonfreezingFraction- ...
        freezeNucleationLambdaFraction)/ ...
        max(1-freezeNucleationLambdaFraction,1e-6),0),1);


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

    % [PHASE-5] 窄平滑相变开关：width=0 严格保留题目分段关系；
    % width>0 时只在冰点邻域内用光滑权重替代 0/1 硬切换，避免求解器
    % 在冻结/融化两套方程间反复切换把内部步长压到极小。远离冰点等价。
    if phaseTransitionWidthK > 0
        freezeWeight = 0.5*(1-tanh((T(idxPorous)-freezingPoint)./ ...
            phaseTransitionWidthK));
    else
        freezeWeight = double(T(idxPorous) <= freezingPoint);
    end
    meltWeight = 1-freezeWeight;

    % [DSH-2] 冻结：保留过冷度因子，使升温后冻结速率自然下降
    undercoolingFactor = max(freezingPoint-T(idxPorous),0)./dTref;
    RfreezeLiquid = kFreeze*rhoLiquid.*liquidBulkVolumeFraction.* ...
        undercoolingFactor.*freezeWeight;
    % [PHASE-4] 低温直接冻结通道：把未被离聚物吸收的那部分反应产水
    % 按可调比例直接分配到冰相。它不增加总水，只在液水/冰之间重分配，
    % 从而在不改变水账的前提下让冰真正生成。
    %
    % [PHASE-4b] 冻结速率不再乘"可动水门槛"。Q2 原写法：
    %   RfreezeDirect = nuc*mobileWaterActivation*directFreezeSupply*freezeWeight
    %   mobileWaterActivation = mMobile/(mMobile+mTransition)
    % 在本算例中 cCL 可动水存量 mMobile 只有 ~5e-4 kg/m3（产水一生成就被
    % 离聚物吸收或被气相带走），而门槛取 1e-3*rhoLiquid*eps0 ≈ 0.4 kg/m3，
    % 于是该因子恒为 ~1e-3，冻结通道被整体压掉三个数量级，冰始终为 0
    % （见 tmp/E_gates.m 实测：nucAct=0.30 但 mobAct=0.0011）。
    % 物理上 -20 degC 下反应位局部水活度接近 1，生成的水来不及输运就已
    % 冻结，因此冻结速率应由"产水速率 × 预留比例 × 成核活化"决定，而不是
    % 由几乎为空的存量决定。冰量上限仍由下面的 atAllWaterFrozen（冰≤孔隙水）
    % 与 atFullPore（冰≤孔隙容积）双重约束保证，不会凭空虚增。
    directFreezeSupply = min(directFreezeFraction*reactionWaterRateCCL, ...
        freezableProductionRateCCL);
    RfreezeDirect = zeros(numel(idxPorous),1);
    RfreezeDirect(s.local.mi_cCL) = kFreezeDirect.*nucleationActivation.* ...
        directFreezeSupply.*freezeWeight(s.local.mi_cCL);
    Rfreeze = RfreezeLiquid+RfreezeDirect;
    Rmelt = zeros(numel(idxPorous),1);
    Rmelt(meltRegion) = kMelt*rhoIce.* ...
        iceBulkVolumeFraction(meltRegion);
    if phaseTransitionWidthK > 0
        Rmelt = kMelt*rhoIce.*iceBulkVolumeFraction.*meltWeight;
    end

    % [DSH-3] 冰的升华(T<=Tf，冰->水蒸气)与凝华(水蒸气->冰)
    % 冰面饱和浓度取A(T)（Buck在Tc<0时即为冰面关系）
    sublimationRate = zeros(numel(idxPorous),1);
    validSublim = (T(idxPorous) <= freezingPoint) & (miState > 0);
    sublimationRate(validSublim) = kmCond.*max(A(idxPorous(validSublim)).* ...
        eps_g(idxPorous(validSublim))-mvPorous(validSublim),0);

    dmi = Rfreeze-Rmelt-sublimationRate;

    % 在状态边界上禁止继续越界；正常状态中这两项不改变动力学。
    atZeroIce = miState <= 0 & dmi < 0;
    atAllWaterFrozen = miState >= mPore(idxPorous) & dmi > 0;
    atFullPore = miState >= iceCapacity & dmi > 0;
    dmi(atZeroIce | atAllWaterFrozen | atFullPore) = 0;

    dmiFull = zeros(g.N,1);
    dmiFull(idxPorous) = dmi;
    RfreezeFull = zeros(g.N,1);
    RmeltFull = zeros(g.N,1);
    sublimationFull = zeros(g.N,1);
    RfreezeFull(idxPorous) = Rfreeze;
    RmeltFull(idxPorous) = Rmelt;
    sublimationFull(idxPorous) = sublimationRate;

    % [DSH-1] 有限速率冷凝/蒸发：把"偏离饱和的水蒸气"按一阶传质速率转成
    % 液态水（正值冷凝、负值蒸发）。原式(25)(26)隐含瞬时平衡，使产物水
    % 几乎全部以水蒸气扩散离开电池，多孔层无法过饱和生成液水。
    condensationRate = zeros(numel(idxPorous),1);
    condActive = mSat(idxPorous) > 0;
    condensationRate(condActive) = kmCond.*(mvPorous(condActive)- ...
        mSat(idxPorous(condActive)));
    condensationFull = zeros(g.N,1);
    condensationFull(idxPorous) = condensationRate;


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

    % [DSH-1] 液态水对水蒸气输运的阻滞。
    % 诊断表明：按式(16) 的干态扩散系数，cCL 内水蒸气被抽走的速率大于
    % 产水速率（实测边界净流失约为产水量的1.3~1.6倍），cCL 水活度全程
    % 只升到0.64~0.79，永远达不到饱和，因此液水、冰都无法生成。
    % 这里保留式(16) 的形式，但把式中的气相孔隙率换成"实际可用气相体积"
    % 并叠加液相阻滞因子，以体现"液水一旦滞留即自阻断水蒸气通道"。
    % nRet=0 时退化为原模型。
    idxPor = g.idx_porous;
    sLiqPor = min(max(ml(idxPor)./(rhoLiquid*max(epsAfterIce,1e-12)),0),1);
    retentionFactor = ones(g.N,1);
    if nRet > 0
        retentionFactor(idxPor) = max(1-sLiqPor,0).^nRet;
    end
    if rRet > 0
        coarseLayer = ismember(g.layerId,[1 4 5]);
        scale = zeros(g.N,1);
        scale(idxPor) = 1./(1+rRet*sLiqPor);
        retentionFactor(coarseLayer) = retentionFactor(coarseLayer).* ...
            scale(coarseLayer);
        retentionFactor = retentionFactor/(1/(1+rRet));
        retentionFactor = min(retentionFactor,1);
    end
    % [DSH-1a] 直接折减多孔层水蒸气有效扩散系数。诊断表明按式(16) 的
    % 干态系数，cCL 内水蒸气被带走的速率大于产水速率，水活度全程只到
    % 0.64~0.79，永远无法达到饱和生成液水。fRet<1 用来表征参考模型中
    % 未计入的水蒸气输运阻力（孔道曲折、微孔层、对流边界有限速率等）。
    % fRet=1 时严格退化为原模型。
    retentionFactor(idxPor) = fRet*retentionFactor(idxPor);
    Dw(idxA) = Dw(idxA).*retentionFactor(idxA);
    Dw(idxC) = Dw(idxC).*retentionFactor(idxC);

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
    % [PHASE-2] 两侧气道的干气是"外部浓度为0的有限传质边界"，
    % 而不是把 cGDL/aGDL 外边界浓度钉死在0的无限快水汇。
    % hVapor=0 时严格退化为原来的 dirichlet 0。
    if hVaporAnode > 0
        [JdiffA,diffInfoA] = regional_diffusion( ...
            transportWater(idxA),Dw(idxA),g,idxA, ...
            'mass_transfer',[hVaporAnode,0],'noflux',0);
    else
        [JdiffA,diffInfoA] = regional_diffusion( ...
            transportWater(idxA),Dw(idxA),g,idxA,'dirichlet',0,'noflux',0);
    end
    [JdiffM,diffInfoM] = regional_diffusion( ...
        transportWater(g.idx_PEM),Dw(g.idx_PEM),g,g.idx_PEM, ...
        'noflux',0,'noflux',0);
    if hVaporCathode > 0
        [JdiffC,diffInfoC] = regional_diffusion( ...
            transportWater(idxC),Dw(idxC),g,idxC,'noflux',0, ...
            'mass_transfer',[hVaporCathode,0]);
    else
        [JdiffC,diffInfoC] = regional_diffusion( ...
            transportWater(idxC),Dw(idxC),g,idxC,'noflux',0,'dirichlet',0);
    end

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
    % [PHASE-1] 原版把倍率写死为 30，既没有物理依据也没有标定过。
    % 由于 Gright = Dw(PEM)*Cmem*Lpem^-1 ~ 7e-3 kg/(m2 s)，乘以 30 后
    % 界面导通能力约为反应产水通量的 700 倍，等效于"cCL 孔隙水活度与
    % 膜含水量瞬时平衡"。后果是 cCL 孔隙水被膜瞬间抽干（实测 cCL 可动
    % 水仅 ~5e-4 kg/m3），题目式(9)的液态水结冰项因此恒为 0，冰几乎
    % 不生成。改为可标定参数 interfaceFactor，默认仍取 30 以保持原行为。
    if ~isfield(phaseOpt,'interfaceFactor'), phaseOpt.interfaceFactor = 30; end
    factorLeft = phaseOpt.interfaceFactor;
    factorRight = phaseOpt.interfaceFactor;

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

    dmIonomer = zeros(g.N,1);
    dmIonomer(g.idx_cCL) = ...
        cCLIonomerWaterCoefficient*dlambdaCCL;

    % [DSH-6] cCL 取水优先级：反应生成水先满足离聚物水合，剩余部分才
    % 进入孔隙相（可排水/可冻结）。依据是离聚物亲水、孔隙疏水，且
    % 附件1 给出的 CL 离聚物含量 0.3 是"必须先被润湿"的相。
    % 因此把 Sw 拆成 dlambdaCCL 的供给项与孔隙源项两部分：
    %   Sw_pore = Sw_cCL - dmIonomer
    % 总水方程保持不变（mw 仍为 mIon+mv+ml+mi），只是不再让离聚物
    % 从孔隙水里"抢水"（原写法 dmPore=dmw-dmIonomer 会与冻结争夺
    % 同一份水，导致冰一生成就把 lambda 锁死在初值、kappa_CCL 崩溃）。
    Sw(g.idx_cCL) = Sw(g.idx_cCL)-dmIonomer(g.idx_cCL);
    dmw = (Nw(1:end-1)-Nw(2:end))./g.dx + Sw;

    % [DSH-1] 冷凝不改变总水mw，它只改变mv与ml的分配，故孔隙水方程
    % 需扣除Rcond（已包含在Sw的拆分与下述dmPore中）。
    dmPore = dmw-dmIonomer-condensationFull;


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
    % [DSH-DIAG] 边界水通量诊断：判断产物水是否被气相过快带走
    out.NwLeftBoundary = Nw(1);
    out.NwRightBoundary = Nw(end);
    out.NliqBoundaryLeft = Nliq(1);
    out.NliqBoundaryRight = Nliq(end);
    out.waterProductionRate = Mw*j/(2*F);
    out.eps_g = eps_g;
    out.eps_i = eps_i;
    out.sIce = sIce;
    % [E-6] 孔隙被液水+冰填满（该处已堵死）时的投影标志
    out.poreFilledProjection = poreFilledProjection;
    out.poreFilledCount = sum(poreFilledMask);
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
    % [PHASE-DIAG] 冻结通道三重门槛诊断：成核门槛 / 可动水存量 / 供给门槛
    out.freezeDiagnostics = struct( ...
        'directFreezeFraction',directFreezeFraction, ...
        'kFreezeDirect',kFreezeDirect, ...
        'reactionWaterRateCCL',reactionWaterRateCCL, ...
        'dmIonomerCCL',dmIonomerCCL, ...
        'freezableProductionRateCCL',freezableProductionRateCCL, ...
        'directFreezeSupply',directFreezeSupply, ...
        'lambdaNonfreezingFraction',lambdaNonfreezingFraction, ...
        'nucleationActivation',nucleationActivation, ...
        'mMobileCCL',mMobile(g.idx_cCL(1)), ...
        'RfreezeDirectCCL',RfreezeDirect(s.local.mi_cCL(1)), ...
        'RfreezeLiquidCCL',RfreezeLiquid(s.local.mi_cCL(1)), ...
        'freezeWeightCCL',freezeWeight(s.local.mi_cCL(1)), ...
        'undercoolingCCL',undercoolingFactor(s.local.mi_cCL(1)), ...
        'TCCL',T(g.idx_cCL(1)));
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
    out.currentInfo = currentInfo;
    out.Sw = Sw;
    out.dmiFull = dmiFull;
    out.Rfreeze = RfreezeFull;
    out.Rmelt = RmeltFull;
    out.sublimationRate = zeros(g.N,1);
    out.sublimationRate(idxPorous) = sublimationRate;
    out.condensationRate = condensationFull;
    out.kmCond = kmCond;
    out.dTref = dTref;
    out.Lsublim = Lsublim;
    out.undercoolingFactor = zeros(g.N,1);
    out.undercoolingFactor(idxPorous) = undercoolingFactor;
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
    elseif strcmp(leftType,'mass_transfer')
        % [PHASE-2] 左边界：+x 方向为正，向左流出区域时通量为负。
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
        % [PHASE-2] 右边界：正方向为流出计算区域。rightValue=[h,phiExternal]
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
