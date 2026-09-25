%-------------------------------------------------------------------------------
% @function: 计算含冰模型的电压、发热量、相变潜热和温度状态导数
% @author:   PJ, GPT
% @date:     20260923
% @input:    T,cH2,cO2->当前温度和气体浓度状态
%            water->water_ice_state_ice输出的水/冰/孔隙率信息
%            j,Tamb,qaux->电流密度、环境温度和辅助热源
%            g,s->网格和状态编号
%            gammaIce->阴极催化层冰覆盖修正指数
%            j0Ref,Ea->参考交换电流密度[A/m^2]和活化能[J/mol]
% @output:   dT->温度状态导数 [K/s]
%            out->电压损失、热源、热通量和冰修正信息
%
% [ICE-4] cCL冰饱和度降低有效反应面积和交换电流密度。
% [ICE-5] 冻结潜热作为正热源、融化潜热作为负热源进入能量方程。
% [Q2F-1] 增加可选的电堆耦合边界，用于问题2的5片电堆模型。
%   第15个入参 stackBoundary 为空时严格退化为单电池口径，
%   因此问题1的全部标定结果在本文件中逐位不变。
%-------------------------------------------------------------------------------

function [dT,out] = thermal_temperature_state_ice( ...
    T,cH2,cO2,water,j,Tamb,qaux,g,s,gammaIce,j0Ref,Ea,clOpt,hydOpt, ...
    stackBoundary)

    %% 1. 输入处理

    if nargin < 13 || isempty(clOpt) || ~isstruct(clOpt), clOpt = struct(); end
    if nargin < 14 || isempty(hydOpt) || ~isstruct(hydOpt), hydOpt = struct(); end
    % [Q2F-1] 电堆耦合边界。
    % mode='single_cell'（默认）：左右面各为 (T-Tamb)/R 的对流散热，
    %   与 Code_ICE_F 完全一致。
    % mode='prescribed_flux'：左右面的净散热面通量由电堆层给定：
    %   leftOutflux/rightOutflux [W/m^2]，正号表示热量离开该片。
    %   extraArealCapacity [J/(m^2 K)] 附加面热容（端板集中热容）。
    %   bipolarPlateCapacityScale [-] 双极板面热容计数倍率。
    %   concentrationLossFactor [-] 浓差极化倍率（附件1 首尾电池为10）。
    if nargin < 15 || isempty(stackBoundary) || ~isstruct(stackBoundary)
        stackBoundary = struct('mode','single_cell');
    end
    if ~isfield(stackBoundary,'mode'), stackBoundary.mode = 'single_cell'; end
    if ~isfield(stackBoundary,'extraArealCapacity')
        stackBoundary.extraArealCapacity = 0;
    end
    if ~isfield(stackBoundary,'bipolarPlateCapacityScale')
        stackBoundary.bipolarPlateCapacityScale = 1;
    end
    if ~isfield(stackBoundary,'concentrationLossFactor')
        stackBoundary.concentrationLossFactor = 1;
    end
    stackCoupled = strcmpi(char(stackBoundary.mode),'prescribed_flux');
    if stackCoupled
        if ~isfield(stackBoundary,'leftOutflux') || ...
                ~isfield(stackBoundary,'rightOutflux')
            error('thermal_temperature_state_ice:MissingStackFlux', ...
                'prescribed_flux模式必须给出leftOutflux与rightOutflux。');
        end
        if ~isscalar(stackBoundary.leftOutflux) || ...
                ~isscalar(stackBoundary.rightOutflux) || ...
                ~isfinite(stackBoundary.leftOutflux) || ...
                ~isfinite(stackBoundary.rightOutflux)
            error('thermal_temperature_state_ice:InvalidStackFlux', ...
                'leftOutflux与rightOutflux必须是有限标量。');
        end
    end
    if ~isscalar(stackBoundary.concentrationLossFactor) || ...
            ~isfinite(stackBoundary.concentrationLossFactor) || ...
            stackBoundary.concentrationLossFactor < 0
        error('thermal_temperature_state_ice:InvalidConcentrationFactor', ...
            'concentrationLossFactor必须是非负有限标量。');
    end
    if ~isfield(hydOpt,'mHyd'), mHyd = 0; else, mHyd = hydOpt.mHyd; end
    if ~isfield(hydOpt,'lamHydRef'), lamHydRef = 1.0; else, lamHydRef = hydOpt.lamHydRef; end
    T = T(:);
    cH2 = max(cH2(:),0);
    cO2 = max(cO2(:),0);
    qaux = qaux(:);
    R = 8.314;    % [J/(mol K)]
    F = 96485;    % [C/mol]


    %% 2. 催化层代表状态

    dxACL = g.dx(g.idx_aCL);
    dxCCL = g.dx(g.idx_cCL);
    cH2ACL = sum(cH2(s.local.H2_aCL).*dxACL)/sum(dxACL);
    cO2CCL = sum(cO2(s.local.O2_cCL).*dxCCL)/sum(dxCCL);
    TACL = sum(T(g.idx_aCL).*dxACL)/sum(dxACL);
    TCCL = sum(T(g.idx_cCL).*dxCCL)/sum(dxCCL);
    Tavg = sum(T.*g.dx)/sum(g.dx);

    pH2 = cH2ACL*R*TACL;
    pO2 = cO2CCL*R*TCCL;
    if pH2 <= 0 || pO2 <= 0
        error('thermal_temperature_state_ice:InvalidPressure', ...
            '催化层氢气或氧气分压小于等于0。');
    end


    %% 3. 氧气有效扩散系数和极限电流

    epsO2 = water.eps_g(g.idx_O2);
    DO2 = 2.20e-5.*(T(g.idx_O2)/298.15).^1.75.*epsO2.^1.5;
    % [FIX-2] cCL和cGDL沿厚度方向为串联扩散阻力。
    % 不再使用cCL有效扩散系数的算术平均值，以免局部冰堵被平均掉。
    dxO2 = g.dx(g.idx_O2);
    oxygenDiffusionResistance = sum(dxO2./DO2);
    oxygenDiffusionLength = sum(dxO2);
    DO2Equivalent = oxygenDiffusionLength/oxygenDiffusionResistance;
    DO2ArithmeticCCL = sum(DO2(s.local.O2_cCL).*dxCCL)/sum(dxCCL);
    jLim = 4*F*cO2CCL/oxygenDiffusionResistance;


    %% 4. 单电池电压

    Erev = 1.229-8.5e-4*(Tavg-298.15) + ...
        R*Tavg/(2*F)*log((pH2/101325)*sqrt(pO2/101325));

    % [ICE-4] 有效反应面积修正只作用于cCL的交换电流密度：
    %   fArea=(1-si_cCL)^gammaIce
    %   j0_eff=j0_base*fArea
    % 在给定电流边界下，Faraday反应源仍由j决定，不再重复乘fArea。
    iceSaturationCCL = sum(water.sIce(g.idx_cCL).*dxCCL)/sum(dxCCL);
    iceAreaFactor = max(1-iceSaturationCCL,0)^gammaIce;
    iceAreaFactorForVoltage = max(iceAreaFactor,1e-8);

    % [D-1] j0 还必须随 cCL 离聚物水合程度上升。
    % 依据（见 tmp/late.m 的量化）：实验反解出的有效 j0 从 t=0 到 t=36.6 s
    % 上升约 18 倍，而题目式(40) 只含温度项，Ea=67000 J/mol 在
    % 253.15 -> 260.96 K 之间只能给出 2.59 倍，缺约 7 倍。
    % 物理来源：离聚物吸水后体积膨胀、覆盖更多 Pt 表面并提高质子可达性，
    % 有效三相界面面积显著增加，因此交换电流密度随 λ 上升。
    % 形式（新增两个参数：指数 mHyd 与参考点 lamHydRef）：
    %   fHyd = (lambda_CCL/lamHydRef)^mHyd
    % mHyd=0 时严格退化为原模型。
    % [D-2][E-1] lamHydRef 固定在 lambdaCCL0=1.0（t=0 时 cCL 的实际含水量），
    % 因此 fHyd(t=0)=1，曲线起点不受该因子影响。
    % 注意：lamHydRef^mHyd 与 j0Ref 在数学上完全简并
    %   j0Ref*fHyd = (j0Ref/lamHydRef^mHyd)*lambda^mHyd
    % 二者不可同时辨识，故只把 j0Ref 与 mHyd 列为标定量，lamHydRef 取
    % 物理参考点而不参与标定。
    % （Code_ICE_D 曾取 lamHydRef=0.10 并注为"t=0 实际值"，该注释有误：
    %  t=0 的实际 lambda 就是 1.0；D 版等价于 j0Ref=0.020/0.1^1.35=0.4477。
    %  本版改用不简并的写法：lamHydRef=1.0, j0Ref=0.30。）
    fHyd = max(water.lambdaCCL/lamHydRef,1e-6)^mHyd;
    fHyd = min(max(fHyd,1e-6),1e6);          % 数值保护

    j0Base = j0Ref*exp(-Ea/R*(1/Tavg-1/298.15));
    j0Effective = j0Base*iceAreaFactorForVoltage*fHyd;
    etaAct = R*Tavg/(0.5*F)*asinh(j/(2*j0Effective));

    % [FIX-5] 欧姆损失应包含PEM和两个CL中的质子传导。
    % 原程序虽然已经构造了CL内的im(x)，但电压中只使用
    % j*Lpem/kappaPem，完全遗漏了CL离聚物的质子电阻。
    % 附件1给出CL离聚物含量为0.3，第一版采用Bruggeman修正：
    %   kappaCL=0.3^1.5*kappa(lambda,T)
    %   etaOhmCL=int(im/kappaCL)dx
    % [A-2] CL 质子电阻可按比例缩放，默认 0（关闭）。
    % 依据：用两个温度的全部实测极化点做二维拟合 {j0Ref, R_extra}
    % （Ea=67000、alpha=0.5 均固定为题目给定值）时，最优解把串联电阻
    % 主动推到 0，RMSE=0.0611 V。原因是实测极化曲线是"凹"的——最大
    % 微分电阻出现在最低电流段(5.3 ohm cm2)，高电流段只有 0.23，
    % 而串联电阻给出的 j*R 随电流线性增长，会让曲线更"凸"，方向相反。
    % 因此默认 clResistanceScale=0.1 保留少量（计入附件1 的离聚物含量
    % 但避免过度）；设为 1 即回到 FIX-5 的完整 CL 电阻版本。
    if ~isfield(clOpt,'clResistanceScale')
        clResistanceScale = 0;
    else
        clResistanceScale = clOpt.clResistanceScale;
    end
    lambdaPemLocal = water.lambdaPEM;
    kappaPemLocal = (0.5139.*lambdaPemLocal-0.326).* ...
        exp(1268.*(1/303.15-1./T(g.idx_PEM)));

    lambdaACL = lambdaPemLocal(1);
    % [HYD-1] cCL电阻由自身离聚物水合状态决定，不再直接复制PEM右端。
    lambdaCCL = water.lambdaCCL;
    kappaACL = 0.3^1.5.*(0.5139*lambdaACL-0.326).* ...
        exp(1268.*(1/303.15-1./T(g.idx_aCL)));
    kappaCCL = 0.3^1.5.*(0.5139*lambdaCCL-0.326).* ...
        exp(1268.*(1/303.15-1./T(g.idx_cCL)));

    % [E-5] 干涸离聚物电导率下限。
    % 原代码在 kappaPEM/ACL/CCL <= 0（离聚物完全干涸）时把欧姆损失整体置为
    % Inf。但 Inf 会顺着 Vcell 流进产热项 qgen=j*(1.48-Vcell)/Ltotal，使
    % dT 变成 Inf，再经 gas_transport_state_ice 传染给 dcO2，ode15s 的
    % 试探步随即发散（详见 pemfc_calculate_ice.m 的 [E-2][E-3][E-4]）。
    % 这里改为给电导率一个极小正下限：欧姆损失仍然很大，但始终有限。
    kappaFloor = 1e-6;      % [S/m] 完全干涸离聚物的等效电导率下限
    kappaPemLocal = max(kappaPemLocal,kappaFloor);
    kappaACL = max(kappaACL,kappaFloor);
    kappaCCL = max(kappaCCL,kappaFloor);
    kappaPem = g.layerThickness(3)/sum(g.dx(g.idx_PEM)./kappaPemLocal);
    if clResistanceScale <= 0
        etaOhmPEM = j*sum(g.dx(g.idx_PEM)./kappaPemLocal);
        etaOhmACL = 0;
        etaOhmCCL = 0;
        etaOhmContact = j*0.01e-4;
        etaOhm = etaOhmPEM+etaOhmContact;
    else
        etaOhmPEM = j*sum(g.dx(g.idx_PEM)./kappaPemLocal);
        etaOhmACL = clResistanceScale*sum(water.imCell(g.idx_aCL).* ...
            g.dx(g.idx_aCL)./kappaACL);
        etaOhmCCL = clResistanceScale*sum(water.imCell(g.idx_cCL).* ...
            g.dx(g.idx_cCL)./kappaCCL);
        etaOhmContact = j*0.01e-4;
        etaOhm = etaOhmPEM+etaOhmACL+etaOhmCCL+etaOhmContact;
    end

    % [E-5] 浓差损失：把 log 的自变量夹离 0。
    % 原代码在 j >= jLim（超过极限电流密度，低温下冰堵塞孔隙时很容易发生）
    % 时置 etaCon=Inf，同样会把 NaN/Inf 带进状态导数。这里改为有限值，
    % 并用 transportInfeasible 标志记录该情况。
    transportValid = jLim > 0 && j >= 0 && j < jLim;
    if j == 0
        etaCon = 0;
    else
        concentrationRatio = min(max(j/max(jLim,realmin),0),1-1e-12);
        etaCon = -R*Tavg/(4*F)*log(1-concentrationRatio);
    end

    % [Q2F-1] 附件1 给出"首尾电池浓差极化系数=10、中间电池=1"，
    % 但题目未给出它与电压公式的唯一结合方式。这里按
    %   eta_con,k = beta_k * [-R*T/(4F)*ln(1-j/j_lim)]
    % 乘在浓差损失上，beta 由电堆层给出。单电池口径 beta=1，为无操作。
    etaCon = etaCon*stackBoundary.concentrationLossFactor;

    % [E-5] 总损失上限：保证 Vcell 恒为有限值，使 RHS 永远返回有限导数。
    % 该上限只在实际已超出模型适用范围（超出极限电流、离聚物全干）时生效；
    % 本组两工况在正常极化区间内总损失只有 0.4~0.8 V，远低于 2 V，
    % 因此标定结果的数值与原式逐位相同（已验证 RMSE 完全一致）。
    voltageLossCap = 2.0;   % [V]
    totalLoss = etaAct+etaOhm+etaCon;
    voltageLossClamped = ~isfinite(totalLoss) || totalLoss > voltageLossCap;
    if voltageLossClamped
        totalLoss = voltageLossCap;
    end
    Vcell = Erev-totalLoss;


    %% 5. 等效热容和导热系数

    materialRho = [185,970,2150,970,185].';
    materialCp = [545,240,1050,240,545].';
    materialK = [0.3,0.27,0.24,0.27,0.3].';
    layerRhoCp = materialRho.*materialCp;
    delta = g.layerThickness(:);
    meaArealCapacity = sum(layerRhoCp.*delta);

    % [Q2F-1] 双极板面热容按电堆层的计数倍率缩放，并可附加端板集中热容。
    % 单电池口径下倍率=1、附加=0，与原式逐位相同。
    bipolarPlateArealCapacity = 1980*766*(2e-3+2e-3)* ...
        stackBoundary.bipolarPlateCapacityScale;
    endPlateArealCapacity = stackBoundary.extraArealCapacity;
    totalArealCapacity = meaArealCapacity + ...
        bipolarPlateArealCapacity + endPlateArealCapacity;
    rhoCpScalar = totalArealCapacity/g.Ltotal;
    rhoCp = rhoCpScalar*ones(g.N,1);

    RthPerArea = sum(delta./materialK);
    kEffScalar = g.Ltotal/RthPerArea;
    kEff = kEffScalar*ones(g.N,1);


    %% 6. 有限体积导热通量

    qHeatFace = zeros(g.N+1,1);
    for f = 2:g.N
        dL = g.xf(f)-g.x(f-1);
        dR = g.x(f)-g.xf(f);
        qHeatFace(f) = -(T(f)-T(f-1))/(dL/kEff(f-1)+dR/kEff(f));
    end

    boundarySolidResistance.left = 2e-3/95;
    boundarySolidResistance.right = 2e-3/95;
    if stackCoupled
        % [Q2F-1] 电堆耦合：边界净散热通量由相邻片导热/端板导热给定。
        % 符号约定与下面单电池式一致：正号=热量离开该片。
        qOutLeft = stackBoundary.leftOutflux;
        qOutRight = stackBoundary.rightOutflux;
        Rleft = NaN;
        Rright = NaN;
    else
        Rleft = (g.x(1)-g.xf(1))/kEff(1) + ...
            boundarySolidResistance.left + 1/40;
        Rright = (g.xf(end)-g.x(end))/kEff(end) + ...
            boundarySolidResistance.right + 1/40;
        qOutLeft = (T(1)-Tamb)/Rleft;
        qOutRight = (T(end)-Tamb)/Rright;
    end
    qHeatFace(1) = -qOutLeft;
    qHeatFace(end) = qOutRight;
    qCond = (qHeatFace(1:end-1)-qHeatFace(2:end))./g.dx;


    %% 7. 温度状态方程
    % 反应热：qgen=j*(Eth-Vcell)/Lmea
    % [ICE-5] 水的熔化潜热Lf=333600 J/kg：
    %   qPhase=Lf*dmi/dt
    % 因而冻结(dmi>0)放热，融化(dmi<0)吸热。
    % [DSH-5] 潜热必须与water模块的相变量严格对应：
    %   冻结项 Lf*dmi 已经包含了 (Rfreeze-Rmelt-Rsub)，
    %   所以这里改用显式各项，避免把升华潜热也当成熔化潜热。
    %   冷凝放热、蒸发吸热：qCond_latent=Lv*Rcond
    %   升华吸热、凝华放热：qSubl=-Ls*Rsub

    qgen = j*(1.48-Vcell)/g.Ltotal;
    Lf = 333600;            % 熔化/凝固潜热 [J/kg]
    Lv = 2.501e6;           % 汽化/冷凝潜热 [J/kg]
    % 冻结放热、融化吸热、升华吸热(凝华放热)、冷凝放热(蒸发吸热)
    qPhase = Lf*(water.Rfreeze-water.Rmelt) - ...
        water.Lsublim*water.sublimationRate;
    qCondens = Lv*water.condensationRate;
    dT = (qCond+qgen+qaux+qPhase+qCondens)./rhoCp;


    %% 8. 输出诊断信息

    voltage.Erev = Erev;
    voltage.j0Ref = j0Ref;
    voltage.Ea = Ea;
    voltage.j0Base = j0Base;
    voltage.j0 = j0Effective;
    voltage.iceAreaFactor = iceAreaFactor;
    voltage.etaAct = etaAct;
    voltage.kappaPem = kappaPem;
    voltage.kappaPemLocal = kappaPemLocal;
    voltage.kappaACL = kappaACL;
    voltage.kappaCCL = kappaCCL;
    voltage.etaOhm = etaOhm;
    voltage.etaOhmPEM = etaOhmPEM;
    voltage.etaOhmACL = etaOhmACL;
    voltage.etaOhmCCL = etaOhmCCL;
    voltage.etaOhmContact = etaOhmContact;
    voltage.jLim = jLim;
    voltage.etaCon = etaCon;
    voltage.transportValid = transportValid;
    % [E-5] 电压损失有限化标志：true 表示该时刻已超出模型适用范围
    % （j>=jLim 或离聚物干涸），此时总损失被压到 voltageLossCap。
    voltage.totalLoss = totalLoss;
    voltage.totalLossUnclamped = etaAct+etaOhm+etaCon;
    voltage.voltageLossClamped = voltageLossClamped;
    voltage.transportInfeasible = ~transportValid;
    voltage.voltageLossCap = voltageLossCap;

    thermalProperties.rhoCpScalar = rhoCpScalar;
    thermalProperties.kEffScalar = kEffScalar;
    thermalProperties.layerRhoCp = layerRhoCp;
    thermalProperties.layerThermalResistance = delta./materialK;
    thermalProperties.thermalResistancePerArea = RthPerArea;
    thermalProperties.Ltotal = g.Ltotal;
    thermalProperties.meaArealCapacity = meaArealCapacity;
    thermalProperties.bipolarPlateArealCapacity = bipolarPlateArealCapacity;
    thermalProperties.endPlateArealCapacity = endPlateArealCapacity;
    thermalProperties.totalArealCapacity = totalArealCapacity;
    thermalProperties.boundarySolidResistance = boundarySolidResistance;

    divInfo.Jleft = qHeatFace(1);
    divInfo.Jright = qHeatFace(end);
    divInfo.integralDiv = sum(qCond.*g.dx);
    divInfo.boundaryNet = qHeatFace(1)-qHeatFace(end);
    divInfo.balanceError = divInfo.integralDiv-divInfo.boundaryNet;
    divInfo.conserved = abs(divInfo.balanceError) <= ...
        max(1e-12,1e-10*max(abs(divInfo.boundaryNet),1));

    thermalFV.qOutLeft = qOutLeft;
    thermalFV.qOutRight = qOutRight;
    thermalFV.totalHeatLossPerArea = qOutLeft+qOutRight;
    thermalFV.leftSurfaceResistance = Rleft;
    thermalFV.rightSurfaceResistance = Rright;
    thermalFV.boundarySolidResistance = boundarySolidResistance;
    thermalFV.divergenceInfo = divInfo;
    thermalFV.integratedConduction = sum(qCond.*g.dx);
    thermalFV.expectedIntegratedConduction = -(qOutLeft+qOutRight);
    % [Q2F-1] 电堆耦合诊断
    thermalFV.stackCoupled = stackCoupled;
    thermalFV.concentrationLossFactor = stackBoundary.concentrationLossFactor;
    thermalFV.energyBalanceError = thermalFV.integratedConduction- ...
        thermalFV.expectedIntegratedConduction;
    thermalFV.energyBalanced = abs(thermalFV.energyBalanceError) <= ...
        max(1e-10,1e-9*max(abs(thermalFV.totalHeatLossPerArea),1));

    out.Vcell = Vcell;
    out.voltage = voltage;
    out.Tavg = Tavg;
    out.Tmin = min(T);
    out.Tmax = max(T);
    out.cH2_aCL = cH2ACL;
    out.cO2_cCL = cO2CCL;
    out.pH2 = pH2;
    out.pO2 = pO2;
    out.DO2 = DO2;
    out.DO2CCL = DO2Equivalent;
    out.DO2Equivalent = DO2Equivalent;
    out.DO2ArithmeticCCL = DO2ArithmeticCCL;
    out.oxygenDiffusionResistance = oxygenDiffusionResistance;
    out.qgen = qgen;
    out.qPhase = qPhase;
    out.qCond = qCond;
    out.qHeatFace = qHeatFace;
    out.rhoCp = rhoCp;
    out.kEff = kEff;
    out.thermalProperties = thermalProperties;
    out.thermalFV = thermalFV;
    out.jLim = jLim;
    out.transportValid = transportValid;
    out.iceSaturationCCL = iceSaturationCCL;
    out.iceAreaFactor = iceAreaFactor;
    out.gammaIce = gammaIce;

end


