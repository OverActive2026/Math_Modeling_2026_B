%-------------------------------------------------------------------------------
% @function: 计算含冰模型的电压、发热量、相变潜热和温度状态导数
% @author:   PJ, GPT
% @date:     20260923
% @input:    T,cH2,cO2->当前温度和气体浓度状态
%            water->water_ice_state_ice输出的水/冰/孔隙率信息
%            j,Tamb,qaux->电流密度、环境温度和辅助热源
%            g,s->网格和状态编号
%            gammaIce->阴极催化层冰覆盖修正指数
%            j0Ref,Ea->充分水合时的参考交换电流密度[A/m^2]和活化能[J/mol]
%            fHydDry,lambdaHydOn,lambdaHydWet,nHyd->cCL水合反应面积参数
%            ionomerBruggemanExponent->CL离聚物有效质子电导指数
%            heatBoundary->问题2电堆热边界（可选）：
%                .mode='prescribed_flux'
%                .leftFlux/.rightFlux 按+x方向为正 [W/m^2]
%                .extraArealCapacity 附加面热容 [J/(m^2 K)]
%                .bipolarPlateCapacityScale 双极板热容倍率 [-]
%                .concentrationLossFactor 电堆位置浓差损失系数 [-]
% @output:   dT->温度状态导数 [K/s]
%            out->电压损失、热源、热通量和冰修正信息
%
% [ICE-4] cCL冰饱和度降低有效反应面积和交换电流密度。
% [ICE-5] 冻结潜热作为正热源、融化潜热作为负热源进入能量方程。
% [Q2-1] 五片电堆由外层程序给定左右热通量，单电池内部
%   仍使用原有有限体积导热方程。不传入heatBoundary时完全保持问题1边界。
%-------------------------------------------------------------------------------

function [dT,out] = thermal_temperature_state_ice( ...
    T,cH2,cO2,water,j,Tamb,qaux,g,s,gammaIce,j0Ref,Ea, ...
    fHydDry,lambdaHydOn,lambdaHydWet,nHyd, ...
    ionomerBruggemanExponent,heatBoundary)

    %% 1. 输入处理

    T = T(:);
    cH2 = max(cH2(:),0);
    cO2 = max(cO2(:),0);
    qaux = qaux(:);
    R = 8.314;    % [J/(mol K)]
    F = 96485;    % [C/mol]

    if nargin < 18 || isempty(heatBoundary)
        heatBoundary = struct('mode','standalone', ...
            'leftFlux',NaN,'rightFlux',NaN,'extraArealCapacity',0);
    end
    if ~isscalar(ionomerBruggemanExponent) || ...
            ~isfinite(ionomerBruggemanExponent) || ...
            ionomerBruggemanExponent <= 0
        error('thermal_temperature_state_ice:InvalidIonomerExponent', ...
            'ionomerBruggemanExponent必须是正有限标量。');
    end
    if ~isfield(heatBoundary,'mode'), heatBoundary.mode = 'standalone'; end
    if ~isfield(heatBoundary,'leftFlux'), heatBoundary.leftFlux = NaN; end
    if ~isfield(heatBoundary,'rightFlux'), heatBoundary.rightFlux = NaN; end
    if ~isfield(heatBoundary,'extraArealCapacity')
        heatBoundary.extraArealCapacity = 0;
    end
    if ~isfield(heatBoundary,'bipolarPlateCapacityScale')
        heatBoundary.bipolarPlateCapacityScale = 1;
    end
    if ~isfield(heatBoundary,'concentrationLossFactor')
        heatBoundary.concentrationLossFactor = 1;
    end
    if ~isscalar(heatBoundary.extraArealCapacity) || ...
            ~isfinite(heatBoundary.extraArealCapacity) || ...
            heatBoundary.extraArealCapacity < 0
        error('thermal_temperature_state_ice:InvalidExtraHeatCapacity', ...
            'extraArealCapacity必须是非负有限标量。');
    end
    if ~isscalar(heatBoundary.bipolarPlateCapacityScale) || ...
            ~isfinite(heatBoundary.bipolarPlateCapacityScale) || ...
            heatBoundary.bipolarPlateCapacityScale < 0
        error('thermal_temperature_state_ice:InvalidBipolarCapacityScale', ...
            'bipolarPlateCapacityScale必须是非负有限标量。');
    end
    if ~isscalar(heatBoundary.concentrationLossFactor) || ...
            ~isfinite(heatBoundary.concentrationLossFactor) || ...
            heatBoundary.concentrationLossFactor <= 0
        error('thermal_temperature_state_ice:InvalidConcentrationFactor', ...
            'concentrationLossFactor必须是正有限标量。');
    end


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

    % [ICE-4] 冰覆盖修正只作用于cCL的交换电流密度：
    %   fIce=(1-si_cCL)^gammaIce
    %
    % [HYD-3] 将“离聚物导电能力”和“催化剂的质子可达性”分开。
    % lambdaCCL仍通过kappaCCL决定质子欧姆损失；此处另用有界函数
    % 表示干态时只有部分催化表面接入连续质子通道：
    %   thetaHyd=clip((lambdaCCL-lambdaHydOn)/
    %                 (lambdaHydWet-lambdaHydOn),0,1)
    %   fHyd=fHydDry+(1-fHydDry)*thetaHyd^nHyd
    %   j0_eff=j0_wet(T)*fHyd*fIce
    % 因此启动初期的干CL惩罚被保留，充分水合后fHyd最多恢复到1，
    % 不使用时间开关，也不会无限放大交换电流。
    iceSaturationCCL = sum(water.sIce(g.idx_cCL).*dxCCL)/sum(dxCCL);
    iceAreaFactor = max(1-iceSaturationCCL,0)^gammaIce;
    iceAreaFactorForVoltage = max(iceAreaFactor,1e-8);
    lambdaCCL = water.lambdaCCL;
    % [FIT-HYD-T] lambdaHydWet是常温下充分水合阈值；低温时
    % 离聚物可保持的非冻结水上限本来就会下降。若仍固定用8
    % 作分母，会在-25 degC工况中对质子可达面积重复惩罚。
    % 因此充分水合阈值取标定值与当前非冻结上限的较小者。
    lambdaHydWetEffective = min(lambdaHydWet, ...
        water.lambdaSaturationCCL);
    lambdaHydWetEffective = max(lambdaHydWetEffective, ...
        lambdaHydOn+1e-6);
    hydrationDegree = min(max((lambdaCCL-lambdaHydOn)/ ...
        (lambdaHydWetEffective-lambdaHydOn),0),1);
    hydrationAreaFactor = fHydDry+(1-fHydDry)*hydrationDegree^nHyd;
    hydrationAreaFactorForVoltage = max(hydrationAreaFactor,1e-8);

    % [ACT-2] j0Ref改为充分水合状态的参考值。干态零时刻的
    % 有效值为j0Ref*fHydDry，避免把初始干燥性永久吸收进本征j0。
    j0Base = j0Ref*exp(-Ea/R*(1/Tavg-1/298.15));
    j0Effective = j0Base*hydrationAreaFactorForVoltage* ...
        iceAreaFactorForVoltage;
    etaAct = R*Tavg/(0.5*F)*asinh(j/(2*j0Effective));

    % [FIX-5] 欧姆损失应包含PEM和两个CL中的质子传导。
    % 原程序虽然已经构造了CL内的im(x)，但电压中只使用
    % j*Lpem/kappaPem，完全遗漏了CL离聚物的质子电阻。
    % 附件1给出CL离聚物含量为0.3，采用Bruggeman修正：
    %   kappaCL=0.3^b_ion*kappa(lambda,T)
    % [FIT-CL] 原程序将b_ion=1.5写死。该指数反映离聚物网络
    % 的曲折度/连通性，附件未给定，因此作为一个可辨识参数。
    %   etaOhmCL=int(im/kappaCL)dx
    lambdaPemLocal = water.lambdaPEM;
    kappaPemLocal = (0.5139.*lambdaPemLocal-0.326).* ...
        exp(1268.*(1/303.15-1./T(g.idx_PEM)));

    lambdaACL = lambdaPemLocal(1);
    % cCL电阻由整层离聚物的平均lambda状态决定。
    kappaACL = 0.3^ionomerBruggemanExponent.* ...
        (0.5139*lambdaACL-0.326).* ...
        exp(1268.*(1/303.15-1./T(g.idx_aCL)));
    kappaCCL = 0.3^ionomerBruggemanExponent.* ...
        (0.5139*lambdaCCL-0.326).* ...
        exp(1268.*(1/303.15-1./T(g.idx_cCL)));

    kappaPem = g.layerThickness(3)/sum(g.dx(g.idx_PEM)./kappaPemLocal);
    if any(kappaPemLocal <= 0) || any(kappaACL <= 0) || ...
            any(kappaCCL <= 0) || ~isfinite(kappaPem)
        etaOhm = Inf;
        etaOhmPEM = Inf;
        etaOhmACL = Inf;
        etaOhmCCL = Inf;
        etaOhmContact = Inf;
    else
        etaOhmPEM = j*sum(g.dx(g.idx_PEM)./kappaPemLocal);
        etaOhmACL = sum(water.imCell(g.idx_aCL).* ...
            g.dx(g.idx_aCL)./kappaACL);
        etaOhmCCL = sum(water.imCell(g.idx_cCL).* ...
            g.dx(g.idx_cCL)./kappaCCL);
        etaOhmContact = j*0.01e-4;
        etaOhm = etaOhmPEM+etaOhmACL+etaOhmCCL+etaOhmContact;
    end

    transportValid = jLim > 0 && j >= 0 && j < jLim;
    if j == 0
        etaCon = 0;
    elseif transportValid
        etaConBase = -R*Tavg/(4*F)*log(1-j/jLim);
        % [Q2-4] 附件1给出首尾电池浓差极化系数10、中间电池1。
        % 这里将其解释为对基础浓差过电位的乘法修正，而不修改
        % 实际极限电流；后续可通过位置分辨电压数据检验这一解释。
        etaCon = heatBoundary.concentrationLossFactor*etaConBase;
    else
        etaCon = Inf;
    end
    Vcell = Erev-etaAct-etaOhm-etaCon;


    %% 5. 等效热容和导热系数

    materialRho = [185,970,2150,970,185].';
    materialCp = [545,240,1050,240,545].';
    materialK = [0.3,0.27,0.24,0.27,0.3].';
    layerRhoCp = materialRho.*materialCp;
    delta = g.layerThickness(:);
    meaArealCapacity = sum(layerRhoCp.*delta);

    bipolarPlateArealCapacityBase = 1980*766*(2e-3+2e-3);
    bipolarPlateArealCapacity = bipolarPlateArealCapacityBase* ...
        heatBoundary.bipolarPlateCapacityScale;
    % [Q2-2] lumped兼容模式下，上层可将端板热容传入
    % extraArealCapacity；默认separate模式的端板由上层独立状态求解。
    endPlateArealCapacity = heatBoundary.extraArealCapacity;
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

    if strcmpi(heatBoundary.mode,'prescribed_flux')
        % [Q2-1] 热通量的符号与内部有限体积面通量一致：
        % 沿+x方向为正。因此左端向环境散热时leftFlux<0，
        % 右端向环境散热时rightFlux>0。
        if ~isscalar(heatBoundary.leftFlux) || ...
                ~isfinite(heatBoundary.leftFlux) || ...
                ~isscalar(heatBoundary.rightFlux) || ...
                ~isfinite(heatBoundary.rightFlux)
            error('thermal_temperature_state_ice:InvalidBoundaryFlux', ...
                '指定热通量模式下左右通量必须为有限标量。');
        end
        qHeatFace(1) = heatBoundary.leftFlux;
        qHeatFace(end) = heatBoundary.rightFlux;
        qOutLeft = -heatBoundary.leftFlux;
        qOutRight = heatBoundary.rightFlux;
        boundarySolidResistance.left = 0;
        boundarySolidResistance.right = 0;
        Rleft = NaN;
        Rright = NaN;
    elseif strcmpi(heatBoundary.mode,'standalone')
        boundarySolidResistance.left = 2e-3/95;
        boundarySolidResistance.right = 2e-3/95;
        Rleft = (g.x(1)-g.xf(1))/kEff(1) + ...
            boundarySolidResistance.left + 1/40;
        Rright = (g.xf(end)-g.x(end))/kEff(end) + ...
            boundarySolidResistance.right + 1/40;
        qOutLeft = (T(1)-Tamb)/Rleft;
        qOutRight = (T(end)-Tamb)/Rright;
        qHeatFace(1) = -qOutLeft;
        qHeatFace(end) = qOutRight;
    else
        error('thermal_temperature_state_ice:InvalidBoundaryMode', ...
            'heatBoundary.mode只能为standalone或prescribed_flux。');
    end
    qCond = (qHeatFace(1:end-1)-qHeatFace(2:end))./g.dx;


    %% 7. 温度状态方程
    % 反应热：qgen=j*(Eth-Vcell)/Lmea
    % [ICE-5] 水的熔化潜热Lf=333600 J/kg：
    %   qPhase=Lf*dmi/dt
    % 因而冻结(dmi>0)放热，融化(dmi<0)吸热。

    qgen = j*(1.48-Vcell)/g.Ltotal;
    qPhase = 333600*water.dmiFull;
    dT = (qCond+qgen+qaux+qPhase)./rhoCp;


    %% 8. 输出诊断信息

    voltage.Erev = Erev;
    voltage.j0Ref = j0Ref;
    voltage.Ea = Ea;
    voltage.j0Base = j0Base;
    voltage.j0 = j0Effective;
    voltage.hydrationDegree = hydrationDegree;
    voltage.hydrationAreaFactor = hydrationAreaFactor;
    voltage.lambdaHydWetEffective = lambdaHydWetEffective;
    voltage.ionomerBruggemanExponent = ionomerBruggemanExponent;
    voltage.fHydDry = fHydDry;
    voltage.lambdaHydOn = lambdaHydOn;
    voltage.lambdaHydWet = lambdaHydWet;
    voltage.nHyd = nHyd;
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
    voltage.concentrationLossFactor = ...
        heatBoundary.concentrationLossFactor;
    voltage.transportValid = transportValid;

    thermalProperties.rhoCpScalar = rhoCpScalar;
    thermalProperties.kEffScalar = kEffScalar;
    thermalProperties.layerRhoCp = layerRhoCp;
    thermalProperties.layerThermalResistance = delta./materialK;
    thermalProperties.thermalResistancePerArea = RthPerArea;
    thermalProperties.Ltotal = g.Ltotal;
    thermalProperties.meaArealCapacity = meaArealCapacity;
    thermalProperties.bipolarPlateArealCapacity = bipolarPlateArealCapacity;
    thermalProperties.bipolarPlateArealCapacityBase = ...
        bipolarPlateArealCapacityBase;
    thermalProperties.bipolarPlateCapacityScale = ...
        heatBoundary.bipolarPlateCapacityScale;
    thermalProperties.endPlateArealCapacity = endPlateArealCapacity;
    thermalProperties.totalArealCapacity = totalArealCapacity;
    thermalProperties.boundarySolidResistance = boundarySolidResistance;
    thermalProperties.boundaryMode = heatBoundary.mode;

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
    out.hydrationDegree = hydrationDegree;
    out.hydrationAreaFactor = hydrationAreaFactor;
    out.gammaIce = gammaIce;

end
