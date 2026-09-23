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
%-------------------------------------------------------------------------------

function [dT,out] = thermal_temperature_state_ice( ...
    T,cH2,cO2,water,j,Tamb,qaux,g,s,gammaIce,j0Ref,Ea)

    %% 1. 输入处理

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
    % [ACT-1] j0Ref和Ea由-20/-25 degC两个零时刻电压联合标定。
    j0Base = j0Ref*exp(-Ea/R*(1/Tavg-1/298.15));
    j0Effective = j0Base*iceAreaFactorForVoltage;
    etaAct = R*Tavg/(0.5*F)*asinh(j/(2*j0Effective));

    kappaPem = (0.5139*water.lambdaMean-0.326)* ...
        exp(1268*(1/303.15-1/Tavg));
    if kappaPem <= 0
        etaOhm = Inf;
    else
        etaOhm = j*(g.layerThickness(3)/kappaPem+0.01e-4);
    end

    transportValid = jLim > 0 && j >= 0 && j < jLim;
    if j == 0
        etaCon = 0;
    elseif transportValid
        etaCon = -R*Tavg/(4*F)*log(1-j/jLim);
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

    bipolarPlateArealCapacity = 1980*766*(2e-3+2e-3);
    endPlateArealCapacity = 0;
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
    Rleft = (g.x(1)-g.xf(1))/kEff(1) + ...
        boundarySolidResistance.left + 1/40;
    Rright = (g.xf(end)-g.x(end))/kEff(end) + ...
        boundarySolidResistance.right + 1/40;
    qOutLeft = (T(1)-Tamb)/Rleft;
    qOutRight = (T(end)-Tamb)/Rright;
    qHeatFace(1) = -qOutLeft;
    qHeatFace(end) = qOutRight;
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
    voltage.iceAreaFactor = iceAreaFactor;
    voltage.etaAct = etaAct;
    voltage.kappaPem = kappaPem;
    voltage.etaOhm = etaOhm;
    voltage.jLim = jLim;
    voltage.etaCon = etaCon;
    voltage.transportValid = transportValid;

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
