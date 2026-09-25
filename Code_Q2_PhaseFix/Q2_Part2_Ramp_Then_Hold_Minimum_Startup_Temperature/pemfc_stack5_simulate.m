%-------------------------------------------------------------------------------
% @function: 求解问题2的5片串联PEMFC自冷启动模型
% @author:   PJ, GPT
% @date:     20260924
% @input:    strategy->q2_current_strategy定义的电流加载策略
%            initialTemperatureC->电堆和环境初始温度 [degC]
%            simOpt->可选求解、热物性和单电池参数
% @output:   result->5片温度、电压、冰体积分数和启动判定
%            model->网格、状态编号、参数和RHS函数句柄
%
% [Q2-1] 相邻电池之间的热通量严格成对出现，保证内部换热守恒。
% [Q2-2] 可比较“端板与端部电池同温”和“端板独立温度节点”。
% [Q2-4] 附件1的首尾/中间电池浓差损失系数分别取10/1。
% [Q2-5] 累积电荷作为额外状态，dq/dt=j，单位为C/cm^2。
%-------------------------------------------------------------------------------

function [result,model] = pemfc_stack5_simulate( ...
    strategy,initialTemperatureC,simOpt)

    %% 1. 输入与默认参数

    if nargin < 1 || isempty(strategy)
        strategy = struct('type','linear','rampRateAcm2s',0.005, ...
            'plateauAcm2',0.3);
    end
    if nargin < 2 || isempty(initialTemperatureC)
        initialTemperatureC = -10;
    end
    if nargin < 3 || isempty(simOpt)
        simOpt = struct();
    end
    if ~isstruct(simOpt) || ~isscalar(simOpt)
        error('pemfc_stack5_simulate:InvalidOptions', ...
            'simOpt必须是标量结构体。');
    end
    if ~isscalar(initialTemperatureC) || ...
            ~isfinite(initialTemperatureC) || initialTemperatureC <= -273.15
        error('pemfc_stack5_simulate:InvalidInitialTemperature', ...
            '初始温度必须是高于绝对零度的有限标量。');
    end

    simOpt = set_default(simOpt,'nCells',5);
    simOpt = set_default(simOpt,'nCellLayer',[6 13 19 19 6]);
    simOpt = set_default(simOpt,'ambientTemperatureC',initialTemperatureC);
    simOpt = set_default(simOpt,'tMax',180);
    simOpt = set_default(simOpt,'outputStep',0.2);
    simOpt = set_default(simOpt,'RelTol',1e-6);
    simOpt = set_default(simOpt,'AbsTol',1e-8);
    simOpt = set_default(simOpt,'MaxStep',0.02);
    simOpt = set_default(simOpt,'verbose',true);
    simOpt = set_default(simOpt,'plot',false);
    simOpt = set_default(simOpt,'progressIntervalS',0);
    simOpt = set_default(simOpt,'progressWallIntervalS',30);
    simOpt = set_default(simOpt,'maxWallTimeS',Inf);
    simOpt = set_default(simOpt,'useJacobianPattern',true);

    % 题目硬约束
    simOpt = set_default(simOpt,'qMaxCcm2',20);
    simOpt = set_default(simOpt,'jMaxAcm2',0.5);
    simOpt = set_default(simOpt,'minimumVoltageV',0.30);
    simOpt = set_default(simOpt,'maximumIceVolumeFraction',0.99);
    % 水+冰总占用率达到98%时先判定孔隙堵塞。2%的剩余孔隙
    % 是给ode15s Newton试探步留的数值安全裕度，避免在事件根被
    % 接受前先跨入扩散系数奇异区。
    simOpt = set_default(simOpt,'maximumPoreOccupancy',0.98);

    % 附件1与题干给定的电堆热参数
    simOpt = set_default(simOpt,'cellAreaM2',25e-4);
    simOpt = set_default(simOpt,'endConvectionCoefficient',40);
    simOpt = set_default(simOpt,'bipolarPlateThicknessM',2e-3);
    simOpt = set_default(simOpt,'bipolarPlateConductivity',95);
    simOpt = set_default(simOpt,'endPlateThicknessM',10e-3);
    simOpt = set_default(simOpt,'endPlateDensity',7900);
    simOpt = set_default(simOpt,'endPlateHeatCapacity',500);
    simOpt = set_default(simOpt,'endPlateConductivity',15);
    simOpt = set_default(simOpt,'endPlateModel','separate');
    simOpt = set_default(simOpt,'bipolarPlateCounting','shared_stack');
    simOpt = set_default(simOpt,'concentrationFactorEnd',10);
    simOpt = set_default(simOpt,'concentrationFactorMiddle',1);

    % 当前Code_ICE单电池标定参数
    simOpt = set_default(simOpt,'j0Ref',0.527582014454);
    simOpt = set_default(simOpt,'Ea',8000);
    simOpt = set_default(simOpt,'tauHyd',10);
    simOpt = set_default(simOpt,'hVaporCathode',0.01);
    simOpt = set_default(simOpt,'hVaporAnode',0.01);
    simOpt = set_default(simOpt,'interfaceRateFactor',0.01);
    simOpt = set_default(simOpt,'fHydDry',0.00705033754574);
    simOpt = set_default(simOpt,'lambdaHydOn',3.5);
    simOpt = set_default(simOpt,'lambdaHydWet',8.53122413758);
    simOpt = set_default(simOpt,'nHyd',3.0);
    simOpt = set_default(simOpt,'ionomerBruggemanExponent',1.77972166611);
    simOpt = set_default(simOpt,'kFreeze',0.05);
    simOpt = set_default(simOpt,'kMelt',1);
    % 问题2默认在冰点附近使用0.2 K窄平滑层，抑制相变开关抖振。
    simOpt = set_default(simOpt,'phaseTransitionWidthK',0.2);
    simOpt = set_default(simOpt,'directFreezeFraction',0.01);
    simOpt = set_default(simOpt,'freezeNucleationLambdaFraction',0.5);
    simOpt = set_default(simOpt,'gammaIce',0.1);
    if ~isfield(simOpt,'init') || isempty(simOpt.init)
        simOpt.init = struct('lambdaCCL0',2.0);
    end

    validate_stack_options(simOpt);
    if simOpt.nCells ~= 5
        error('pemfc_stack5_simulate:FiveCellsRequired', ...
            '问题2固定为5片串联单电池。');
    end


    %% 2. 建立五个相同的单电池状态

    initialKelvin = initialTemperatureC+273.15;
    initOpt = simOpt.init;
    initOpt.T0 = initialKelvin;
    [g,s,x0Cell,init] = pemfc_setup_ice(simOpt.nCellLayer,initOpt);

    nCells = simOpt.nCells;
    cellStateIndex = cell(nCells,1);
    for k = 1:nCells
        cellStateIndex{k} = (k-1)*s.Nx+(1:s.Nx);
    end
    chargeIndex = nCells*s.Nx+1;
    if strcmpi(char(simOpt.endPlateModel),'separate')
        endPlateIndex = chargeIndex+(1:2);
        x0 = [repmat(x0Cell,nCells,1);0;initialKelvin;initialKelvin];
    else
        endPlateIndex = [];
        x0 = [repmat(x0Cell,nCells,1);0];
    end


    %% 3. 电堆热耦合参数
    % 对式(2)两边同除以面积A，可得进入一维能量方程的面热通量：
    %   q''_k,k+1=(k_eff/delta)*(T_k-T_k+1)。
    % 默认将一个单电池节距作为两个平均温度节点间距，
    % 其热阻由MEA五层和阴/阳极双极板串联得到。

    materialConductivity = [0.3,0.27,0.24,0.27,0.3];
    meaThermalResistance = sum( ...
        g.layerThickness./materialConductivity);
    if strcmpi(char(simOpt.bipolarPlateCounting),'shared_stack')
        % 5片堆共6块双极板：端部电池分配1.5块，
        % 中间电池分配1块，合计6块。基准单片热容含2块。
        bipolarPlateCapacityScale = [0.75;0.5;0.5;0.5;0.75];
        intercellBipolarPlateCount = 1;
    else
        bipolarPlateCapacityScale = ones(nCells,1);
        intercellBipolarPlateCount = 2;
    end
    bipolarPlateThermalResistance = ...
        intercellBipolarPlateCount*simOpt.bipolarPlateThicknessM/ ...
        simOpt.bipolarPlateConductivity;
    intercellThermalResistance = ...
        meaThermalResistance+bipolarPlateThermalResistance;
    intercellDistance = g.Ltotal+intercellBipolarPlateCount* ...
        simOpt.bipolarPlateThicknessM;
    stackEffectiveConductivity = ...
        intercellDistance/intercellThermalResistance;
    intercellConductancePerArea = 1/intercellThermalResistance;

    stackThermal.cellAreaM2 = simOpt.cellAreaM2;
    stackThermal.ambientTemperatureK = ...
        simOpt.ambientTemperatureC+273.15;
    stackThermal.endConvectionCoefficient = ...
        simOpt.endConvectionCoefficient;
    stackThermal.endPlateArealCapacity = ...
        simOpt.endPlateDensity*simOpt.endPlateHeatCapacity* ...
        simOpt.endPlateThicknessM;
    stackThermal.endPlateModel = char(simOpt.endPlateModel);
    stackThermal.bipolarPlateCounting = ...
        char(simOpt.bipolarPlateCounting);
    stackThermal.bipolarPlateCapacityScale = ...
        bipolarPlateCapacityScale;
    stackThermal.endPlateConductivity = simOpt.endPlateConductivity;
    stackThermal.cellToEndPlateThermalResistance = ...
        0.5*meaThermalResistance + ...
        simOpt.bipolarPlateThicknessM/simOpt.bipolarPlateConductivity + ...
        0.5*simOpt.endPlateThicknessM/simOpt.endPlateConductivity;
    stackThermal.cellToEndPlateConductancePerArea = ...
        1/stackThermal.cellToEndPlateThermalResistance;
    stackThermal.endPlateExternalResistance = ...
        0.5*simOpt.endPlateThicknessM/simOpt.endPlateConductivity + ...
        1/simOpt.endConvectionCoefficient;
    stackThermal.endPlateExternalConductancePerArea = ...
        1/stackThermal.endPlateExternalResistance;
    stackThermal.intercellDistance = intercellDistance;
    stackThermal.meaThermalResistance = meaThermalResistance;
    stackThermal.bipolarPlateThermalResistance = ...
        bipolarPlateThermalResistance;
    stackThermal.intercellThermalResistance = ...
        intercellThermalResistance;
    stackThermal.stackEffectiveConductivity = ...
        stackEffectiveConductivity;
    stackThermal.intercellConductancePerArea = ...
        intercellConductancePerArea;

    ice = struct('kFreeze',simOpt.kFreeze,'kMelt',simOpt.kMelt, ...
        'phaseTransitionWidthK',simOpt.phaseTransitionWidthK, ...
        'directFreezeFraction',simOpt.directFreezeFraction, ...
        'freezeNucleationLambdaFraction', ...
        simOpt.freezeNucleationLambdaFraction, ...
        'gammaIce',simOpt.gammaIce);
    hydration = struct('tauHyd',simOpt.tauHyd);
    transport = struct('hVaporCathode',simOpt.hVaporCathode, ...
        'hVaporAnode',simOpt.hVaporAnode, ...
        'interfaceRateFactor',simOpt.interfaceRateFactor);
    activation = struct('j0Ref',simOpt.j0Ref,'Ea',simOpt.Ea, ...
        'fHydDry',simOpt.fHydDry,'lambdaHydOn',simOpt.lambdaHydOn, ...
        'lambdaHydWet',simOpt.lambdaHydWet,'nHyd',simOpt.nHyd, ...
        'ionomerBruggemanExponent',simOpt.ionomerBruggemanExponent);
    constraints = struct('qMaxCcm2',simOpt.qMaxCcm2, ...
        'jMaxAcm2',simOpt.jMaxAcm2, ...
        'minimumVoltageV',simOpt.minimumVoltageV, ...
        'maximumIceVolumeFraction',simOpt.maximumIceVolumeFraction, ...
        'maximumPoreOccupancy',simOpt.maximumPoreOccupancy);


    %% 4. 调用刚性求解器并用事件判定启动结果

    tOutput = (0:simOpt.outputStep:simOpt.tMax).';
    if tOutput(end) < simOpt.tMax
        tOutput(end+1,1) = simOpt.tMax;
    end

    nonNegativeCell = [s.idx.T(:);s.idx.H2(:);s.idx.O2(:);s.idx.mw(:); ...
        s.idx.mi(:);s.idx.lambdaCCL(:)];
    nonNegativeStates = chargeIndex;
    for k = 1:nCells
        nonNegativeStates = [nonNegativeStates; ...
            (k-1)*s.Nx+nonNegativeCell]; %#ok<AGROW>
    end
    nonNegativeStates = sort(nonNegativeStates);

    rhs = @(t,x) stack_rhs_q2(t,x,strategy,ice,hydration, ...
        activation,transport,stackThermal,simOpt,g,s,cellStateIndex, ...
        chargeIndex,endPlateIndex);
    eventFunction = @(t,x) stack_events_q2(t,x,rhs,constraints, ...
        chargeIndex);
    odeOpt = odeset('RelTol',simOpt.RelTol,'AbsTol',simOpt.AbsTol, ...
        'MaxStep',simOpt.MaxStep,'NonNegative',nonNegativeStates, ...
        'Events',eventFunction);
    jacobianPattern = [];
    if simOpt.useJacobianPattern
        % [NUM-JAC] 单片内部取保守稠密块，只把真实存在的相邻片
        % 导热耦合和端板耦合放到块外。它不改变方程，只减少ode15s
        % 为数值雅可比分组扰动状态时重复调用整套五片RHS的次数。
        jacobianPattern = build_stack_jacobian_pattern( ...
            s,cellStateIndex,chargeIndex,endPlateIndex,numel(x0));
        odeOpt = odeset(odeOpt,'JPattern',jacobianPattern);
    end
    if simOpt.progressIntervalS > 0
        progressFunction = @(t,y,flag) stack_progress_output( ...
            t,y,flag,simOpt.progressIntervalS,simOpt.tMax);
        odeOpt = odeset(odeOpt,'OutputFcn',progressFunction);
    end
    stack_wallclock_heartbeat('reset',0,simOpt.tMax, ...
        simOpt.progressWallIntervalS,simOpt.maxWallTimeS);

    if simOpt.verbose
        fprintf('\n开始计算5片电堆：T0=%.3f degC，策略=%s\n', ...
            initialTemperatureC,char(strategy.type));
        if strcmpi(char(strategy.type),'constant')
            fprintf('恒定电流密度=%.6g A/cm^2\n', ...
                strategy.currentAcm2);
        end
        fprintf(['相邻节距=%.6g m，k_eff=%.6g W/(m K)，', ...
            'G/A=%.6g W/(m^2 K)\n'],intercellDistance, ...
            stackEffectiveConductivity,intercellConductancePerArea);
        fprintf('端板模式=%s，双极板计数=%s\n', ...
            stackThermal.endPlateModel,stackThermal.bipolarPlateCounting);
    end

    % 恒流策略可能在t=0就使电压或电流越界。MATLAB事件
    % 只检测“穿过零点”，因此必须在积分前单独判定初值，
    % 否则初始已经低于0.30 V的工况会被错误地继续积分。
    [~,initialOut] = rhs(0,x0);
    initialEventIndex = [];
    if initialOut.jAcm2 > constraints.jMaxAcm2+1e-12
        initialEventIndex = 5;
    elseif min(initialOut.Vcell) < constraints.minimumVoltageV-1e-10
        initialEventIndex = 3;
    elseif max(initialOut.maxIceVolumeFractionCell) >= ...
            constraints.maximumIceVolumeFraction
        initialEventIndex = 4;
    elseif max(initialOut.maxPoreOccupancyCell) >= ...
            constraints.maximumPoreOccupancy
        initialEventIndex = 6;
    end

    tic;
    if isempty(initialEventIndex)
        % 恢复原目录的调用方式：全时域一次ode15s调用，不按阶梯重启。
        [tSol,xSol,tEvent,~,eventIndex] = ...
            ode15s(rhs,tOutput,x0,odeOpt);
    else
        tSol = 0;
        xSol = x0.';
        tEvent = 0;
        eventIndex = initialEventIndex;
    end
    cpuTime = toc;
    stack_wallclock_heartbeat('done',tSol(end),simOpt.tMax, ...
        simOpt.progressWallIntervalS,simOpt.maxWallTimeS);


    %% 5. 重算各片代数输出

    nTime = numel(tSol);
    TavgC = zeros(nTime,nCells);
    Vcell = zeros(nTime,nCells);
    maxIceVolumeFractionCell = zeros(nTime,nCells);
    % [DIAG-ICE] 冰填充率sIce=epsIce/eps0，表示原始孔隙中
    % 被冰占据的比例。与题目式(9)的总体积冰分数分开保存。
    maxIceFillingRatioCell = zeros(nTime,nCells);
    maxPoreOccupancyCell = zeros(nTime,nCells);
    iceVolumeFractionCCL = zeros(nTime,nCells);
    lambdaCCL = zeros(nTime,nCells);
    limitingCurrentAcm2 = zeros(nTime,nCells);
    maximumLiquidSaturationCell = zeros(nTime,nCells);
    etaActCell = zeros(nTime,nCells);
    etaOhmCell = zeros(nTime,nCells);
    etaOhmPEMCell = zeros(nTime,nCells);
    etaOhmACLCell = zeros(nTime,nCells);
    etaOhmCCLCell = zeros(nTime,nCells);
    etaConCell = zeros(nTime,nCells);
    lambdaPEMMean = zeros(nTime,nCells);
    interfaceHeatFlux = zeros(nTime,nCells-1);
    endConvectionHeatFlux = zeros(nTime,2);
    endPlateTemperatureC = zeros(nTime,2);
    currentDensityAcm2 = zeros(nTime,1);

    for n = 1:nTime
        [~,out] = rhs(tSol(n),xSol(n,:).');
        TavgC(n,:) = out.Tavg.'-273.15;
        Vcell(n,:) = out.Vcell.';
        maxIceVolumeFractionCell(n,:) = ...
            out.maxIceVolumeFractionCell.';
        maxIceFillingRatioCell(n,:) = ...
            out.maxIceFillingRatioCell.';
        maxPoreOccupancyCell(n,:) = out.maxPoreOccupancyCell.';
        iceVolumeFractionCCL(n,:) = out.iceVolumeFractionCCL.';
        lambdaCCL(n,:) = out.lambdaCCL.';
        for k = 1:nCells
            limitingCurrentAcm2(n,k) = out.cell{k}.thermal.jLim/1e4;
            maximumLiquidSaturationCell(n,k) = max( ...
                out.cell{k}.water.liquidSaturationCapillary);
            etaActCell(n,k) = out.cell{k}.thermal.voltage.etaAct;
            etaOhmCell(n,k) = out.cell{k}.thermal.voltage.etaOhm;
            etaOhmPEMCell(n,k) = out.cell{k}.thermal.voltage.etaOhmPEM;
            etaOhmACLCell(n,k) = out.cell{k}.thermal.voltage.etaOhmACL;
            etaOhmCCLCell(n,k) = out.cell{k}.thermal.voltage.etaOhmCCL;
            etaConCell(n,k) = out.cell{k}.thermal.voltage.etaCon;
            lambdaPEMMean(n,k) = out.cell{k}.water.lambdaMean;
        end
        interfaceHeatFlux(n,:) = out.interfaceHeatFlux.';
        endConvectionHeatFlux(n,:) = ...
            [out.leftConvectionOutflux,out.rightConvectionOutflux];
        endPlateTemperatureC(n,:) = out.endPlateTemperatureK.'-273.15;
        currentDensityAcm2(n) = out.j/1e4;
    end


    %% 6. 根据题目式(4)-(9)统一判定成功与失败

    chargeCcm2 = xSol(:,chargeIndex);
    minimumCellVoltage = min(Vcell(:));
    maximumIceVolumeFraction = max(maxIceVolumeFractionCell(:));
    maximumIceFillingRatio = max(maxIceFillingRatioCell(:));
    maximumPoreOccupancy = max(maxPoreOccupancyCell(:));
    maximumCurrentDensity = max(currentDensityAcm2);
    temperatureReached = all(TavgC(end,:) >= -1e-7);
    voltageValid = minimumCellVoltage >= ...
        constraints.minimumVoltageV-1e-8;
    iceValid = maximumIceVolumeFraction < ...
        constraints.maximumIceVolumeFraction;
    poreOccupancyValid = maximumPoreOccupancy < ...
        constraints.maximumPoreOccupancy;
    currentValid = maximumCurrentDensity <= ...
        constraints.jMaxAcm2+1e-12;
    chargeValid = chargeCcm2(end) <= constraints.qMaxCcm2+1e-8;
    successEvent = any(eventIndex == 1);
    requestedEndReached = tSol(end) >= tOutput(end)- ...
        100*eps(max(abs(tOutput(end)),1));
    solverTerminatedUnexpectedly = isempty(eventIndex) && ...
        ~requestedEndReached;
    success = ~solverTerminatedUnexpectedly && successEvent && ...
        temperatureReached && voltageValid && ...
        iceValid && poreOccupancyValid && currentValid && chargeValid;

    if success
        successTime = tEvent(find(eventIndex == 1,1));
        stopReason = '五片平均温度同时达到0 degC';
    else
        successTime = NaN;
        if solverTerminatedUnexpectedly
            stopReason = '求解器提前终止（通常为相变切换导致步长过小）';
        elseif any(eventIndex == 2)
            stopReason = '累积电荷达到上限';
        elseif any(eventIndex == 3)
            stopReason = '单片电压下降到0.30 V';
        elseif any(eventIndex == 4)
            stopReason = '局部冰体积分数达到上限';
        elseif any(eventIndex == 5) || ~currentValid
            stopReason = '电流密度达到或超过上限';
        elseif any(eventIndex == 6) || ~poreOccupancyValid
            stopReason = '局部孔隙被水和冰占用至上限';
        elseif successEvent
            stopReason = '到0 degC时存在历史或加载约束违反';
        else
            stopReason = '在最大仿真时间内未启动成功';
        end
    end

    [~,minimumVoltageLinearIndex] = min(Vcell(:));
    [minimumVoltageTimeIndex,minimumVoltageCellIndex] = ...
        ind2sub(size(Vcell),minimumVoltageLinearIndex);
    [~,maximumIceLinearIndex] = max(maxIceVolumeFractionCell(:));
    [maximumIceTimeIndex,maximumIceCellIndex] = ...
        ind2sub(size(maxIceVolumeFractionCell),maximumIceLinearIndex);
    [~,coldestCellIndex] = min(TavgC(end,:));

    result.initialTemperatureC = initialTemperatureC;
    result.strategy = strategy;
    result.t = tSol;
    result.x = xSol;
    result.currentDensityAcm2 = currentDensityAcm2;
    result.chargeCcm2 = chargeCcm2;
    result.TavgC = TavgC;
    result.Vcell = Vcell;
    result.stackVoltage = sum(Vcell,2);
    result.maxIceVolumeFractionCell = maxIceVolumeFractionCell;
    result.maximumIceVolumeFraction = maximumIceVolumeFraction;
    result.maxIceFillingRatioCell = maxIceFillingRatioCell;
    result.maximumIceFillingRatio = maximumIceFillingRatio;
    result.maxPoreOccupancyCell = maxPoreOccupancyCell;
    result.maximumPoreOccupancy = maximumPoreOccupancy;
    result.iceVolumeFractionCCL = iceVolumeFractionCCL;
    result.lambdaCCL = lambdaCCL;
    % [DIAG-1] 用于区分“冰堵”、“液水淹没”和“CL质子电阻”。
    result.limitingCurrentDensityAcm2 = limitingCurrentAcm2;
    result.maximumLiquidSaturationCell = maximumLiquidSaturationCell;
    result.etaActCell = etaActCell;
    result.etaOhmCell = etaOhmCell;
    result.etaOhmPEMCell = etaOhmPEMCell;
    result.etaOhmACLCell = etaOhmACLCell;
    result.etaOhmCCLCell = etaOhmCCLCell;
    result.etaConCell = etaConCell;
    result.lambdaPEMMean = lambdaPEMMean;
    result.interfaceHeatFluxWm2 = interfaceHeatFlux;
    result.interfaceHeatPowerW = ...
        interfaceHeatFlux*simOpt.cellAreaM2;
    result.endConvectionHeatFluxWm2 = endConvectionHeatFlux;
    result.endPlateTemperatureC = endPlateTemperatureC;
    result.success = success;
    result.successTimeS = successTime;
    result.stopTimeS = tSol(end);
    result.stopReason = stopReason;
    result.chargeUsedCcm2 = chargeCcm2(end);
    result.maximumCurrentDensityAcm2 = maximumCurrentDensity;
    result.minimumCellVoltageV = minimumCellVoltage;
    result.minimumVoltageCellIndex = minimumVoltageCellIndex;
    result.minimumVoltageTimeS = tSol(minimumVoltageTimeIndex);
    result.maximumIceCellIndex = maximumIceCellIndex;
    result.maximumIceTimeS = tSol(maximumIceTimeIndex);
    result.coldestCellIndexAtStop = coldestCellIndex;
    result.constraints = constraints;
    result.constraintCheck = struct('temperatureReached',temperatureReached, ...
        'voltageValid',voltageValid,'iceValid',iceValid, ...
        'poreOccupancyValid',poreOccupancyValid, ...
        'currentValid',currentValid,'chargeValid',chargeValid);
    result.solverTerminatedUnexpectedly = solverTerminatedUnexpectedly;
    result.event.time = tEvent;
    result.event.index = eventIndex;
    result.cpuTime = cpuTime;
    result.options = simOpt;
    result.integrationScheme = 'original_ode15s_single_call_v1';
    result.stackThermal = stackThermal;

    model.g = g;
    model.s = s;
    model.init = init;
    model.x0Cell = x0Cell;
    model.x0 = x0;
    model.cellStateIndex = cellStateIndex;
    model.chargeIndex = chargeIndex;
    model.endPlateIndex = endPlateIndex;
    model.strategy = strategy;
    model.options = simOpt;
    model.stackThermal = stackThermal;
    model.rhs = rhs;
    model.jacobianPattern = jacobianPattern;

    if simOpt.verbose
        fprintf('终止原因：%s\n',stopReason);
        fprintf(['时间=%.4f s，电荷=%.4f C/cm^2，最低电压=%.4f V，', ...
            '最大冰体积分数=%.5f，最大冰填充率=%.5f，', ...
            '最大孔隙占用率=%.5f\n'], ...
            tSol(end),chargeCcm2(end),minimumCellVoltage, ...
            maximumIceVolumeFraction,maximumIceFillingRatio, ...
            maximumPoreOccupancy);
    end


    %% 7. 可选绘图

    if simOpt.plot
        plot_stack_result(result);
    end

end


%% ========================================================================
% 局部函数：五片电堆RHS
% ========================================================================

function [dx,out] = stack_rhs_q2(t,x,strategy,ice,hydration,activation, ...
    transport,stackThermal,simOpt,g,s,cellStateIndex,chargeIndex, ...
    endPlateIndex)

    timedOut = stack_wallclock_heartbeat('update',t,simOpt.tMax, ...
        simOpt.progressWallIntervalS,simOpt.maxWallTimeS);
    if timedOut
        error('pemfc_stack5_simulate:WallTimeLimit', ...
            'Single black-box simulation exceeded %.1f s wall time.', ...
            simOpt.maxWallTimeS);
    end
    nCells = numel(cellStateIndex);
    [j,jAcm2] = q2_current_strategy(t,strategy);
    if ~isscalar(j) || ~isfinite(j) || j < 0
        error('pemfc_stack5_simulate:InvalidCurrent', ...
            '加载策略必须返回非负有限标量电流。');
    end

    Tavg = zeros(nCells,1);
    for k = 1:nCells
        xCell = x(cellStateIndex{k});
        Tcell = xCell(s.idx.T);
        Tavg(k) = sum(Tcell.*g.dx)/sum(g.dx);
    end
    interfaceHeatFlux = stackThermal.intercellConductancePerArea.* ...
        (Tavg(1:end-1)-Tavg(2:end));

    if strcmpi(stackThermal.endPlateModel,'separate')
        endPlateTemperatureK = x(endPlateIndex);
        leftCellToEndPlateFlux = ...
            stackThermal.cellToEndPlateConductancePerArea* ...
            (Tavg(1)-endPlateTemperatureK(1));
        rightCellToEndPlateFlux = ...
            stackThermal.cellToEndPlateConductancePerArea* ...
            (Tavg(end)-endPlateTemperatureK(2));
        leftConvectionOutflux = ...
            stackThermal.endPlateExternalConductancePerArea* ...
            (endPlateTemperatureK(1)-stackThermal.ambientTemperatureK);
        rightConvectionOutflux = ...
            stackThermal.endPlateExternalConductancePerArea* ...
            (endPlateTemperatureK(2)-stackThermal.ambientTemperatureK);
        leftBoundaryOutflux = leftCellToEndPlateFlux;
        rightBoundaryOutflux = rightCellToEndPlateFlux;
    else
        endPlateTemperatureK = [Tavg(1);Tavg(end)];
        leftConvectionOutflux = stackThermal.endConvectionCoefficient* ...
            (Tavg(1)-stackThermal.ambientTemperatureK);
        rightConvectionOutflux = stackThermal.endConvectionCoefficient* ...
            (Tavg(end)-stackThermal.ambientTemperatureK);
        leftBoundaryOutflux = leftConvectionOutflux;
        rightBoundaryOutflux = rightConvectionOutflux;
    end

    dx = zeros(size(x));
    cellOut = cell(nCells,1);
    for k = 1:nCells
        if k == 1
            leftFlux = -leftBoundaryOutflux;
        else
            leftFlux = interfaceHeatFlux(k-1);
        end
        if k == nCells
            rightFlux = rightBoundaryOutflux;
        else
            rightFlux = interfaceHeatFlux(k);
        end

        isEndCell = k == 1 || k == nCells;
        heatBoundary.mode = 'prescribed_flux';
        heatBoundary.leftFlux = leftFlux;
        heatBoundary.rightFlux = rightFlux;
        if strcmpi(stackThermal.endPlateModel,'lumped')
            heatBoundary.extraArealCapacity = ...
                double(isEndCell)*stackThermal.endPlateArealCapacity;
        else
            heatBoundary.extraArealCapacity = 0;
        end
        heatBoundary.bipolarPlateCapacityScale = ...
            stackThermal.bipolarPlateCapacityScale(k);
        if isEndCell
            heatBoundary.concentrationLossFactor = ...
                simOpt.concentrationFactorEnd;
        else
            heatBoundary.concentrationLossFactor = ...
                simOpt.concentrationFactorMiddle;
        end

        [dxCell,cellOut{k}] = cell_rhs_q2(t,x(cellStateIndex{k}),j, ...
            stackThermal.ambientTemperatureK,heatBoundary,ice,hydration, ...
            activation,transport,g,s);
        dx(cellStateIndex{k}) = dxCell;
    end
    % [Q2-5] A/cm^2 = C/(s cm^2)。
    dx(chargeIndex) = jAcm2;
    if strcmpi(stackThermal.endPlateModel,'separate')
        dx(endPlateIndex(1)) = (leftCellToEndPlateFlux- ...
            leftConvectionOutflux)/stackThermal.endPlateArealCapacity;
        dx(endPlateIndex(2)) = (rightCellToEndPlateFlux- ...
            rightConvectionOutflux)/stackThermal.endPlateArealCapacity;
    end

    if nargout > 1
        out.t = t;
        out.j = j;
        out.jAcm2 = jAcm2;
        out.Tavg = Tavg;
        out.Vcell = zeros(nCells,1);
        out.maxIceVolumeFractionCell = zeros(nCells,1);
        out.maxIceFillingRatioCell = zeros(nCells,1);
        out.maxPoreOccupancyCell = zeros(nCells,1);
        out.iceVolumeFractionCCL = zeros(nCells,1);
        out.lambdaCCL = zeros(nCells,1);
        out.cell = cellOut;
        for k = 1:nCells
            out.Vcell(k) = cellOut{k}.Vcell;
            out.maxIceVolumeFractionCell(k) = ...
                max(cellOut{k}.water.iceVolumeFraction);
            out.maxIceFillingRatioCell(k) = ...
                max(cellOut{k}.water.sIce);
            out.maxPoreOccupancyCell(k) = ...
                max(cellOut{k}.water.poreOccupancyFraction);
            out.iceVolumeFractionCCL(k) = sum( ...
                cellOut{k}.water.iceVolumeFraction(g.idx_cCL).* ...
                g.dx(g.idx_cCL))/g.layerThickness(4);
            out.lambdaCCL(k) = cellOut{k}.lambdaCCL;
        end
        out.interfaceHeatFlux = interfaceHeatFlux;
        out.leftConvectionOutflux = leftConvectionOutflux;
        out.rightConvectionOutflux = rightConvectionOutflux;
        out.endPlateTemperatureK = endPlateTemperatureK;
        out.leftBoundaryOutflux = leftBoundaryOutflux;
        out.rightBoundaryOutflux = rightBoundaryOutflux;
    end

end


%% ========================================================================
% 局部函数：单片RHS，物理内核仍只调用水-冰、热、气体三个模块
% ========================================================================

function [dx,out] = cell_rhs_q2(t,x,j,Tamb,heatBoundary,ice,hydration, ...
    activation,transport,g,s)

    TState = x(s.idx.T);
    cH2State = x(s.idx.H2);
    cO2State = x(s.idx.O2);
    mwState = x(s.idx.mw);
    miState = x(s.idx.mi);
    lambdaCCLState = x(s.idx.lambdaCCL);

    % 仅投影求解器Newton试探的非物理低温，正常解不受影响。
    T = max(TState,150);
    % 防止Newton试探点在分压对数中产生log(0)。
    cH2 = max(cH2State,1e-12);
    cO2 = max(cO2State,1e-12);
    mw = max(mwState,0);
    mi = max(miState,0);
    lambdaCCL = max(lambdaCCLState,0);

    [dmw,dmi,dlambdaCCL,water] = water_ice_state_ice( ...
        T,mw,mi,lambdaCCL,j,ice.kFreeze,ice.kMelt, ...
        ice.phaseTransitionWidthK, ...
        ice.directFreezeFraction,ice.freezeNucleationLambdaFraction, ...
        hydration.tauHyd, ...
        transport.hVaporCathode,transport.hVaporAnode, ...
        transport.interfaceRateFactor,g,s);
    [dT,thermal] = thermal_temperature_state_ice( ...
        T,cH2,cO2,water,j,Tamb,zeros(g.N,1),g,s,ice.gammaIce, ...
        activation.j0Ref,activation.Ea,activation.fHydDry, ...
        activation.lambdaHydOn,activation.lambdaHydWet, ...
        activation.nHyd,activation.ionomerBruggemanExponent,heatBoundary);
    [dcH2,dcO2,gas] = gas_transport_state_ice( ...
        T,cH2,cO2,mw,dmw,dT,water,thermal,j,g,s);

    dx = zeros(s.Nx,1);
    dx(s.idx.T) = dT;
    dx(s.idx.H2) = dcH2;
    dx(s.idx.O2) = dcO2;
    dx(s.idx.mw) = dmw;
    dx(s.idx.mi) = dmi;
    dx(s.idx.lambdaCCL) = dlambdaCCL;

    if nargout > 1
        out.t = t;
        out.Vcell = thermal.Vcell;
        out.Tavg = thermal.Tavg;
        out.lambdaCCL = water.lambdaCCL;
        out.water = water;
        out.thermal = thermal;
        out.gas = gas;
    end

end


%% ========================================================================
% 局部函数：启动成功和三个失败边界
% ========================================================================

function [value,isTerminal,direction] = stack_events_q2( ...
    t,x,rhs,constraints,chargeIndex)

    [~,out] = rhs(t,x);
    minimumTemperatureC = min(out.Tavg)-273.15;
    chargeMargin = constraints.qMaxCcm2-x(chargeIndex);
    voltageMargin = min(out.Vcell)-constraints.minimumVoltageV;
    iceMargin = constraints.maximumIceVolumeFraction- ...
        max(out.maxIceVolumeFractionCell);
    poreOccupancyMargin = constraints.maximumPoreOccupancy- ...
        max(out.maxPoreOccupancyCell);
    % j=jmax本身是允许的；事件面略高于上限，避免
    % 平台刚好取jmax时被误判为失败。
    currentEventLimit = constraints.jMaxAcm2+ ...
        max(1e-10,1e-8*constraints.jMaxAcm2);
    currentMargin = currentEventLimit-out.jAcm2;

    value = [minimumTemperatureC;chargeMargin;voltageMargin;iceMargin; ...
        currentMargin;poreOccupancyMargin];
    isTerminal = [1;1;1;1;1;1];
    direction = [1;-1;-1;-1;-1;-1];

end


function output = set_default(input,fieldName,defaultValue)

    output = input;
    if ~isfield(output,fieldName) || isempty(output.(fieldName))
        output.(fieldName) = defaultValue;
    end

end


function status = stack_progress_output(t,~,flag,intervalS,tMax)

    % 仅按仿真物理时间间隔打印，不在每个内部步长刷屏。
    persistent nextPrintTime
    persistent lastSimulationTime
    persistent wallClock
    status = 0;

    if strcmp(flag,'init')
        nextPrintTime = intervalS*(floor(t(1)/intervalS)+1);
        lastSimulationTime = t(1);
        wallClock = tic;
        print_progress(lastSimulationTime,tMax,0);
    elseif isempty(flag)
        currentTime = t(end);
        lastSimulationTime = currentTime;
        if isempty(nextPrintTime)
            nextPrintTime = intervalS;
        end
        if currentTime+10*eps(max(abs(currentTime),1)) >= nextPrintTime
            print_progress(currentTime,tMax,toc(wallClock));
            nextPrintTime = nextPrintTime+intervalS;
            while nextPrintTime <= ...
                    currentTime+10*eps(max(abs(currentTime),1))
                nextPrintTime = nextPrintTime+intervalS;
            end
        end
    elseif strcmp(flag,'done')
        if ~isempty(lastSimulationTime)
            print_progress(lastSimulationTime,tMax,toc(wallClock));
            fprintf('当前积分区间结束。\n');
        end
        nextPrintTime = [];
        lastSimulationTime = [];
        wallClock = [];
    end

end


function print_progress(simulationTime,tMax,wallSeconds)

    percent = 100*min(max(simulationTime/max(tMax,eps),0),1);
    fprintf(1,['仿真进度：t = %8.3f / %8.3f s  ', ...
        '(%6.2f%%)，实际耗时 = %8.2f s\n'], ...
        simulationTime,tMax,percent,wallSeconds);
    % 使桌面MATLAB命令窗口及时刷新，而不是等本次仿真结束后集中显示。
    drawnow limitrate;

end


function timedOut = stack_wallclock_heartbeat( ...
    action,simulationTime,tMax,intervalS,maxWallTimeS)

    % 刚性求解时物理时间可能很久不推进。每隔一段墙钟时间打印心跳，
    % 使用户能区分“仍在Newton迭代”和“MATLAB已经停止”。
    persistent wallClock lastPrintWall callCount
    timedOut = false;
    if strcmp(action,'reset')
        wallClock = tic;
        lastPrintWall = 0;
        callCount = 0;
        return;
    elseif strcmp(action,'done')
        wallClock = [];
        lastPrintWall = [];
        callCount = [];
        return;
    end
    if intervalS<=0 || isempty(wallClock)
        return;
    end
    callCount = callCount+1;
    % 避免在每次RHS调用时都查询墙钟；每500次检查一次即可。
    if mod(callCount,500)~=0
        return;
    end
    elapsed = toc(wallClock);
    timedOut = elapsed >= maxWallTimeS;
    if elapsed-lastPrintWall >= intervalS
        fprintf(1,['计算心跳：求解器当前尝试 t = %.6f s，', ...
            '目标上限 = %.3f s，实际耗时 = %.1f s\n'], ...
            simulationTime,tMax,elapsed);
        drawnow limitrate;
        lastPrintWall = elapsed;
    end

end


function pattern = build_stack_jacobian_pattern( ...
    s,cellStateIndex,chargeIndex,endPlateIndex,nState)

    % 单片RHS内各物理场通过局部电压、热源和有效参数相互耦合。
    % 这里用稠密单片块作为安全上界，不冒险漏掉非零偏导。
    nCells = numel(cellStateIndex);
    rowParts = cell(3*nCells+4,1);
    columnParts = cell(3*nCells+4,1);
    partIndex = 0;
    for k = 1:nCells
        currentCell = cellStateIndex{k};
        [rowGrid,columnGrid] = ndgrid(currentCell,currentCell);
        partIndex = partIndex+1;
        rowParts{partIndex} = rowGrid(:);
        columnParts{partIndex} = columnGrid(:);

        % 跨片只通过相邻单片温度场的平均温度发生导热耦合。
        temperatureRows = currentCell(s.idx.T);
        if k > 1
            previousTemperature = cellStateIndex{k-1}(s.idx.T);
            [rowGrid,columnGrid] = ndgrid( ...
                temperatureRows,previousTemperature);
            partIndex = partIndex+1;
            rowParts{partIndex} = rowGrid(:);
            columnParts{partIndex} = columnGrid(:);
        end
        if k < nCells
            nextTemperature = cellStateIndex{k+1}(s.idx.T);
            [rowGrid,columnGrid] = ndgrid(temperatureRows,nextTemperature);
            partIndex = partIndex+1;
            rowParts{partIndex} = rowGrid(:);
            columnParts{partIndex} = columnGrid(:);
        end
    end

    % 电荷状态只依赖当前电流策略和自身；若策略以后扩展为
    % 状态反馈，单片块仍然是保守的，届时只需补充相应列。
    partIndex = partIndex+1;
    rowParts{partIndex} = chargeIndex;
    columnParts{partIndex} = chargeIndex;

    if ~isempty(endPlateIndex)
        firstTemperature = cellStateIndex{1}(s.idx.T);
        lastTemperature = cellStateIndex{end}(s.idx.T);
        patternBlocks = {firstTemperature,endPlateIndex(1); ...
            lastTemperature,endPlateIndex(2); ...
            endPlateIndex(1),firstTemperature; ...
            endPlateIndex(2),lastTemperature; ...
            endPlateIndex,endPlateIndex};
        for blockIndex = 1:size(patternBlocks,1)
            [rowGrid,columnGrid] = ndgrid(patternBlocks{blockIndex,1}, ...
                patternBlocks{blockIndex,2});
            partIndex = partIndex+1;
            rowParts{partIndex} = rowGrid(:);
            columnParts{partIndex} = columnGrid(:);
        end
    end

    patternRows = vertcat(rowParts{1:partIndex});
    patternColumns = vertcat(columnParts{1:partIndex});
    pattern = sparse(patternRows,patternColumns,1,nState,nState);
    pattern = spones(pattern);

end


function validate_stack_options(opt)

    positiveScalarFields = {'nCells','tMax','outputStep','RelTol','AbsTol', ...
        'MaxStep','qMaxCcm2','jMaxAcm2','minimumVoltageV', ...
        'maximumIceVolumeFraction','maximumPoreOccupancy','cellAreaM2', ...
        'endConvectionCoefficient','bipolarPlateThicknessM', ...
        'bipolarPlateConductivity','endPlateThicknessM', ...
        'endPlateDensity','endPlateHeatCapacity','endPlateConductivity', ...
        'concentrationFactorEnd', ...
        'concentrationFactorMiddle','j0Ref','tauHyd','lambdaHydWet','nHyd', ...
        'ionomerBruggemanExponent'};
    for k = 1:numel(positiveScalarFields)
        fieldName = positiveScalarFields{k};
        value = opt.(fieldName);
        if ~isscalar(value) || ~isfinite(value) || value <= 0
            error('pemfc_stack5_simulate:InvalidPositiveOption', ...
                '参数%s必须是正有限标量。',fieldName);
        end
    end
    nonnegativeScalarFields = {'Ea','hVaporCathode','hVaporAnode', ...
        'interfaceRateFactor','lambdaHydOn','kFreeze','kMelt','gammaIce', ...
        'phaseTransitionWidthK','progressIntervalS','progressWallIntervalS'};
    for k = 1:numel(nonnegativeScalarFields)
        fieldName = nonnegativeScalarFields{k};
        value = opt.(fieldName);
        if ~isscalar(value) || ~isfinite(value) || value < 0
            error('pemfc_stack5_simulate:InvalidNonnegativeOption', ...
                '参数%s必须是非负有限标量。',fieldName);
        end
    end
    if ~isscalar(opt.maxWallTimeS) || isnan(opt.maxWallTimeS) || ...
            opt.maxWallTimeS <= 0
        error('pemfc_stack5_simulate:InvalidWallTimeLimit', ...
            'maxWallTimeS must be a positive scalar or Inf.');
    end
    if ~isscalar(opt.fHydDry) || opt.fHydDry <= 0 || opt.fHydDry > 1
        error('pemfc_stack5_simulate:InvalidHydrationArea', ...
            'fHydDry必须位于(0,1]。');
    end
    if ~isscalar(opt.useJacobianPattern) || ...
            ~ismember(double(opt.useJacobianPattern),[0 1])
        error('pemfc_stack5_simulate:InvalidJacobianPatternSwitch', ...
            'useJacobianPattern必须是逻辑标量或0/1。');
    end
    if ~isscalar(opt.directFreezeFraction) || ...
            ~isfinite(opt.directFreezeFraction) || ...
            opt.directFreezeFraction < 0 || opt.directFreezeFraction > 1
        error('pemfc_stack5_simulate:InvalidDirectFreezeFraction', ...
            'directFreezeFraction必须位于[0,1]。');
    end
    if ~isscalar(opt.freezeNucleationLambdaFraction) || ...
            ~isfinite(opt.freezeNucleationLambdaFraction) || ...
            opt.freezeNucleationLambdaFraction < 0 || ...
            opt.freezeNucleationLambdaFraction > 1
        error('pemfc_stack5_simulate:InvalidNucleationThreshold', ...
            'freezeNucleationLambdaFraction必须位于[0,1]。');
    end
    if opt.lambdaHydWet <= opt.lambdaHydOn
        error('pemfc_stack5_simulate:InvalidHydrationThreshold', ...
            'lambdaHydWet必须大于lambdaHydOn。');
    end
    if ~isstruct(opt.init) || ~isscalar(opt.init)
        error('pemfc_stack5_simulate:InvalidInitialOptions', ...
            'init必须是标量结构体。');
    end
    if ~isscalar(opt.ambientTemperatureC) || ...
            ~isfinite(opt.ambientTemperatureC) || ...
            opt.ambientTemperatureC <= -273.15
        error('pemfc_stack5_simulate:InvalidAmbientTemperature', ...
            '环境温度必须高于绝对零度。');
    end
    if opt.maximumIceVolumeFraction > 1
        error('pemfc_stack5_simulate:InvalidIceConstraint', ...
            '最大冰体积分数不能大于1。');
    end
    if opt.maximumPoreOccupancy > 1
        error('pemfc_stack5_simulate:InvalidPoreOccupancyConstraint', ...
            '最大孔隙占用率不能大于1。');
    end
    validEndPlateModels = {'lumped','separate'};
    if ~ischar(opt.endPlateModel) && ~isstring(opt.endPlateModel) || ...
            ~ismember(lower(char(opt.endPlateModel)),validEndPlateModels)
        error('pemfc_stack5_simulate:InvalidEndPlateModel', ...
            'endPlateModel只能为lumped或separate。');
    end
    validBipolarCounting = {'per_cell_pair','shared_stack'};
    if (~ischar(opt.bipolarPlateCounting) && ...
            ~isstring(opt.bipolarPlateCounting)) || ...
            ~ismember(lower(char(opt.bipolarPlateCounting)), ...
            validBipolarCounting)
        error('pemfc_stack5_simulate:InvalidBipolarPlateCounting', ...
            ['bipolarPlateCounting只能为per_cell_pair', ...
             '或shared_stack。']);
    end

end


function plot_stack_result(result)

    figure('Name','Question 2: five-cell stack cold start');
    tiledlayout(2,2);

    nexttile;
    yyaxis left;
    plot(result.t,result.currentDensityAcm2,'LineWidth',1.5);
    ylabel('j / A cm^{-2}');
    yyaxis right;
    plot(result.t,result.chargeCcm2,'LineWidth',1.5);
    ylabel('q_{use} / C cm^{-2}');
    xlabel('Time / s'); title('Current and charge'); grid on;

    nexttile;
    plot(result.t,result.TavgC,'LineWidth',1.2);
    yline(0,'k--');
    xlabel('Time / s'); ylabel('Temperature / ^\circC');
    title('Cell-average temperature'); grid on;

    nexttile;
    plot(result.t,result.Vcell,'LineWidth',1.2);
    yline(result.constraints.minimumVoltageV,'k--');
    xlabel('Time / s'); ylabel('Voltage / V');
    title('Individual-cell voltage'); grid on;

    nexttile;
    iceLines = plot(result.t,result.maxIceFillingRatioCell,'LineWidth',1.2);
    yline(1,'k--');
    xlabel('Time / s'); ylabel('Maximum local ice filling ratio');
    title('Pore ice filling ratio'); ylim([0,1.05]); grid on;

    legendText = arrayfun(@(k)sprintf('Cell %d',k),1:5, ...
        'UniformOutput',false);
    legend(iceLines,legendText,'Location','best');

end
