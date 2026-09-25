%-------------------------------------------------------------------------------
% @function: 问题2——5片串联电堆给定电荷约束下的自冷启动仿真
% @author:   PJ, DSH
% @date:     20260925
%
% 电池模型：Code_ICE_F 的单电池物理内核（水-冰、热、气体三个模块原样调用），
%   本文件只增加电堆层：五片热耦合、端板独立热节点、共同电流加载、
%   累积电荷状态、题目约束的事件判定。水/冰/气体/电压方程全部不重写。
%
% @input:  strategy->q2_current_strategy 定义的加载策略
%          initialTemperatureC->电堆初始温度 [degC]，题(1)取 -10
%          simOpt->可选设置，见下方默认值
% @output: result->启动结果、时间历程和各约束余量
%          model ->网格、状态编号、热耦合参数和RHS句柄
%
% 状态空间：
%   X = [x_1;x_2;x_3;x_4;x_5; q_use; T_endL; T_endR]
%   其中 x_k 与问题1完全一致（见 pemfc_setup_ice.m），
%   q_use 为累积电荷 [C/cm^2]，T_endL/T_endR 为两块端板温度 [K]。
%
% [Q2F-3] 相邻片导热（题目式2）：对式(2)两边除以电池面积 A，
%   得到进入一维能量方程的面热通量 q''= (k_eff/delta)*(T_k-T_{k+1})。
%   同一个 q'' 作为第 k 片的右边界散热和第 k+1 片的左边界吸热，
%   因此内部换热对电堆总能严格守恒（见 verify 字段）。
% [Q2F-4] 端板作为独立状态求解，通过"半MEA+一块双极板+半端板"的串联热阻
%   与第1/5片换热，端板外侧再与环境对流。不把端板与端部电池强制同温。
% [Q2F-5] 5片堆共6块双极板，端部片计1.5块、中间片计1块，
%   相对"每片各计两块"的基准倍率为 [0.75,0.5,0.5,0.5,0.75]。
% [Q2F-6] 附件1 的首尾电池浓差极化系数（10/1）乘在 eta_con 上，
%   由 thermal_temperature_state_ice 的 concentrationLossFactor 实现。
% [Q2F-7] 累积电荷作为状态积分：dq_use/dt = j[A/cm^2]，单位自然为 C/cm^2。
% [Q2F-8] 事件判定：五片最低温度达 0 degC 为成功候选；
%   q_use>=20 C/cm^2、最低电压<0.30 V、最大冰体积分数>=0.99、
%   j>0.5 A/cm^2、孔隙占用率>=0.98 为失败边界。
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
        error('pemfc_stack5_simulate:InvalidOptions','simOpt必须是标量结构体。');
    end
    if ~isscalar(initialTemperatureC) || ~isfinite(initialTemperatureC) || ...
            initialTemperatureC <= -273.15
        error('pemfc_stack5_simulate:InvalidInitialTemperature', ...
            '初始温度必须是高于绝对零度的有限标量。');
    end
    % 先校验策略本身，避免把参数错误拖到积分中途才暴露。
    q2_current_strategy(0,strategy);

    simOpt = set_default(simOpt,'nCells',5);
    simOpt = set_default(simOpt,'nCellLayer',[6 13 19 19 6]);
    simOpt = set_default(simOpt,'ambientTemperatureC',initialTemperatureC);
    simOpt = set_default(simOpt,'tMax',300);
    simOpt = set_default(simOpt,'outputStep',0.5);
    simOpt = set_default(simOpt,'RelTol',1e-6);
    simOpt = set_default(simOpt,'AbsTol',1e-8);
    simOpt = set_default(simOpt,'MaxStep',0.05);
    simOpt = set_default(simOpt,'verbose',true);
    simOpt = set_default(simOpt,'plot',false);
    simOpt = set_default(simOpt,'useJacobianPattern',true);
    % [Q2F-12] 求解进度打印间隔 [s]（按**仿真物理时间**计），0 = 关闭。
    % 五片全阶模型单次求解可达数十秒到数分钟，没有进度输出时无法判断
    % 是在正常推进还是卡在某个刚性区。见 stack_progress_output。
    simOpt = set_default(simOpt,'progressIntervalS',0);

    % 题目硬约束（问题2题干）
    simOpt = set_default(simOpt,'qMaxCcm2',20);       % 累积电荷上限
    simOpt = set_default(simOpt,'jMaxAcm2',0.5);      % 最大电流密度
    simOpt = set_default(simOpt,'minimumVoltageV',0.30);
    simOpt = set_default(simOpt,'maximumIceVolumeFraction',0.99);
    % 水+冰总占用率达到 98% 时先判孔隙堵塞。留 2% 残余孔隙是给
    % ode15s 的 Newton 试探步留数值裕度，避免事件根被接受前先跨入
    % 扩散系数奇异区（与 Code_Q2_PhaseFix 的处理一致）。
    simOpt = set_default(simOpt,'maximumPoreOccupancy',0.98);

    % 题干与附件1 给定的电堆热参数
    simOpt = set_default(simOpt,'cellAreaM2',25e-4);           % A=25 cm^2
    simOpt = set_default(simOpt,'endConvectionCoefficient',40); % h=40 W/(m^2 K)
    simOpt = set_default(simOpt,'bipolarPlateThicknessM',2e-3); % 双极板 2 mm
    simOpt = set_default(simOpt,'bipolarPlateConductivity',95);
    simOpt = set_default(simOpt,'endPlateThicknessM',10e-3);    % 端板 10 mm
    simOpt = set_default(simOpt,'endPlateDensity',7900);
    simOpt = set_default(simOpt,'endPlateHeatCapacity',500);
    simOpt = set_default(simOpt,'endPlateConductivity',15);
    simOpt = set_default(simOpt,'endPlateModel','separate');
    simOpt = set_default(simOpt,'bipolarPlateCounting','shared_stack');
    % 附件1：首尾电池浓差极化系数 10，中间电池 1
    simOpt = set_default(simOpt,'concentrationFactorEnd',10);
    simOpt = set_default(simOpt,'concentrationFactorMiddle',1);

    % [Q2F-9] 电化学与产冰参数直接沿用 Code_ICE_F（问题1 的标定结果），
    % 以保证"电堆层只是加载问题1的单电池模型"。
    % 需要注意：F 版的 j0Ref/mHyd/clResistanceScale 是由 -20/-25 degC 两条
    % 曲线标定的，本文在 -10 degC 使用属于温度外推，论文中必须说明。
    simOpt = set_default(simOpt,'j0Ref',0.28);
    simOpt = set_default(simOpt,'Ea',67000);
    simOpt = set_default(simOpt,'mHyd',1.70);
    simOpt = set_default(simOpt,'lamHydRef',1.0);
    simOpt = set_default(simOpt,'clResistanceScale',0.66);
    simOpt = set_default(simOpt,'tauHyd',40);
    simOpt = set_default(simOpt,'gammaIce',3.5);
    simOpt = set_default(simOpt,'directFreezeFraction',0.15);
    simOpt = set_default(simOpt,'freezeNucleationLambdaFraction',0.25);
    simOpt = set_default(simOpt,'kFreezeDirect',1);
    simOpt = set_default(simOpt,'kFreeze',0.4);
    % 问题2 的温度会升到 0 degC 以上，冰必须能融化，故 kMelt=1（不可辨识，
    % 取文献量级；问题1 全程 T<0 时 kMelt 不起作用）。
    simOpt = set_default(simOpt,'kMelt',1);
    % 问题2 会跨过冰点，冻结/融化硬开关会让 ode15s 在冰点附近反复抖振
    % （Code_Q2_PhaseFix 实测在 58 s 附近触发最小步长失败），
    % 故在冰点邻域启用 0.2 K 窄平滑层。远离冰点时与原分段式等价。
    simOpt = set_default(simOpt,'phaseTransitionWidthK',0.2);
    simOpt = set_default(simOpt,'interfaceFactor',30);
    simOpt = set_default(simOpt,'hVaporAnode',0);
    simOpt = set_default(simOpt,'hVaporCathode',0);
    simOpt = set_default(simOpt,'kmCond',0);
    simOpt = set_default(simOpt,'dTref',20);
    simOpt = set_default(simOpt,'nRet',0);
    simOpt = set_default(simOpt,'rRet',0);
    simOpt = set_default(simOpt,'fRet',1);
    simOpt = set_default(simOpt,'DwFilmRef',0);
    simOpt = set_default(simOpt,'nSorp',3);
    simOpt = set_default(simOpt,'Kcov',0);
    if ~isfield(simOpt,'init') || isempty(simOpt.init)
        simOpt.init = struct('lambdaCCL0',1.0);
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

    % 多孔层的总孔隙率（附件1 五层取值），用于题目式(9) 与孔隙占用率
    eps0Layer = [0.8,0.3916,0,0.4207,0.8];
    eps0 = eps0Layer(g.layerId(:));
    eps0 = eps0(:);
    idxPorous = g.idx_porous;


    %% 3. 电堆热耦合参数

    materialConductivity = [0.3,0.27,0.24,0.27,0.3];
    meaThermalResistance = sum(g.layerThickness./materialConductivity);
    if strcmpi(char(simOpt.bipolarPlateCounting),'shared_stack')
        % [Q2F-5] 5片堆共6块双极板：端部片1.5块、中间片1块，合计6块，
        % 而基准单片热容已含2块，故倍率为 [0.75,0.5,0.5,0.5,0.75]。
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
    stackThermal.ambientTemperatureK = simOpt.ambientTemperatureC+273.15;
    stackThermal.endConvectionCoefficient = simOpt.endConvectionCoefficient;
    stackThermal.endPlateArealCapacity = ...
        simOpt.endPlateDensity*simOpt.endPlateHeatCapacity* ...
        simOpt.endPlateThicknessM;
    stackThermal.endPlateModel = char(simOpt.endPlateModel);
    stackThermal.bipolarPlateCounting = char(simOpt.bipolarPlateCounting);
    stackThermal.bipolarPlateCapacityScale = bipolarPlateCapacityScale;
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
    stackThermal.bipolarPlateThermalResistance = bipolarPlateThermalResistance;
    stackThermal.intercellThermalResistance = intercellThermalResistance;
    stackThermal.stackEffectiveConductivity = stackEffectiveConductivity;
    stackThermal.intercellConductancePerArea = intercellConductancePerArea;

    ice = struct('kFreeze',simOpt.kFreeze,'kMelt',simOpt.kMelt, ...
        'gammaIce',simOpt.gammaIce);
    phaseOpt = struct('kmCond',simOpt.kmCond,'dTref',simOpt.dTref, ...
        'nRet',simOpt.nRet,'rRet',simOpt.rRet,'fRet',simOpt.fRet, ...
        'hVaporAnode',simOpt.hVaporAnode, ...
        'hVaporCathode',simOpt.hVaporCathode, ...
        'directFreezeFraction',simOpt.directFreezeFraction, ...
        'freezeNucleationLambdaFraction', ...
        simOpt.freezeNucleationLambdaFraction, ...
        'kFreezeDirect',simOpt.kFreezeDirect, ...
        'phaseTransitionWidthK',simOpt.phaseTransitionWidthK, ...
        'interfaceFactor',simOpt.interfaceFactor);
    clOpt = struct('clResistanceScale',simOpt.clResistanceScale);
    hydOpt = struct('mHyd',simOpt.mHyd,'lamHydRef',simOpt.lamHydRef);
    hydration = struct('tauHyd',simOpt.tauHyd);
    activation = struct('j0Ref',simOpt.j0Ref,'Ea',simOpt.Ea);
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

    rhsArgs = {ice,phaseOpt,clOpt,hydOpt,hydration,activation, ...
        stackThermal,simOpt,g,s,eps0,cellStateIndex,chargeIndex, ...
        endPlateIndex};
    rhs = @(t,x) stack_rhs_q2(t,x,strategy,rhsArgs{:});
    eventFunction = @(t,x) stack_events_q2(t,x,rhs,constraints,chargeIndex);
    odeOpt = odeset('RelTol',simOpt.RelTol,'AbsTol',simOpt.AbsTol, ...
        'MaxStep',simOpt.MaxStep,'NonNegative',nonNegativeStates, ...
        'Events',eventFunction);
    if simOpt.useJacobianPattern
        % [Q2F-10] 只声明真实存在的耦合：片内取稠密块，片间只有相邻片的
        % 导热（经由 Tavg 耦合到全部温度节点），端板只与第1/5片耦合。
        % 这不是新的物理方程，只减少 ode15s 为数值雅可比分组扰动状态时
        % 重复调用整套五片 RHS 的次数。
        odeOpt = odeset(odeOpt,'JPattern',build_stack_jacobian_pattern( ...
            s,cellStateIndex,chargeIndex,endPlateIndex,numel(x0)));
    end
    if simOpt.progressIntervalS > 0
        % [Q2F-12] 进度输出：每推进 progressIntervalS 秒**仿真时间**打印一行，
        % 同时给出实际等待时间。ode15s 的 OutputFcn 只在被接受的步之后调用，
        % 因此不会在 Newton/RHS 内部刷屏，对刚性求解的额外开销可忽略。
        outputFunction = @(t,y,flag) stack_progress_output( ...
            t,y,flag,simOpt.progressIntervalS,simOpt.tMax);
        odeOpt = odeset(odeOpt,'OutputFcn',outputFunction);
    end

    if simOpt.verbose
        fprintf('\n========== 问题2 5片电堆自冷启动 ==========\n');
        fprintf('策略=%s，T0=%.3f degC，环境=%.3f degC\n', ...
            char(strategy.type),initialTemperatureC,simOpt.ambientTemperatureC);
        if strcmpi(char(strategy.type),'constant')
            fprintf('恒定电流密度=%.6g A/cm^2\n',strategy.currentAcm2);
        elseif strcmpi(char(strategy.type),'linear')
            fprintf('升载速率=%.6g A/(cm^2 s)，平台=%.6g A/cm^2\n', ...
                strategy.rampRateAcm2s,strategy.plateauAcm2);
        else
            fprintf('阶梯电流=[%s] A/cm^2，切换时刻=[%s] s\n', ...
                num2str(strategy.levelsAcm2(:).'), ...
                num2str(strategy.switchTimesS(:).'));
        end
        fprintf(['相邻节距=%.6g m，k_eff=%.6g W/(m K)，', ...
            'G/A=%.6g W/(m^2 K)\n'],intercellDistance, ...
            stackEffectiveConductivity,intercellConductancePerArea);
        fprintf('端板模式=%s，双极板计数=%s，浓差系数(首尾/中间)=%g/%g\n', ...
            stackThermal.endPlateModel,stackThermal.bipolarPlateCounting, ...
            simOpt.concentrationFactorEnd,simOpt.concentrationFactorMiddle);
    end

    % 恒流策略可能在 t=0 就使电流或电压越界。MATLAB 事件只检测"穿过零点"，
    % 因此必须在积分前单独判定初值。
    [~,initialOut] = rhs(0,x0);
    initialEventIndex = [];
    if initialOut.jAcm2 > constraints.jMaxAcm2+max(1e-10,1e-8*constraints.jMaxAcm2)
        initialEventIndex = 5;
    elseif min(initialOut.Vcell) < constraints.minimumVoltageV-1e-10
        initialEventIndex = 3;
    elseif max(initialOut.maxIceVolumeFractionCell) >= ...
            constraints.maximumIceVolumeFraction
        initialEventIndex = 4;
    elseif max(initialOut.maxPoreOccupancyCell) >= constraints.maximumPoreOccupancy
        initialEventIndex = 6;
    end

    tic;
    if isempty(initialEventIndex)
        [tSol,xSol,tEvent,~,eventIndex] = ode15s(rhs,tOutput,x0,odeOpt);
    else
        tSol = 0;
        xSol = x0.';
        tEvent = 0;
        eventIndex = initialEventIndex;
    end
    cpuTime = toc;


    %% 5. 重算各片代数输出

    nTime = numel(tSol);
    TavgC = zeros(nTime,nCells);
    Vcell = zeros(nTime,nCells);
    minimumCellVoltageV = zeros(nTime,nCells);
    iceVolumeFractionMax = zeros(nTime,nCells);
    iceVolumeFractionCCL = zeros(nTime,nCells);
    iceSaturationCCL = zeros(nTime,nCells);
    poreOccupancyMax = zeros(nTime,nCells);
    lambdaCCL = zeros(nTime,nCells);
    endPlateTemperatureC = zeros(nTime,2);
    interfaceHeatFlux = zeros(nTime,nCells-1);
    currentDensityAcm2 = zeros(nTime,1);
    chargeState = zeros(nTime,1);
    cellOut = cell(nTime,nCells);

    for n = 1:nTime
        [~,out] = rhs(tSol(n),xSol(n,:).');
        TavgC(n,:) = (out.Tavg-273.15).';
        Vcell(n,:) = out.Vcell.';
        minimumCellVoltageV(n,:) = out.minimumVoltageCell.';
        iceVolumeFractionMax(n,:) = out.maxIceVolumeFractionCell.';
        iceVolumeFractionCCL(n,:) = out.iceVolumeFractionCCL.';
        iceSaturationCCL(n,:) = out.iceSaturationCCL.';
        poreOccupancyMax(n,:) = out.maxPoreOccupancyCell.';
        lambdaCCL(n,:) = out.lambdaCCL.';
        endPlateTemperatureC(n,:) = (out.endPlateTemperatureK-273.15).';
        interfaceHeatFlux(n,:) = out.interfaceHeatFlux.';
        currentDensityAcm2(n) = out.jAcm2;
        chargeState(n) = xSol(n,chargeIndex);
        cellOut(n,:) = out.cell;
    end


    %% 6. 按题目式(4)-(9) 统一判定启动成功与失败

    % 注意各约束的方向不同（这是最容易写错的地方）：
    %   温度：要求"五片平均温度**同时达到** 0 degC"，
    %         即逐时刻取五片最小值，再对时间取**最大**；
    %         不能用全程最低温度判定，否则初温 -10 degC 必然使判据永假。
    %   电压/电流/冰：要求**全程**不越界，故对时间取最不利值。
    %   电荷：只看末端累积值。
    cellMinimumTemperatureC = min(TavgC,[],2);      % 逐时刻五片最低温
    startupMinimumTemperatureC = max(cellMinimumTemperatureC);
    historyMinimumTemperatureC = min(TavgC,[],'all');
    historyMinimumVoltageV = min(minimumCellVoltageV,[],'all');
    historyMaximumIce = max(iceVolumeFractionMax,[],'all');
    historyMaximumPoreOccupancy = max(poreOccupancyMax,[],'all');
    historyMaximumCurrent = max(currentDensityAcm2);
    chargeUsedCcm2 = chargeState(end);

    reachedZeroC = startupMinimumTemperatureC >= 0;
    [successTimeS,successTimeIndex] = startup_time(TavgC,tSol);

    failureReasons = {};
    if chargeUsedCcm2 > constraints.qMaxCcm2+1e-9
        failureReasons{end+1} = ...
            sprintf('累积电荷 %.4g C/cm^2 超过上限 %.4g', ...
            chargeUsedCcm2,constraints.qMaxCcm2);
    end
    if historyMaximumCurrent > constraints.jMaxAcm2+1e-9
        failureReasons{end+1} = ...
            sprintf('最大电流密度 %.4g A/cm^2 超过上限 %.4g', ...
            historyMaximumCurrent,constraints.jMaxAcm2);
    end
    if historyMinimumVoltageV < constraints.minimumVoltageV
        failureReasons{end+1} = ...
            sprintf('最低单片电压 %.4g V 低于下限 %.4g', ...
            historyMinimumVoltageV,constraints.minimumVoltageV);
    end
    if historyMaximumIce >= constraints.maximumIceVolumeFraction
        failureReasons{end+1} = ...
            sprintf('最大冰体积分数 %.4g 达到阈值 %.4g', ...
            historyMaximumIce,constraints.maximumIceVolumeFraction);
    end
    if historyMaximumPoreOccupancy >= constraints.maximumPoreOccupancy
        failureReasons{end+1} = ...
            sprintf('最大孔隙占用率 %.4g 达到阈值 %.4g', ...
            historyMaximumPoreOccupancy,constraints.maximumPoreOccupancy);
    end
    if ~reachedZeroC
        failureReasons{end+1} = sprintf( ...
            '五片平均温度未能同时达到 0 degC（最高仅 %.4g degC）', ...
            startupMinimumTemperatureC);
    end
    success = isempty(failureReasons);

    % 约束余量：>0 表示该约束仍有裕量
    margin.qMaxCcm2 = constraints.qMaxCcm2-chargeUsedCcm2;
    margin.jMaxAcm2 = constraints.jMaxAcm2-historyMaximumCurrent;
    margin.minimumVoltageV = historyMinimumVoltageV-constraints.minimumVoltageV;
    margin.maximumIceVolumeFraction = ...
        constraints.maximumIceVolumeFraction-historyMaximumIce;
    margin.maximumPoreOccupancy = ...
        constraints.maximumPoreOccupancy-historyMaximumPoreOccupancy;
    margin.temperatureC = startupMinimumTemperatureC;

    result.strategy = strategy;
    result.initialTemperatureC = initialTemperatureC;
    result.ambientTemperatureC = simOpt.ambientTemperatureC;
    result.success = success;
    result.successTimeS = successTimeS;
    result.failureReasons = failureReasons;
    result.terminalEventIndex = eventIndex;
    result.terminalEventTimeS = tEvent;
    result.chargeUsedCcm2 = chargeUsedCcm2;
    result.maximumCurrentDensityAcm2 = historyMaximumCurrent;
    result.minimumCellVoltageV = historyMinimumVoltageV;
    result.maximumIceVolumeFraction = historyMaximumIce;
    result.maximumPoreOccupancy = historyMaximumPoreOccupancy;
    result.minimumCellTemperatureC = historyMinimumTemperatureC;
    result.startupMinimumTemperatureC = startupMinimumTemperatureC;
    result.margin = margin;

    result.t = tSol;
    result.x = xSol;
    result.TavgC = TavgC;
    result.Vcell = Vcell;
    result.minimumVoltagePerCell = minimumCellVoltageV;
    result.iceVolumeFractionMax = iceVolumeFractionMax;
    result.iceVolumeFractionCCL = iceVolumeFractionCCL;
    result.iceSaturationCCL = iceSaturationCCL;
    result.poreOccupancyMax = poreOccupancyMax;
    result.lambdaCCL = lambdaCCL;
    result.endPlateTemperatureC = endPlateTemperatureC;
    result.interfaceHeatFlux = interfaceHeatFlux;
    result.currentDensityAcm2 = currentDensityAcm2;
    result.chargeStateCcm2 = chargeState;
    result.cellOut = cellOut;
    result.cpuTime = cpuTime;
    result.simOpt = simOpt;

    model.g = g;
    model.s = s;
    model.x0 = x0;
    model.init = init;
    model.ice = ice;
    model.phaseOpt = phaseOpt;
    model.clOpt = clOpt;
    model.hydOpt = hydOpt;
    model.hydration = hydration;
    model.activation = activation;
    model.stackThermal = stackThermal;
    model.constraints = constraints;
    model.cellStateIndex = cellStateIndex;
    model.chargeIndex = chargeIndex;
    model.endPlateIndex = endPlateIndex;
    model.eps0 = eps0;
    model.rhs = rhs;

    if simOpt.verbose
        if success
            fprintf(['启动成功：t_s=%.3f s，累积电荷=%.4f C/cm^2，', ...
                '最低电压=%.5f V，最大冰体积分数=%.4g\n'], ...
                successTimeS,chargeUsedCcm2,historyMinimumVoltageV, ...
                historyMaximumIce);
        else
            fprintf('启动失败：\n');
            for k = 1:numel(failureReasons)
                fprintf('   - %s\n',failureReasons{k});
            end
        end
        fprintf('耗时%.2f s，积分步数=%d\n',cpuTime,nTime);
    end

    if simOpt.plot
        plot_stack_result(result);
    end

end


%% ========================================================================
% 局部函数：五片电堆RHS
% ========================================================================

function [dx,out] = stack_rhs_q2(t,x,strategy,ice,phaseOpt,clOpt,hydOpt, ...
    hydration,activation,stackThermal,simOpt,g,s,eps0,cellStateIndex, ...
    chargeIndex,endPlateIndex)

    nCells = numel(cellStateIndex);
    [j,jAcm2] = q2_current_strategy(t,strategy);
    if ~isscalar(j) || ~isfinite(j) || j < 0
        error('pemfc_stack5_simulate:InvalidCurrent', ...
            '加载策略必须返回非负有限标量电流。');
    end

    % 各片体积平均温度
    Tavg = zeros(nCells,1);
    for k = 1:nCells
        Tcell = x(cellStateIndex{k}(s.idx.T));
        Tavg(k) = sum(Tcell.*g.dx)/sum(g.dx);
    end

    % [Q2F-3] 相邻片导热面热通量（题目式2 除以面积 A）
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
            leftOutflux = leftBoundaryOutflux;
        else
            leftOutflux = -interfaceHeatFlux(k-1);
        end
        if k == nCells
            rightOutflux = rightBoundaryOutflux;
        else
            rightOutflux = interfaceHeatFlux(k);
        end

        isEndCell = k == 1 || k == nCells;
        stackBoundary.mode = 'prescribed_flux';
        stackBoundary.leftOutflux = leftOutflux;
        stackBoundary.rightOutflux = rightOutflux;
        if strcmpi(stackThermal.endPlateModel,'lumped')
            stackBoundary.extraArealCapacity = ...
                double(isEndCell)*stackThermal.endPlateArealCapacity;
        else
            stackBoundary.extraArealCapacity = 0;
        end
        stackBoundary.bipolarPlateCapacityScale = ...
            stackThermal.bipolarPlateCapacityScale(k);
        if isEndCell
            stackBoundary.concentrationLossFactor = ...
                simOpt.concentrationFactorEnd;
        else
            stackBoundary.concentrationLossFactor = ...
                simOpt.concentrationFactorMiddle;
        end

        [dxCell,cellOut{k}] = cell_rhs_q2(t,x(cellStateIndex{k}),j, ...
            stackThermal.ambientTemperatureK,stackBoundary,ice,phaseOpt, ...
            clOpt,hydOpt,hydration,activation,g,s,eps0);
        dx(cellStateIndex{k}) = dxCell;
    end

    % [Q2F-7] A/cm^2 = C/(s cm^2)
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
        out.minimumVoltageCell = zeros(nCells,1);
        out.maxIceVolumeFractionCell = zeros(nCells,1);
        out.iceVolumeFractionCCL = zeros(nCells,1);
        out.iceSaturationCCL = zeros(nCells,1);
        out.maxPoreOccupancyCell = zeros(nCells,1);
        out.lambdaCCL = zeros(nCells,1);
        out.cell = cellOut;
        for k = 1:nCells
            water = cellOut{k}.water;
            idxPorous = g.idx_porous;
            out.Vcell(k) = cellOut{k}.Vcell;
            out.minimumVoltageCell(k) = cellOut{k}.Vcell;
            out.maxIceVolumeFractionCell(k) = max(water.eps_i(idxPorous));
            out.iceVolumeFractionCCL(k) = sum( ...
                water.eps_i(g.idx_cCL).*g.dx(g.idx_cCL))/g.layerThickness(4);
            out.iceSaturationCCL(k) = max(water.sIce(g.idx_cCL));
            out.maxPoreOccupancyCell(k) = max( ...
                (water.ml(idxPorous)+water.eps_i(idxPorous))./ ...
                eps0(idxPorous));
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
% 局部函数：单片RHS，物理内核只调用水-冰、热、气体三个模块
% ========================================================================

function [dx,out] = cell_rhs_q2(t,x,j,Tamb,stackBoundary,ice,phaseOpt, ...
    clOpt,hydOpt,hydration,activation,g,s,eps0)

    T = x(s.idx.T);
    cH2State = x(s.idx.H2);
    cO2State = x(s.idx.O2);
    mwState = x(s.idx.mw);
    miState = x(s.idx.mi);
    lambdaCCLState = x(s.idx.lambdaCCL);

    cH2 = max(cH2State,0);
    cO2 = max(cO2State,0);
    mw = max(mwState,0);
    mi = max(miState,0);
    lambdaCCL = max(lambdaCCLState,0);

    [dmw,dmi,dlambdaCCL,water] = water_ice_state_ice( ...
        T,mw,mi,lambdaCCL,j,ice.kFreeze,ice.kMelt, ...
        hydration.tauHyd,phaseOpt,g,s);
    [dT,thermal] = thermal_temperature_state_ice( ...
        T,cH2,cO2,water,j,Tamb,zeros(g.N,1),g,s,ice.gammaIce, ...
        activation.j0Ref,activation.Ea,clOpt,hydOpt,stackBoundary);
    [dcH2,dcO2,gas] = gas_transport_state_ice( ...
        T,cH2,cO2,mw,dmw,dT,water,thermal,j,g,s);

    dx = zeros(s.Nx,1);
    dx(s.idx.T) = dT;
    dx(s.idx.H2) = dcH2;
    dx(s.idx.O2) = dcO2;
    dx(s.idx.mw) = dmw;
    dx(s.idx.mi) = dmi;
    dx(s.idx.lambdaCCL) = dlambdaCCL;

    % [Q2F-11] 状态导数必须是有限值。若某分量出现 NaN/Inf，ode15s 会把
    % 它带进下一次试探，最后在 water_ice_state_ice 里报成"水扩散系数非法"，
    % 完全掩盖真正原因（问题2 中端部浓差系数=10，更容易把 j 推过 j_lim）。
    if any(~isfinite(dx))
        componentNames = {'T','H2','O2','mw','mi','lambdaCCL'};
        componentIndex = {s.idx.T,s.idx.H2,s.idx.O2,s.idx.mw, ...
            s.idx.mi,s.idx.lambdaCCL};
        badComponents = '';
        for componentK = 1:numel(componentNames)
            if any(~isfinite(dx(componentIndex{componentK})))
                badComponents = [badComponents ' ' componentNames{componentK}]; %#ok<AGROW>
            end
        end
        error('pemfc_stack5_simulate:NonFiniteDerivative', ...
            ['t=%.6g s 处状态导数出现 NaN/Inf，分量：%s。\n' ...
            'j=%.6g A/m^2, Vcell=%.6g, Erev=%.6g, etaAct=%.6g, ' ...
            'etaOhm=%.6g, etaCon=%.6g, jLim=%.6g, 浓差倍率=%g。'], ...
            t,badComponents,j,thermal.Vcell,thermal.voltage.Erev, ...
            thermal.voltage.etaAct,thermal.voltage.etaOhm, ...
            thermal.voltage.etaCon,thermal.jLim, ...
            stackBoundary.concentrationLossFactor);
    end

    if nargout > 1
        out.t = t;
        out.Vcell = thermal.Vcell;
        out.Tavg = thermal.Tavg;
        out.lambdaCCL = water.lambdaCCL;
        out.water = water;
        out.thermal = thermal;
        out.gas = gas;
        out.jLim = thermal.jLim;
    end

end


%% ========================================================================
% 局部函数：启动成功与五个失败边界的事件函数
% ========================================================================

function [value,isTerminal,direction] = stack_events_q2( ...
    t,x,rhs,constraints,chargeIndex)

    [~,out] = rhs(t,x);
    minimumTemperatureC = min(out.Tavg)-273.15;
    chargeMargin = constraints.qMaxCcm2-x(chargeIndex);
    voltageMargin = min(out.minimumVoltageCell)-constraints.minimumVoltageV;
    iceMargin = constraints.maximumIceVolumeFraction- ...
        max(out.maxIceVolumeFractionCell);
    poreOccupancyMargin = constraints.maximumPoreOccupancy- ...
        max(out.maxPoreOccupancyCell);
    % j=jmax 本身是允许的；事件面略高于上限，避免平台刚好取 jmax 时误判。
    currentEventLimit = constraints.jMaxAcm2+ ...
        max(1e-10,1e-8*constraints.jMaxAcm2);
    currentMargin = currentEventLimit-out.jAcm2;

    value = [minimumTemperatureC;chargeMargin;voltageMargin;iceMargin; ...
        currentMargin;poreOccupancyMargin];
    isTerminal = [1;1;1;1;1;1];
    direction = [1;-1;-1;-1;-1;-1];

end


%% ========================================================================
% 辅助函数
% ========================================================================

function [successTimeS,successTimeIndex] = startup_time(TavgC,tSol)

    % 启动成功时刻 = 五片平均温度同时达到 0 degC 的最早时刻。
    reached = all(TavgC >= 0,2);
    successTimeIndex = find(reached,1,'first');
    if isempty(successTimeIndex)
        successTimeS = NaN;
        successTimeIndex = NaN;
    else
        successTimeS = tSol(successTimeIndex);
    end

end


function output = set_default(input,fieldName,defaultValue)

    output = input;
    if ~isfield(output,fieldName) || isempty(output.(fieldName))
        output.(fieldName) = defaultValue;
    end

end


function pattern = build_stack_jacobian_pattern( ...
    s,cellStateIndex,chargeIndex,endPlateIndex,nState)

    % [Q2F-10] 保守稀疏结构：
    %   片内         -> 稠密块（水/冰/热/气体在片内全耦合）
    %   相邻片       -> 只有温度通过 Tavg 耦合，但 Tavg 用到全部温度节点，
    %                   故取"两片温度节点"的稠密块作为上界
    %   端板         -> 只与第1/5片的全部温度节点耦合
    %   电荷状态     -> 只依赖时间，无状态耦合
    nCells = numel(cellStateIndex);
    rowIndex = [];
    columnIndex = [];
    for k = 1:nCells
        rows = cellStateIndex{k};
        rowIndex = [rowIndex;repmat(rows(:),numel(rows),1)]; %#ok<AGROW>
        columnIndex = [columnIndex;reshape(repmat(rows(:).', ...
            numel(rows),1),[],1)]; %#ok<AGROW>
    end
    for k = 1:nCells-1
        rows = [cellStateIndex{k}(s.idx.T);cellStateIndex{k+1}(s.idx.T)];
        rowIndex = [rowIndex;repmat(rows(:),numel(rows),1)]; %#ok<AGROW>
        columnIndex = [columnIndex;reshape(repmat(rows(:).', ...
            numel(rows),1),[],1)]; %#ok<AGROW>
    end
    if ~isempty(endPlateIndex)
        for side = 1:2
            if side == 1
                temperatureRows = cellStateIndex{1}(s.idx.T);
            else
                temperatureRows = cellStateIndex{end}(s.idx.T);
            end
            endPlateRow = endPlateIndex(side);
            rowIndex = [rowIndex;repmat(endPlateRow, ...
                numel(temperatureRows),1)]; %#ok<AGROW>
            columnIndex = [columnIndex;temperatureRows(:)];
            rowIndex = [rowIndex;temperatureRows(:)];
            columnIndex = [columnIndex;repmat(endPlateRow, ...
                numel(temperatureRows),1)]; %#ok<AGROW>
        end
    end
    rowIndex = [rowIndex;chargeIndex];
    columnIndex = [columnIndex;chargeIndex];
    pattern = sparse(rowIndex,columnIndex,true,nState,nState);

end


function validate_stack_options(simOpt)

    positiveScalarFields = {'tMax','outputStep','RelTol','AbsTol','MaxStep', ...
        'qMaxCcm2','jMaxAcm2','minimumVoltageV','cellAreaM2', ...
        'endConvectionCoefficient','bipolarPlateThicknessM', ...
        'bipolarPlateConductivity','endPlateThicknessM','endPlateDensity', ...
        'endPlateHeatCapacity','endPlateConductivity','j0Ref','Ea','tauHyd', ...
        'gammaIce','kFreeze','mHyd','lamHydRef','clResistanceScale'};
    for k = 1:numel(positiveScalarFields)
        fieldName = positiveScalarFields{k};
        value = simOpt.(fieldName);
        if ~isscalar(value) || ~isfinite(value) || value <= 0
            error('pemfc_stack5_simulate:InvalidOption', ...
                '选项%s必须是正有限标量。',fieldName);
        end
    end
    if ~isscalar(simOpt.nCells) || simOpt.nCells ~= round(simOpt.nCells) || ...
            simOpt.nCells < 2
        error('pemfc_stack5_simulate:InvalidCellCount', ...
            'nCells必须是不小于2的整数。');
    end
    if ~isnumeric(simOpt.nCellLayer) || numel(simOpt.nCellLayer) ~= 5 || ...
            any(simOpt.nCellLayer <= 0) || any(mod(simOpt.nCellLayer,1) ~= 0)
        error('pemfc_stack5_simulate:InvalidGrid','nCellLayer必须包含5个正整数。');
    end
    if simOpt.maximumIceVolumeFraction <= 0 || ...
            simOpt.maximumIceVolumeFraction > 1
        error('pemfc_stack5_simulate:InvalidIceLimit', ...
            'maximumIceVolumeFraction必须位于(0,1]。');
    end
    if simOpt.maximumPoreOccupancy <= 0 || simOpt.maximumPoreOccupancy > 1
        error('pemfc_stack5_simulate:InvalidPoreLimit', ...
            'maximumPoreOccupancy必须位于(0,1]。');
    end
    % progressIntervalS=0 表示关闭，因此只要求非负。
    if ~isscalar(simOpt.progressIntervalS) || ...
            ~isfinite(simOpt.progressIntervalS) || simOpt.progressIntervalS < 0
        error('pemfc_stack5_simulate:InvalidProgressInterval', ...
            'progressIntervalS必须是非负有限标量（0 表示关闭）。');
    end

end


function status = stack_progress_output(t,~,flag,intervalS,tMax)

    % [Q2F-12] 只按**仿真物理时间**间隔打印，不在每个内部步长刷屏；
    % 同时给出实际等待时间，便于判断是正常推进还是卡在刚性区。
    % 结构沿用 Code_Q2_PhaseFix 的 stack_progress_output，
    % 并合并了其单电池版本 ice_progress_output 的"实际等待"信息。
    persistent nextPrintTime wallClock lastSimulationTime
    status = 0;

    if strcmp(flag,'init')
        wallClock = tic;
        lastSimulationTime = t(1);
        nextPrintTime = t(1)+intervalS;
        fprintf(['[电堆] 求解进度：%8.3f / %8.3f s，', ...
            '实际等待 %8.1f s\n'],t(1),tMax,0);
    elseif isempty(flag)
        currentTime = t(end);
        lastSimulationTime = currentTime;
        if isempty(nextPrintTime)
            nextPrintTime = currentTime+intervalS;
        end
        if currentTime+10*eps(max(abs(currentTime),1)) >= nextPrintTime
            fprintf(['[电堆] 求解进度：%8.3f / %8.3f s，', ...
                '实际等待 %8.1f s\n'],currentTime,tMax,toc(wallClock));
            nextPrintTime = nextPrintTime+intervalS;
            while nextPrintTime <= ...
                    currentTime+10*eps(max(abs(currentTime),1))
                nextPrintTime = nextPrintTime+intervalS;
            end
        end
    elseif strcmp(flag,'done')
        if ~isempty(wallClock)
            fprintf(['[电堆] 积分结束：%8.3f / %8.3f s，', ...
                '实际等待 %8.1f s\n'],lastSimulationTime,tMax, ...
                toc(wallClock));
        end
        nextPrintTime = [];
        wallClock = [];
        lastSimulationTime = [];
    end

end


function plot_stack_result(result)

    figure('Name','PEMFC 5-cell stack cold start');
    tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

    nexttile;
    plot(result.t,result.currentDensityAcm2,'LineWidth',1.5);
    xlabel('Time / s'); ylabel('j / (A cm^{-2})');
    title('Current loading'); grid on;

    nexttile;
    plot(result.t,result.chargeStateCcm2,'LineWidth',1.5);
    xlabel('Time / s'); ylabel('q_{use} / (C cm^{-2})');
    title('Accumulated charge'); grid on;

    nexttile;
    plot(result.t,result.TavgC,'LineWidth',1.2); hold on;
    yline(0,'k--','0 ^\circC');
    xlabel('Time / s'); ylabel('Cell average temperature / ^\circC');
    title('Five cell temperatures'); grid on;

    nexttile;
    plot(result.t,result.Vcell,'LineWidth',1.2); hold on;
    yline(0.30,'k--','0.30 V');
    xlabel('Time / s'); ylabel('Cell voltage / V');
    title('Five cell voltages'); grid on;

    nexttile;
    plot(result.t,result.iceVolumeFractionMax,'LineWidth',1.2);
    xlabel('Time / s'); ylabel('\epsilon_{ice} / -');
    title('Maximum local ice volume fraction (eq.9)'); grid on;

    nexttile;
    plot(result.t,result.endPlateTemperatureC,'LineWidth',1.5);
    xlabel('Time / s'); ylabel('End plate temperature / ^\circC');
    title('End plates'); grid on;

    sgtitle(sprintf('%s strategy, T_0 = %g ^\\circC, success = %d', ...
        char(result.strategy.type),result.initialTemperatureC, ...
        result.success));

end
