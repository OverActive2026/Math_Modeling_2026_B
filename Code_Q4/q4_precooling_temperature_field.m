%-------------------------------------------------------------------------------
% @function: 计算五片PEMFC电堆在低温环境中的预冷温度场
% @author:   PJ, GPT
% @date:     20260924
% @input:    coolingDurationMin->预冷时间 [min]，非负有限标量
%            simOpt->网格、热物性、求解、保存和绘图选项
% @output:   result->完整预冷历程与可直接传给Q4冷启动模型的initialState
%            model->网格、热参数和状态编号
%
% [Q4-PRE-1] 预冷阶段j=0、qaux=0，只求解热传导和端板对流散热。
% [Q4-PRE-2] 与问题3保持同一热模型：片内有限体积温度场、片间平均温度
%             导热、独立端板温度节点和shared_stack双极板热容口径。
% [Q4-PRE-3] initialState保存每片每个网格及两端板温度，供后续优化复用；
%             后续模型会校验网格和热模型元数据，避免错接初值。
%-------------------------------------------------------------------------------

function [result,model] = q4_precooling_temperature_field( ...
    coolingDurationMin,simOpt)

    %% 1. 输入与默认参数

    if nargin < 1 || isempty(coolingDurationMin)
        coolingDurationMin = 20;
    end
    if nargin < 2 || isempty(simOpt)
        simOpt = struct();
    end
    if ~isstruct(simOpt) || ~isscalar(simOpt)
        error('q4_precooling_temperature_field:InvalidOptions', ...
            'simOpt必须是标量结构体。');
    end
    if ~isscalar(coolingDurationMin) || ...
            ~isfinite(coolingDurationMin) || coolingDurationMin < 0
        error('q4_precooling_temperature_field:InvalidDuration', ...
            '预冷时间必须是非负有限标量，单位为min。');
    end

    simOpt = set_default(simOpt,'nCells',5);
    simOpt = set_default(simOpt,'nCellLayer',[6 13 19 19 6]);
    simOpt = set_default(simOpt,'initialTemperatureC',25);
    simOpt = set_default(simOpt,'ambientTemperatureC',-30);
    simOpt = set_default(simOpt,'outputStepS',5);
    simOpt = set_default(simOpt,'RelTol',1e-7);
    simOpt = set_default(simOpt,'AbsTol',1e-9);
    simOpt = set_default(simOpt,'MaxStep',2);
    simOpt = set_default(simOpt,'cellAreaM2',25e-4);
    simOpt = set_default(simOpt,'endConvectionCoefficient',40);
    simOpt = set_default(simOpt,'bipolarPlateThicknessM',2e-3);
    simOpt = set_default(simOpt,'bipolarPlateDensity',1980);
    simOpt = set_default(simOpt,'bipolarPlateHeatCapacity',766);
    simOpt = set_default(simOpt,'bipolarPlateConductivity',95);
    simOpt = set_default(simOpt,'endPlateThicknessM',10e-3);
    simOpt = set_default(simOpt,'endPlateDensity',7900);
    simOpt = set_default(simOpt,'endPlateHeatCapacity',500);
    simOpt = set_default(simOpt,'endPlateConductivity',15);
    simOpt = set_default(simOpt,'endPlateModel','separate');
    simOpt = set_default(simOpt,'bipolarPlateCounting','shared_stack');
    simOpt = set_default(simOpt,'plot',false);
    simOpt = set_default(simOpt,'verbose',true);
    simOpt = set_default(simOpt,'saveFile','');

    validate_options(simOpt);
    if simOpt.nCells ~= 5
        error('q4_precooling_temperature_field:FiveCellsRequired', ...
            '问题4固定为5片串联单电池。');
    end
    if ~strcmpi(char(simOpt.endPlateModel),'separate')
        error('q4_precooling_temperature_field:SeparateEndPlatesRequired', ...
            '预冷温度场必须使用separate独立端板节点，以便完整传递初值。');
    end

    coolingDurationS = 60*coolingDurationMin;
    initialTemperatureK = simOpt.initialTemperatureC+273.15;
    ambientTemperatureK = simOpt.ambientTemperatureC+273.15;


    %% 2. 建立与问题3一致的网格和热参数

    gridInit = struct('T0',initialTemperatureK,'lambdaCCL0',2.0);
    [g,s] = pemfc_setup_ice(simOpt.nCellLayer,gridInit);
    nCells = simOpt.nCells;
    nCellTemperatureState = g.N;
    nTemperatureState = nCells*nCellTemperatureState;
    cellTemperatureIndex = cell(nCells,1);
    for k = 1:nCells
        cellTemperatureIndex{k} = ...
            (k-1)*nCellTemperatureState+(1:nCellTemperatureState);
    end
    endPlateIndex = nTemperatureState+(1:2);

    thermal = build_thermal_parameters(g,simOpt);
    x0 = initialTemperatureK*ones(nTemperatureState+2,1);


    %% 3. 求解纯热预冷过程

    if coolingDurationS == 0
        tSol = 0;
        xSol = x0.';
        cpuTime = 0;
    else
        tOutput = (0:simOpt.outputStepS:coolingDurationS).';
        if tOutput(end) < coolingDurationS
            tOutput(end+1,1) = coolingDurationS;
        end
        rhs = @(t,x) precooling_rhs(t,x,g,thermal, ...
            cellTemperatureIndex,endPlateIndex);
        odeOpt = odeset('RelTol',simOpt.RelTol,'AbsTol',simOpt.AbsTol, ...
            'MaxStep',simOpt.MaxStep, ...
            'NonNegative',(1:numel(x0)).');

        if simOpt.verbose
            fprintf(['\n开始计算问题4预冷温度场：25 degC -> ', ...
                '环境 %.3f degC，预冷 %.3f min\n'], ...
                simOpt.ambientTemperatureC,coolingDurationMin);
        end
        tic;
        [tSol,xSol] = ode15s(rhs,tOutput,x0,odeOpt);
        cpuTime = toc;
    end


    %% 4. 整理温度场、能量校核和标准初值

    nTime = numel(tSol);
    temperatureFieldK = zeros(nTime,g.N,nCells);
    for k = 1:nCells
        temperatureFieldK(:,:,k) = xSol(:,cellTemperatureIndex{k});
    end
    temperatureFieldC = temperatureFieldK-273.15;
    cellAverageTemperatureK = zeros(nTime,nCells);
    spatialWeight = g.dx/sum(g.dx);
    for k = 1:nCells
        cellAverageTemperatureK(:,k) = ...
            temperatureFieldK(:,:,k)*spatialWeight;
    end
    cellAverageTemperatureC = cellAverageTemperatureK-273.15;
    endPlateTemperatureK = xSol(:,endPlateIndex);
    endPlateTemperatureC = endPlateTemperatureK-273.15;

    leftConvectionHeatFluxWm2 = ...
        thermal.endPlateExternalConductancePerArea.* ...
        (endPlateTemperatureK(:,1)-ambientTemperatureK);
    rightConvectionHeatFluxWm2 = ...
        thermal.endPlateExternalConductancePerArea.* ...
        (endPlateTemperatureK(:,2)-ambientTemperatureK);
    externalHeatLossW = simOpt.cellAreaM2.* ...
        (leftConvectionHeatFluxWm2+rightConvectionHeatFluxWm2);
    heatRemovedByConvectionJ = trapz(tSol,externalHeatLossW);

    cellSensibleEnergyDecreaseJ = simOpt.cellAreaM2*sum( ...
        thermal.cellArealCapacityJm2K(:).'.* ...
        (initialTemperatureK-cellAverageTemperatureK(end,:)));
    endPlateSensibleEnergyDecreaseJ = simOpt.cellAreaM2* ...
        thermal.endPlateArealCapacityJm2K*sum( ...
        initialTemperatureK-endPlateTemperatureK(end,:));
    sensibleEnergyDecreaseJ = cellSensibleEnergyDecreaseJ+ ...
        endPlateSensibleEnergyDecreaseJ;
    energyBalanceErrorJ = sensibleEnergyDecreaseJ- ...
        heatRemovedByConvectionJ;
    energyBalanceRelativeError = abs(energyBalanceErrorJ)/ ...
        max(abs(sensibleEnergyDecreaseJ),1);

    finalTemperatureCellK = reshape( ...
        xSol(end,1:nTemperatureState),g.N,nCells);
    finalEndPlateTemperatureK = xSol(end,endPlateIndex).';
    initialState = struct();
    initialState.schemaVersion = 1;
    initialState.kind = 'q4_precooling_temperature_state';
    initialState.temperatureCellK = finalTemperatureCellK;
    initialState.temperatureCellC = finalTemperatureCellK-273.15;
    initialState.endPlateTemperatureK = finalEndPlateTemperatureK;
    initialState.endPlateTemperatureC = ...
        finalEndPlateTemperatureK-273.15;
    initialState.cellAverageTemperatureK = ...
        cellAverageTemperatureK(end,:).';
    initialState.cellAverageTemperatureC = ...
        cellAverageTemperatureC(end,:).';
    initialState.nCells = nCells;
    initialState.nCellLayer = g.nCellLayer;
    initialState.layerNames = g.layerNames;
    initialState.xWithinCellM = g.x;
    initialState.dxM = g.dx;
    initialState.initialTemperatureC = simOpt.initialTemperatureC;
    initialState.ambientTemperatureC = simOpt.ambientTemperatureC;
    initialState.coolingDurationS = coolingDurationS;
    initialState.coolingDurationMin = coolingDurationMin;
    initialState.endPlateModel = char(simOpt.endPlateModel);
    initialState.bipolarPlateCounting = ...
        char(simOpt.bipolarPlateCounting);
    initialState.thermalModelId = ...
        'Q4_PRECOOL_Q3_THERMAL_V1_SEPARATE_END_PLATES';

    % 片内网格按aGDL->cGDL固定编号，左右端片受到相反方向的边界
    % 热通量，因此逐网格温度轮廓并非同编号相等。电堆左右对称性应
    % 由成对单片的体积平均温度及两端板温度检验。
    stackSymmetryErrorK = max(abs( ...
        initialState.cellAverageTemperatureK(1)- ...
        initialState.cellAverageTemperatureK(end)),abs( ...
        initialState.cellAverageTemperatureK(2)- ...
        initialState.cellAverageTemperatureK(end-1)));
    stackSymmetryErrorK = max(stackSymmetryErrorK,abs( ...
        finalEndPlateTemperatureK(1)-finalEndPlateTemperatureK(2)));

    result.tS = tSol;
    result.coolingDurationS = coolingDurationS;
    result.coolingDurationMin = coolingDurationMin;
    result.temperatureFieldK = temperatureFieldK;
    result.temperatureFieldC = temperatureFieldC;
    result.cellAverageTemperatureK = cellAverageTemperatureK;
    result.cellAverageTemperatureC = cellAverageTemperatureC;
    result.endPlateTemperatureK = endPlateTemperatureK;
    result.endPlateTemperatureC = endPlateTemperatureC;
    result.leftConvectionHeatFluxWm2 = leftConvectionHeatFluxWm2;
    result.rightConvectionHeatFluxWm2 = rightConvectionHeatFluxWm2;
    result.heatRemovedByConvectionJ = heatRemovedByConvectionJ;
    result.sensibleEnergyDecreaseJ = sensibleEnergyDecreaseJ;
    result.energyBalanceErrorJ = energyBalanceErrorJ;
    result.energyBalanceRelativeError = energyBalanceRelativeError;
    result.maximumCellAverageSpreadC = max( ...
        max(cellAverageTemperatureC,[],2)- ...
        min(cellAverageTemperatureC,[],2));
    result.finalCellAverageSpreadC = ...
        max(cellAverageTemperatureC(end,:))- ...
        min(cellAverageTemperatureC(end,:));
    result.finalMinimumTemperatureC = min([ ...
        finalTemperatureCellK(:);finalEndPlateTemperatureK])-273.15;
    result.finalMaximumTemperatureC = max([ ...
        finalTemperatureCellK(:);finalEndPlateTemperatureK])-273.15;
    result.symmetryErrorK = stackSymmetryErrorK;
    result.initialState = initialState;
    result.options = simOpt;
    result.cpuTime = cpuTime;

    model.g = g;
    model.s = s;
    model.thermal = thermal;
    model.cellTemperatureIndex = cellTemperatureIndex;
    model.endPlateIndex = endPlateIndex;
    model.x0 = x0;

    if simOpt.verbose
        fprintf('各片预冷终态平均温度 [degC] = %s\n', ...
            mat2str(cellAverageTemperatureC(end,:),6));
        fprintf('两端板预冷终态温度 [degC] = %s\n', ...
            mat2str(endPlateTemperatureC(end,:),6));
        fprintf(['终态片间温差=%.6g degC，温度场对称误差=', ...
            '%.3e K，能量相对误差=%.3e\n'], ...
            result.finalCellAverageSpreadC,result.symmetryErrorK, ...
            result.energyBalanceRelativeError);
    end

    if ~isempty(simOpt.saveFile)
        saveFile = char(simOpt.saveFile);
        saveFolder = fileparts(saveFile);
        if ~isempty(saveFolder) && ~exist(saveFolder,'dir')
            mkdir(saveFolder);
        end
        saveData = struct('initialState',initialState, ...
            'precoolingResult',result);
        save(saveFile,'-struct','saveData','-v7.3');
        if simOpt.verbose
            fprintf('已保存可复用初值：%s\n',saveFile);
        end
    end

    if simOpt.plot
        plot_precooling_result(result,g);
    end

end


%% ========================================================================
% 预冷阶段温度状态方程
% ========================================================================

function dx = precooling_rhs(~,x,g,thermal,cellIndex,endPlateIndex)

    nCells = numel(cellIndex);
    Tavg = zeros(nCells,1);
    spatialWeight = g.dx/sum(g.dx);
    for k = 1:nCells
        Tavg(k) = x(cellIndex{k}).'*spatialWeight;
    end

    interfaceHeatFlux = thermal.intercellConductancePerArea.* ...
        (Tavg(1:end-1)-Tavg(2:end));
    endPlateTemperatureK = x(endPlateIndex);
    leftCellToEndPlateFlux = ...
        thermal.cellToEndPlateConductancePerArea* ...
        (Tavg(1)-endPlateTemperatureK(1));
    rightCellToEndPlateFlux = ...
        thermal.cellToEndPlateConductancePerArea* ...
        (Tavg(end)-endPlateTemperatureK(2));
    leftConvectionOutflux = ...
        thermal.endPlateExternalConductancePerArea* ...
        (endPlateTemperatureK(1)-thermal.ambientTemperatureK);
    rightConvectionOutflux = ...
        thermal.endPlateExternalConductancePerArea* ...
        (endPlateTemperatureK(2)-thermal.ambientTemperatureK);

    dx = zeros(size(x));
    for k = 1:nCells
        if k == 1
            leftFlux = -leftCellToEndPlateFlux;
        else
            leftFlux = interfaceHeatFlux(k-1);
        end
        if k == nCells
            rightFlux = rightCellToEndPlateFlux;
        else
            rightFlux = interfaceHeatFlux(k);
        end

        T = x(cellIndex{k});
        qHeatFace = zeros(g.N+1,1);
        qHeatFace(1) = leftFlux;
        qHeatFace(end) = rightFlux;
        for f = 2:g.N
            dL = g.xf(f)-g.x(f-1);
            dR = g.x(f)-g.xf(f);
            qHeatFace(f) = -(T(f)-T(f-1))/( ...
                dL/thermal.cellEffectiveConductivityWmK(k)+ ...
                dR/thermal.cellEffectiveConductivityWmK(k));
        end
        qCond = (qHeatFace(1:end-1)-qHeatFace(2:end))./g.dx;
        dx(cellIndex{k}) = qCond./thermal.cellVolumetricCapacityJm3K(k);
    end

    dx(endPlateIndex(1)) = (leftCellToEndPlateFlux- ...
        leftConvectionOutflux)/thermal.endPlateArealCapacityJm2K;
    dx(endPlateIndex(2)) = (rightCellToEndPlateFlux- ...
        rightConvectionOutflux)/thermal.endPlateArealCapacityJm2K;

end


function thermal = build_thermal_parameters(g,opt)

    materialRho = [185,970,2150,970,185].';
    materialCp = [545,240,1050,240,545].';
    materialConductivity = [0.3,0.27,0.24,0.27,0.3].';
    meaArealCapacity = sum( ...
        materialRho.*materialCp.*g.layerThickness(:));
    meaThermalResistance = sum( ...
        g.layerThickness(:)./materialConductivity);
    meaEffectiveConductivity = g.Ltotal/meaThermalResistance;

    bipolarPlateArealCapacityBase = ...
        opt.bipolarPlateDensity*opt.bipolarPlateHeatCapacity* ...
        (2*opt.bipolarPlateThicknessM);
    if strcmpi(char(opt.bipolarPlateCounting),'shared_stack')
        bipolarPlateCapacityScale = [0.75;0.5;0.5;0.5;0.75];
        intercellBipolarPlateCount = 1;
    else
        bipolarPlateCapacityScale = ones(opt.nCells,1);
        intercellBipolarPlateCount = 2;
    end
    cellArealCapacity = meaArealCapacity+ ...
        bipolarPlateArealCapacityBase*bipolarPlateCapacityScale;

    bipolarPlateThermalResistance = ...
        intercellBipolarPlateCount*opt.bipolarPlateThicknessM/ ...
        opt.bipolarPlateConductivity;
    intercellThermalResistance = ...
        meaThermalResistance+bipolarPlateThermalResistance;
    intercellDistance = g.Ltotal+intercellBipolarPlateCount* ...
        opt.bipolarPlateThicknessM;

    thermal.ambientTemperatureK = opt.ambientTemperatureC+273.15;
    thermal.meaArealCapacityJm2K = meaArealCapacity;
    thermal.meaThermalResistanceM2KW = meaThermalResistance;
    thermal.bipolarPlateArealCapacityBaseJm2K = ...
        bipolarPlateArealCapacityBase;
    thermal.bipolarPlateCapacityScale = bipolarPlateCapacityScale;
    thermal.cellArealCapacityJm2K = cellArealCapacity;
    thermal.cellVolumetricCapacityJm3K = cellArealCapacity/g.Ltotal;
    thermal.cellEffectiveConductivityWmK = ...
        meaEffectiveConductivity*ones(opt.nCells,1);
    thermal.intercellBipolarPlateCount = intercellBipolarPlateCount;
    thermal.intercellThermalResistanceM2KW = ...
        intercellThermalResistance;
    thermal.intercellDistanceM = intercellDistance;
    thermal.stackEffectiveConductivityWmK = ...
        intercellDistance/intercellThermalResistance;
    thermal.intercellConductancePerArea = ...
        1/intercellThermalResistance;
    thermal.endPlateArealCapacityJm2K = ...
        opt.endPlateDensity*opt.endPlateHeatCapacity* ...
        opt.endPlateThicknessM;
    thermal.cellToEndPlateThermalResistanceM2KW = ...
        0.5*meaThermalResistance+ ...
        opt.bipolarPlateThicknessM/opt.bipolarPlateConductivity+ ...
        0.5*opt.endPlateThicknessM/opt.endPlateConductivity;
    thermal.cellToEndPlateConductancePerArea = ...
        1/thermal.cellToEndPlateThermalResistanceM2KW;
    thermal.endPlateExternalResistanceM2KW = ...
        0.5*opt.endPlateThicknessM/opt.endPlateConductivity+ ...
        1/opt.endConvectionCoefficient;
    thermal.endPlateExternalConductancePerArea = ...
        1/thermal.endPlateExternalResistanceM2KW;

end


function validate_options(opt)

    positiveScalarFields = {'nCells','outputStepS','RelTol','AbsTol', ...
        'MaxStep','cellAreaM2','endConvectionCoefficient', ...
        'bipolarPlateThicknessM','bipolarPlateDensity', ...
        'bipolarPlateHeatCapacity','bipolarPlateConductivity', ...
        'endPlateThicknessM','endPlateDensity','endPlateHeatCapacity', ...
        'endPlateConductivity'};
    for k = 1:numel(positiveScalarFields)
        fieldName = positiveScalarFields{k};
        value = opt.(fieldName);
        if ~isscalar(value) || ~isfinite(value) || value <= 0
            error('q4_precooling_temperature_field:InvalidPositiveOption', ...
                '参数%s必须是正有限标量。',fieldName);
        end
    end
    temperatureFields = {'initialTemperatureC','ambientTemperatureC'};
    for k = 1:numel(temperatureFields)
        fieldName = temperatureFields{k};
        value = opt.(fieldName);
        if ~isscalar(value) || ~isfinite(value) || value <= -273.15
            error('q4_precooling_temperature_field:InvalidTemperature', ...
                '参数%s必须是高于绝对零度的有限标量。',fieldName);
        end
    end
    if ~isnumeric(opt.nCellLayer) || numel(opt.nCellLayer) ~= 5 || ...
            any(~isfinite(opt.nCellLayer)) || any(opt.nCellLayer <= 0) || ...
            any(mod(opt.nCellLayer,1) ~= 0)
        error('q4_precooling_temperature_field:InvalidGrid', ...
            'nCellLayer必须包含5个正整数。');
    end
    validBipolarCounting = {'per_cell_pair','shared_stack'};
    if (~ischar(opt.bipolarPlateCounting) && ...
            ~isstring(opt.bipolarPlateCounting)) || ...
            ~ismember(lower(char(opt.bipolarPlateCounting)), ...
            validBipolarCounting)
        error('q4_precooling_temperature_field:InvalidBipolarCounting', ...
            'bipolarPlateCounting只能为per_cell_pair或shared_stack。');
    end
    if ~isscalar(opt.plot) || ~ismember(double(opt.plot),[0 1]) || ...
            ~isscalar(opt.verbose) || ~ismember(double(opt.verbose),[0 1])
        error('q4_precooling_temperature_field:InvalidSwitch', ...
            'plot和verbose必须是逻辑标量或0/1。');
    end
    if ~(ischar(opt.saveFile) || ...
            (isstring(opt.saveFile) && isscalar(opt.saveFile)))
        error('q4_precooling_temperature_field:InvalidSaveFile', ...
            'saveFile必须是字符向量或字符串标量。');
    end

end


function output = set_default(input,fieldName,defaultValue)

    output = input;
    if ~isfield(output,fieldName) || isempty(output.(fieldName))
        output.(fieldName) = defaultValue;
    end

end


function plot_precooling_result(result,g)

    figure('Name','Question 4: stack precooling temperature field', ...
        'Color','w');
    tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

    nexttile;
    plot(result.tS/60,result.cellAverageTemperatureC,'LineWidth',1.3);
    hold on;
    plot(result.tS/60,result.endPlateTemperatureC,'--','LineWidth',1.1);
    yline(result.options.ambientTemperatureC,'k:');
    xlabel('Cooling time / min');
    ylabel('Temperature / degC');
    title('Cell-average and end-plate temperatures');
    legend('Cell 1','Cell 2','Cell 3','Cell 4','Cell 5', ...
        'Left end plate','Right end plate','Ambient', ...
        'Location','bestoutside');
    grid on;

    nexttile;
    finalField = squeeze(result.temperatureFieldC(end,:,:));
    plot(g.x_um,finalField,'LineWidth',1.2);
    xlabel('Within-cell coordinate / um');
    ylabel('Temperature / degC');
    title(sprintf('Temperature field after %.3g min', ...
        result.coolingDurationMin));
    legend('Cell 1','Cell 2','Cell 3','Cell 4','Cell 5', ...
        'Location','bestoutside');
    grid on;

end
