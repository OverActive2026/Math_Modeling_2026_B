%-------------------------------------------------------------------------------
% @function: 读取附件2的-20/-25 degC工况并求解第一版含冰PEMFC模型
% @author:   PJ, GPT
% @date:     20260923
% @input:    temperatureC->启动温度，只允许-20或-25 [degC]
%            simOpt->可选求解与冰参数
%                    .RelTol   相对误差，默认1e-6
%                    .AbsTol   绝对误差，默认1e-8
%                    .MaxStep  最大步长，默认0.01 s
%                    .tEnd     结束时间，默认附件2全部36.6 s
%                    .plot     是否绘图，默认false
%                    .verbose  是否显示计算信息，默认true
%                    .j0Ref    298.15 K参考交换电流密度，默认0.0080757935 A/m^2
%                    .Ea       活化能，默认22204.0159 J/mol
%                    .kFreeze  冻结系数，默认0.01 [1/(K s)]（待标定）
%                    .kMelt    融化系数，默认0.01 [1/(K s)]（待标定）
%                    .gammaIce 冰覆盖修正指数，默认3.5（待标定）
%                    .init     可选初值结构体，传给pemfc_setup_ice
% @output:   result->时间、状态、电压、温度、冰量和实验对比结果
%            model->网格、状态编号、初值、输入和RHS函数句柄
%
% 使用示例：
%   result20 = pemfc_calculate_ice(-20);
%   result25 = pemfc_calculate_ice(-25,struct('plot',true));
%   opt = struct('kFreeze',0.02,'gammaIce',4,'tEnd',10);
%   result20 = pemfc_calculate_ice(-20,opt);
%-------------------------------------------------------------------------------

function [result,model] = pemfc_calculate_ice(temperatureC,simOpt)

    %% 1. 输入处理

    if nargin < 1 || isempty(temperatureC), temperatureC = -20; end
    if nargin < 2 || isempty(simOpt), simOpt = struct(); end

    if ~isscalar(temperatureC) || ~ismember(temperatureC,[-20,-25])
        error('pemfc_calculate_ice:InvalidTemperature', ...
            'temperatureC只能取-20或-25。');
    end
    if ~isstruct(simOpt) || ~isscalar(simOpt)
        error('pemfc_calculate_ice:InvalidOptions', ...
            'simOpt必须是标量结构体。');
    end

    if ~isfield(simOpt,'RelTol'), simOpt.RelTol = 1e-6; end
    if ~isfield(simOpt,'AbsTol'), simOpt.AbsTol = 1e-8; end
    if ~isfield(simOpt,'MaxStep'), simOpt.MaxStep = 0.01; end
    if ~isfield(simOpt,'plot'), simOpt.plot = false; end
    if ~isfield(simOpt,'verbose'), simOpt.verbose = true; end
    if ~isfield(simOpt,'init') || isempty(simOpt.init), simOpt.init = struct(); end

    % [ACT-1] 由-20/-25 degC两个零时刻电压联合反算，避免冰和后续
    % 传质过程干扰活化参数。二者仍保留为simOpt输入，便于继续标定。
    if ~isfield(simOpt,'j0Ref'), simOpt.j0Ref = 0.00807579349527; end
    if ~isfield(simOpt,'Ea'), simOpt.Ea = 22204.0159133; end

    % 以下三个量是第一版冰模型唯一集中暴露的待标定量，并非已辨识结果。
    % [ICE-2] 冻结/融化相变参数
    if ~isfield(simOpt,'kFreeze'), simOpt.kFreeze = 0.01; end
    if ~isfield(simOpt,'kMelt'), simOpt.kMelt = 0.01; end
    % [ICE-4] cCL冰覆盖对有效反应面积的修正指数
    if ~isfield(simOpt,'gammaIce'), simOpt.gammaIce = 3.5; end

    if ~isscalar(simOpt.RelTol) || simOpt.RelTol <= 0 || ...
            ~isscalar(simOpt.AbsTol) || simOpt.AbsTol <= 0 || ...
            ~isscalar(simOpt.MaxStep) || simOpt.MaxStep <= 0
        error('pemfc_calculate_ice:InvalidSolverOptions', ...
            'RelTol、AbsTol和MaxStep必须是正标量。');
    end
    if ~isscalar(simOpt.j0Ref) || ~isfinite(simOpt.j0Ref) || ...
            simOpt.j0Ref <= 0 || ~isscalar(simOpt.Ea) || ...
            ~isfinite(simOpt.Ea) || simOpt.Ea < 0
        error('pemfc_calculate_ice:InvalidActivationParameters', ...
            'j0Ref必须为正有限标量，Ea必须为非负有限标量。');
    end
    if ~isscalar(simOpt.kFreeze) || ~isfinite(simOpt.kFreeze) || ...
            simOpt.kFreeze < 0 || ~isscalar(simOpt.kMelt) || ...
            ~isfinite(simOpt.kMelt) || simOpt.kMelt < 0 || ...
            ~isscalar(simOpt.gammaIce) || ~isfinite(simOpt.gammaIce) || ...
            simOpt.gammaIce < 0
        error('pemfc_calculate_ice:InvalidIceParameters', ...
            'kFreeze、kMelt和gammaIce必须是非负有限标量。');
    end
    if ~isstruct(simOpt.init) || ~isscalar(simOpt.init)
        error('pemfc_calculate_ice:InvalidInitialOptions', ...
            'simOpt.init必须是标量结构体。');
    end


    %% 2. 读取附件2实验数据

    codeDir = fileparts(mfilename('fullpath'));
    rootDir = fileparts(codeDir);
    attachmentDir = fullfile(rootDir, ...
        '氢燃料电池低温冷启动建模与控制策略研究  附件');
    filePath = fullfile(attachmentDir,'附件2.xlsx');
    sheetName = sprintf('%d℃',temperatureC);

    if ~isfile(filePath)
        error('pemfc_calculate_ice:MissingExperimentFile', ...
            '未找到附件2：%s',filePath);
    end

    raw = readmatrix(filePath,'Sheet',sheetName,'Range','A3:E10000');
    raw = raw(~all(isnan(raw),2),:);
    if size(raw,2) < 5 || isempty(raw)
        error('pemfc_calculate_ice:InvalidExperimentData', ...
            '工作表%s没有读取到5列有效数据。',sheetName);
    end

    valid = all(isfinite(raw(:,1:5)),2);
    raw = raw(valid,1:5);
    profile.time = raw(:,1);
    profile.current = raw(:,2);
    profile.voltage = raw(:,3);
    profile.temperatureC = raw(:,4);
    profile.currentDensityAcm2 = raw(:,5);
    profile.currentDensity = raw(:,5)*1e4;
    profile.sheetName = sheetName;
    profile.filePath = filePath;

    if numel(profile.time) < 2 || any(diff(profile.time) <= 0)
        error('pemfc_calculate_ice:InvalidExperimentTime', ...
            '实验时间必须严格递增且至少包含两个点。');
    end
    if ~isfield(simOpt,'tEnd') || isempty(simOpt.tEnd)
        simOpt.tEnd = profile.time(end);
    end
    if ~isscalar(simOpt.tEnd) || ~isfinite(simOpt.tEnd) || ...
            simOpt.tEnd <= profile.time(1) || simOpt.tEnd > profile.time(end)
        error('pemfc_calculate_ice:InvalidEndTime', ...
            'tEnd必须位于(%.6g,%.6g] s。', ...
            profile.time(1),profile.time(end));
    end


    %% 3. 建立对应温度下的模型初值

    startTemperature = temperatureC+273.15;
    initOpt = simOpt.init;
    initOpt.T0 = startTemperature;
    [g,s,x0,init] = pemfc_setup_ice([],initOpt);


    %% 4. 构造电流输入和冰参数

    u.j = @(timeQuery) interp1(profile.time,profile.currentDensity, ...
        min(max(timeQuery,profile.time(1)),profile.time(end)),'linear');
    u.Tamb = startTemperature;
    u.qaux = 0;

    ice.kFreeze = simOpt.kFreeze;
    ice.kMelt = simOpt.kMelt;
    ice.gammaIce = simOpt.gammaIce;
    activation.j0Ref = simOpt.j0Ref;
    activation.Ea = simOpt.Ea;

    tOutput = profile.time(profile.time <= simOpt.tEnd);
    if tOutput(end) < simOpt.tEnd
        tOutput(end+1,1) = simOpt.tEnd;
    end


    %% 5. 调用刚性求解器

    nonNegativeStates = [s.idx.H2(:);s.idx.O2(:);s.idx.mw(:);s.idx.mi(:)];
    odeOpt = odeset( ...
        'RelTol',simOpt.RelTol, ...
        'AbsTol',simOpt.AbsTol, ...
        'MaxStep',simOpt.MaxStep, ...
        'NonNegative',nonNegativeStates);

    if simOpt.verbose
        fprintf('\n开始计算含冰模型%d degC工况：%.3f-%.3f s\n', ...
            temperatureC,tOutput(1),tOutput(end));
        fprintf(['冰参数（当前为待标定初值）：kFreeze=%.6g, ', ...
            'kMelt=%.6g 1/(K s), gammaIce=%.6g\n'], ...
            ice.kFreeze,ice.kMelt,ice.gammaIce);
        fprintf('活化参数：j0Ref=%.9g A/m^2, Ea=%.6f kJ/mol\n', ...
            activation.j0Ref,activation.Ea/1000);
    end

    rhs = @(t,x) model_rhs_ice(t,x,u,ice,activation,g,s);
    tic;
    [tSol,xSol] = ode15s(rhs,tOutput,x0,odeOpt);
    cpuTime = toc;


    %% 6. 计算代数输出和冰量指标

    Nt = numel(tSol);
    Vcell = zeros(Nt,1);
    Tavg = zeros(Nt,1);
    Tmin = zeros(Nt,1);
    Tmax = zeros(Nt,1);
    lambdaMean = zeros(Nt,1);
    currentDensity = zeros(Nt,1);
    jLim = zeros(Nt,1);
    minH2 = zeros(Nt,1);
    minO2 = zeros(Nt,1);
    minMw = zeros(Nt,1);
    minIce = zeros(Nt,1);
    iceMassPerArea = zeros(Nt,1);
    maxIceSaturation = zeros(Nt,1);
    iceSaturationCCL = zeros(Nt,1);
    iceAreaFactor = zeros(Nt,1);
    latentHeatPerArea = zeros(Nt,1);

    for k = 1:Nt
        [~,out] = model_rhs_ice( ...
            tSol(k),xSol(k,:).',u,ice,activation,g,s);
        Vcell(k) = out.Vcell;
        Tavg(k) = out.Tavg;
        Tmin(k) = out.Tmin;
        Tmax(k) = out.Tmax;
        lambdaMean(k) = out.lambdaMean;
        currentDensity(k) = out.j;
        jLim(k) = out.jLim;
        minH2(k) = out.minH2;
        minO2(k) = out.minO2;
        minMw(k) = out.minMw;
        minIce(k) = out.minIce;
        iceMassPerArea(k) = sum(out.water.miFull.*g.dx);
        maxIceSaturation(k) = max(out.water.sIce);
        iceSaturationCCL(k) = out.iceSaturationCCL;
        iceAreaFactor(k) = out.iceAreaFactor;
        latentHeatPerArea(k) = sum(out.qPhase.*g.dx);
    end

    experimentVoltage = interp1(profile.time,profile.voltage,tSol,'linear');
    experimentTemperatureC = interp1( ...
        profile.time,profile.temperatureC,tSol,'linear');


    %% 7. 汇总计算结果

    result.caseTemperatureC = temperatureC;
    result.t = tSol;
    result.x = xSol;
    result.Vcell = Vcell;
    result.Tavg = Tavg;
    result.TavgC = Tavg-273.15;
    result.Tmin = Tmin;
    result.Tmax = Tmax;
    result.lambdaMean = lambdaMean;
    result.currentDensity = currentDensity;
    result.currentDensityAcm2 = currentDensity/1e4;
    result.jLim = jLim;
    result.minH2 = minH2;
    result.minO2 = minO2;
    result.minMw = minMw;
    result.minIce = minIce;
    result.iceMass = xSol(:,s.idx.mi);
    result.iceMassPerArea = iceMassPerArea;
    result.maxIceSaturation = maxIceSaturation;
    result.iceSaturationCCL = iceSaturationCCL;
    result.iceAreaFactor = iceAreaFactor;
    result.latentHeatPerArea = latentHeatPerArea;
    result.cpuTime = cpuTime;
    result.experiment = profile;
    result.experimentVoltageAtSolutionTime = experimentVoltage;
    result.experimentTemperatureCAtSolutionTime = experimentTemperatureC;
    result.voltageRMSE = sqrt(mean((Vcell-experimentVoltage).^2));
    result.temperatureRMSE = sqrt(mean((result.TavgC- ...
        experimentTemperatureC).^2));
    result.solver = simOpt;
    result.iceParameters = ice;
    result.activationParameters = activation;

    mwPorous = xSol(:,s.idx.mw(g.idx_porous));
    iceState = xSol(:,s.idx.mi);
    result.health.allFinite = all(isfinite(xSol(:))) && all(isfinite(Vcell));
    result.health.nonNegativeGasWaterIce = all(minH2 >= -1e-10) && ...
        all(minO2 >= -1e-10) && all(minMw >= -1e-10) && ...
        all(minIce >= -1e-10);
    result.health.iceNotExceedTotalWater = ...
        all(iceState(:) <= mwPorous(:)+1e-8);
    result.health.iceSaturationValid = ...
        all(maxIceSaturation >= -1e-10 & maxIceSaturation <= 1+1e-8);
    result.health.currentBelowLimit = all(currentDensity < jLim);

    if ~result.health.allFinite
        error('pemfc_calculate_ice:NonFiniteSolution', ...
            '求解结果包含NaN或Inf。');
    end

    model.g = g;
    model.s = s;
    model.x0 = x0;
    model.init = init;
    model.u = u;
    model.iceParameters = ice;
    model.activationParameters = activation;
    model.rhs = rhs;

    if simOpt.verbose
        fprintf('计算完成：耗时%.3f s，最终电压%.6f V，平均温度%.6f degC\n', ...
            cpuTime,Vcell(end),result.TavgC(end));
        fprintf('最终冰面密度%.6e kg/m^2，cCL平均冰饱和度%.6f\n', ...
            iceMassPerArea(end),iceSaturationCCL(end));
        fprintf('实验对比：电压RMSE=%.6e V，温度RMSE=%.6e degC\n', ...
            result.voltageRMSE,result.temperatureRMSE);
    end


    %% 8. 可选绘图

    if simOpt.plot
        figure('Name',sprintf('PEMFC ICE %d degC cold start',temperatureC));
        tiledlayout(2,2);

        nexttile;
        plot(tSol,currentDensity/1e4,'LineWidth',1.5);
        xlabel('Time / s'); ylabel('Current density / A cm^{-2}');
        title('Current input'); grid on;

        nexttile;
        plot(tSol,Vcell,'LineWidth',1.5); hold on;
        plot(tSol,experimentVoltage,'--','LineWidth',1.2);
        xlabel('Time / s'); ylabel('Cell voltage / V');
        title('Cell voltage');
        legend('Model','Experiment','Location','best'); grid on;

        nexttile;
        plot(tSol,result.TavgC,'LineWidth',1.5); hold on;
        plot(tSol,experimentTemperatureC,'--','LineWidth',1.2);
        xlabel('Time / s'); ylabel('Temperature / ^\circC');
        title('Average temperature');
        legend('Model','Experiment','Location','best'); grid on;

        nexttile;
        plot(tSol,iceSaturationCCL,'LineWidth',1.5);
        xlabel('Time / s'); ylabel('cCL ice saturation');
        title('Ice accumulation in cCL'); grid on;
    end

end


%% ========================================================================
% 局部函数：第一版含冰模型RHS
% RHS只调用水-冰、热、气体三个物理模块。
% ========================================================================

function [dx,out] = model_rhs_ice(t,x,u,ice,activation,g,s)

    %% 1. 读取当前输入

    if isa(u.j,'function_handle'), j = u.j(t); else, j = u.j; end
    if isa(u.Tamb,'function_handle'), Tamb = u.Tamb(t); else, Tamb = u.Tamb; end
    if isa(u.qaux,'function_handle'), qaux = u.qaux(t); else, qaux = u.qaux; end

    if ~isscalar(j) || ~isfinite(j) || j < 0
        error('pemfc_calculate_ice:InvalidCurrent', ...
            '电流密度必须是非负有限标量。');
    end
    if ~isscalar(Tamb) || ~isfinite(Tamb) || Tamb <= 0
        error('pemfc_calculate_ice:InvalidAmbientTemperature', ...
            '环境温度必须是正有限标量。');
    end
    if isscalar(qaux), qaux = qaux*ones(g.N,1); else, qaux = qaux(:); end


    %% 2. 拆分状态并处理求解器试探值

    T = x(s.idx.T);
    cH2State = x(s.idx.H2);
    cO2State = x(s.idx.O2);
    mwState = x(s.idx.mw);
    miState = x(s.idx.mi);

    cH2 = max(cH2State,0);
    cO2 = max(cO2State,0);
    mw = max(mwState,0);
    mi = max(miState,0);
    trialProjectionApplied = any(cH2State < 0) || any(cO2State < 0) || ...
        any(mwState < 0) || any(miState < 0);


    %% 3. 依次计算水-冰、热和气体状态

    [dmw,dmi,water] = water_ice_state_ice( ...
        T,mw,mi,j,ice.kFreeze,ice.kMelt,g,s);
    [dT,thermal] = thermal_temperature_state_ice( ...
        T,cH2,cO2,water,j,Tamb,qaux,g,s,ice.gammaIce, ...
        activation.j0Ref,activation.Ea);
    [dcH2,dcO2,gas] = gas_transport_state_ice( ...
        T,cH2,cO2,mw,dmw,dT,water,thermal,j,g,s);


    %% 4. 打包状态导数

    dx = zeros(s.Nx,1);
    dx(s.idx.T) = dT;
    dx(s.idx.H2) = dcH2;
    dx(s.idx.O2) = dcO2;
    dx(s.idx.mw) = dmw;
    dx(s.idx.mi) = dmi;


    %% 5. 汇总求解、标定和控制所需输出

    if nargout > 1
        out.t = t;
        out.j = j;
        out.Tamb = Tamb;
        out.Vcell = thermal.Vcell;
        out.voltage = thermal.voltage;
        out.Tavg = thermal.Tavg;
        out.Tmin = thermal.Tmin;
        out.Tmax = thermal.Tmax;
        out.dT = dT;
        out.cH2_aCL = thermal.cH2_aCL;
        out.cO2_cCL = thermal.cO2_cCL;
        out.pH2 = thermal.pH2;
        out.pO2 = thermal.pO2;
        out.DH2 = gas.DH2;
        out.DO2 = gas.DO2;
        out.DO2Equivalent = thermal.DO2Equivalent;
        out.DO2ArithmeticCCL = thermal.DO2ArithmeticCCL;
        out.oxygenDiffusionResistance = thermal.oxygenDiffusionResistance;
        out.eps_g = gas.eps_g;
        out.depsGdt = gas.depsGdt;
        out.JH2 = gas.JH2;
        out.JO2 = gas.JO2;
        out.water = water;
        out.Nw = water.Nw;
        out.dmw = dmw;
        out.dmi = dmi;
        out.lambdaMean = water.lambdaMean;
        out.reaction = gas.reaction;
        out.sources = gas.source;
        out.qgen = thermal.qgen;
        out.qPhase = thermal.qPhase;
        out.qCond = thermal.qCond;
        out.qHeatFace = thermal.qHeatFace;
        out.thermalProperties = thermal.thermalProperties;
        out.thermalFV = thermal.thermalFV;
        out.jLim = thermal.jLim;
        out.transportValid = thermal.transportValid;
        out.iceSaturationCCL = thermal.iceSaturationCCL;
        out.iceAreaFactor = thermal.iceAreaFactor;
        out.minH2 = min(cH2State);
        out.minO2 = min(cO2State);
        out.minMw = min(mwState);
        out.minIce = min(miState);
        out.trialProjectionApplied = trialProjectionApplied || ...
            water.iceProjectionApplied;
    end

end
