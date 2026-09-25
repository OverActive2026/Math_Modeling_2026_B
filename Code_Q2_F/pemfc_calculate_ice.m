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
%                    .j0Ref    298.15 K参考交换电流密度，默认0.00496453591 A/m^2
%                    .Ea       活化能，默认7028.9809 J/mol
%                    .tauHyd   cCL离聚物水合时间常数，默认10 s
%                    .kFreeze  冻结系数，默认1 [1/s]（文献量级，待标定）
%                    .kMelt    融化系数，默认0 [1/s]（本组低温数据不可辨识）
%                    .gammaIce 冰覆盖修正指数，默认0（约束拟合落在下界）
%                    .init     可选初值结构体，传给pemfc_setup_ice
% @output:   result->时间、状态、电压、温度、冰量和实验对比结果
%            model->网格、状态编号、初值、输入和RHS函数句柄
%
% 使用示例：
%   result20 = pemfc_calculate_ice(-20);
%   result25 = pemfc_calculate_ice(-25,struct('plot',true));
%   opt = struct('kFreeze',1,'gammaIce',4,'tEnd',10);
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
    % [Q2F-12] 求解进度打印间隔 [s]（按**仿真物理时间**计），0 = 关闭。
    if ~isfield(simOpt,'progressIntervalS'), simOpt.progressIntervalS = 0; end
    if ~isscalar(simOpt.progressIntervalS) || ...
            ~isfinite(simOpt.progressIntervalS) || simOpt.progressIntervalS < 0
        error('pemfc_calculate_ice:InvalidProgressInterval', ...
            'progressIntervalS必须是非负有限标量（0 表示关闭）。');
    end
    if ~isfield(simOpt,'init') || isempty(simOpt.init), simOpt.init = struct(); end

    % [ACT-1] 由-20/-25 degC两个零时刻电压联合反算，避免冰和后续
    % [A-2] 默认活化参数 = 用 -20/-25 两个温度的全部实测极化点做二维拟合
    % {j0Ref, R_extra} 的结果；其中 Ea 固定为题目式(40) 给定值 67000 J/mol，
    % alpha 固定为式(39) 的 0.5，串联电阻最优解为 0（故 clResistanceScale=0）。
    % 全曲线 RMSE = 0.0611 V，两端初始点残差 −2.4 / −20.7 mV。
    if ~isfield(simOpt,'j0Ref'), simOpt.j0Ref = 0.28; end
    if ~isfield(simOpt,'Ea'), simOpt.Ea = 67000; end

    % [HYD-1] cCL离聚物由反应水逐步水合，tauHyd控制CL质子电阻
    % 从“启动初期较高”向“水合后较低”过渡的时间尺度。
    if ~isfield(simOpt,'tauHyd'), simOpt.tauHyd = 40; end
    % [DSH-1][DSH-2] 相变参数：kmCond=0 时退化为原瞬时平衡模型
    if ~isfield(simOpt,'kmCond'), simOpt.kmCond = 0; end
    if ~isfield(simOpt,'dTref'),  simOpt.dTref  = 20; end
    if ~isfield(simOpt,'nRet'),   simOpt.nRet   = 0; end
    if ~isfield(simOpt,'rRet'),   simOpt.rRet   = 0; end
    if ~isfield(simOpt,'fRet'),   simOpt.fRet   = 1; end
    % [A-2] CL 离聚物质子电阻缩放（0=关闭，1=完整FIX-5）
    % [B-2] 经扫描确定：clResistanceScale=0.60、tauHyd=26 s 可同时匹配
    % 6.4->13.2 s 的压降幅度（0.2810 vs 实验 0.2800）与温度曲线（T_RMSE 0.23 K）。
    if ~isfield(simOpt,'clResistanceScale'), simOpt.clResistanceScale = 0.66; end
    % [D-1] j0 的水合指数
    if ~isfield(simOpt,'mHyd'), simOpt.mHyd = 1.70; end
    % [E-1] fHyd 的参考 lambda 固定为 1.0 = lambdaCCL0。
    % 说明：fHyd=(lambda/lamHydRef)^mHyd 中 lamHydRef^mHyd 与 j0Ref 完全
    % 简并（fHyd = (j0Ref/lamHydRef^mHyd)*lambda^mHyd），二者不能同时辨识。
    % 因此把 lamHydRef 固定在物理参考点 lambda_CCL(t=0)=1.0，使其成为
    % 非标定量（此时 fHyd(t=0)=1），尺度全部并入 j0Ref 标定。
    if ~isfield(simOpt,'lamHydRef'), simOpt.lamHydRef = 1.0; end

    % 以下三个量是第一版冰模型唯一集中暴露的待标定量，并非已辨识结果。
    % [ICE-2] 冻结/融化相变参数
    % [F-1] kFreeze 由用户指定固定为 0.4（不再作为待标定量）。
    if ~isfield(simOpt,'kFreeze'), simOpt.kFreeze = 0.4; end
    if ~isfield(simOpt,'kMelt'), simOpt.kMelt = 0; end
    % [ICE-4] cCL冰覆盖对有效反应面积的修正指数（附件1 给定 3.5）
    if ~isfield(simOpt,'gammaIce'), simOpt.gammaIce = 3.5; end

    % [PHASE-2] 两侧气道的有限水蒸气传质系数 [m/s]；0 = 零浓度理想汇（原口径）
    if ~isfield(simOpt,'hVaporAnode'),   simOpt.hVaporAnode = 0; end
    if ~isfield(simOpt,'hVaporCathode'), simOpt.hVaporCathode = 0; end
    % [PHASE-3][PHASE-4] 低温直接冻结通道
    if ~isfield(simOpt,'directFreezeFraction'), simOpt.directFreezeFraction = 0.15; end
    if ~isfield(simOpt,'freezeNucleationLambdaFraction')
        simOpt.freezeNucleationLambdaFraction = 0.25;
    end
    % [PHASE-4b] 就地直接冻结速率系数 [1/s]
    if ~isfield(simOpt,'kFreezeDirect'), simOpt.kFreezeDirect = 1; end
    % [PHASE-5] 窄平滑相变宽度 [K]；0 = 严格分段式（题目口径）
    if ~isfield(simOpt,'phaseTransitionWidthK'), simOpt.phaseTransitionWidthK = 0; end
    % [PHASE-1] CL/PEM 界面水交换倍率（原版写死 30）
    if ~isfield(simOpt,'interfaceFactor'), simOpt.interfaceFactor = 30; end

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
    % [A-1] tauHyd<=0 表示启用平衡水合模式（lambda 直接取水活度平衡值）
    if ~isscalar(simOpt.tauHyd) || ~isfinite(simOpt.tauHyd)
        error('pemfc_calculate_ice:InvalidHydrationParameter', ...
            'tauHyd必须是有限标量（<=0 表示平衡水合模式）。');
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
    hydration.tauHyd = simOpt.tauHyd;
    % [DSH-1][DSH-2] 相变参数
    phaseOpt.kmCond = simOpt.kmCond;
    phaseOpt.dTref  = simOpt.dTref;
    phaseOpt.nRet   = simOpt.nRet;
    phaseOpt.rRet   = simOpt.rRet;
    phaseOpt.fRet   = simOpt.fRet;
    % [PHASE-2][PHASE-3][PHASE-4][PHASE-5]
    phaseOpt.hVaporAnode   = simOpt.hVaporAnode;
    phaseOpt.hVaporCathode = simOpt.hVaporCathode;
    phaseOpt.directFreezeFraction = simOpt.directFreezeFraction;
    phaseOpt.freezeNucleationLambdaFraction = simOpt.freezeNucleationLambdaFraction;
    phaseOpt.kFreezeDirect = simOpt.kFreezeDirect;
    phaseOpt.phaseTransitionWidthK = simOpt.phaseTransitionWidthK;
    phaseOpt.interfaceFactor = simOpt.interfaceFactor;
    clOpt.clResistanceScale = simOpt.clResistanceScale;
    hydOpt.mHyd = simOpt.mHyd;
    hydOpt.lamHydRef = simOpt.lamHydRef;

    tOutput = profile.time(profile.time <= simOpt.tEnd);
    if tOutput(end) < simOpt.tEnd
        tOutput(end+1,1) = simOpt.tEnd;
    end


    %% 5. 调用刚性求解器

    nonNegativeStates = [s.idx.H2(:);s.idx.O2(:);s.idx.mw(:);s.idx.mi(:); ...
        s.idx.lambdaCCL(:)];
    odeOpt = odeset( ...
        'RelTol',simOpt.RelTol, ...
        'AbsTol',simOpt.AbsTol, ...
        'MaxStep',simOpt.MaxStep, ...
        'NonNegative',nonNegativeStates);
    if simOpt.progressIntervalS > 0
        % [Q2F-12] 进度输出：每推进 progressIntervalS 秒**仿真时间**打印一行，
        % 同时给出实际等待时间。结构参考 Code_Q2_PhaseFix 的
        % ice_progress_output，便于与电堆版本对照。
        outputFunction = @(t,y,flag) ice_progress_output( ...
            t,y,flag,simOpt.progressIntervalS,simOpt.tEnd,temperatureC);
        odeOpt = odeset(odeOpt,'OutputFcn',outputFunction);
    end

    if simOpt.verbose
        fprintf('\n开始计算含冰模型%d degC工况：%.3f-%.3f s\n', ...
            temperatureC,tOutput(1),tOutput(end));
        fprintf(['冰参数（当前为待标定初值）：kFreeze=%.6g, ', ...
            'kMelt=%.6g 1/s, gammaIce=%.6g\n'], ...
            ice.kFreeze,ice.kMelt,ice.gammaIce);
        fprintf('活化参数：j0Ref=%.9g A/m^2, Ea=%.6f kJ/mol\n', ...
            activation.j0Ref,activation.Ea/1000);
        fprintf('cCL水合参数：tauHyd=%.6g s, lambdaCCL0=%.6g\n', ...
            hydration.tauHyd,init.lambdaCCL0);
    end

    rhs = @(t,x) model_rhs_ice(t,x,u,ice,hydration,activation,phaseOpt,clOpt,hydOpt,g,s);
    tic;
    % 必须在调用前清空 lastwarn，否则 lastwarn 会返回上一次运行残留的
    % 警告，把本来成功的积分误判为失败。
    lastwarn('');
    warningState = warning('query','MATLAB:ode15s:IntegrationTolNotMet');
    warning('off','MATLAB:ode15s:IntegrationTolNotMet');
    [tSol,xSol] = ode15s(rhs,tOutput,x0,odeOpt);
    [lastWarningMessage,lastWarningId] = lastwarn;
    warning(warningState.state,'MATLAB:ode15s:IntegrationTolNotMet');
    cpuTime = toc;

    % [E-4] 显式判定积分是否真正完成。
    % ode15s 在失败时有两种表现：(a) 解被截断（返回的 tSol 短于请求区间）；
    % (b) 解中出现 NaN/Inf。原代码两者都不检查，于是 NaN 会流到第 6 节的
    % 代数输出循环里，报成"水扩散系数非法"，让人误以为物性参数有问题。
    integrationTruncated = numel(tSol) == 0 || tSol(end) < tOutput(end)-1e-9;
    integrationNonFinite = any(~isfinite(tSol)) || any(~isfinite(xSol(:)));
    if integrationTruncated || integrationNonFinite || ...
            strcmp(lastWarningId,'MATLAB:ode15s:IntegrationTolNotMet')
        if any(all(isfinite(xSol),2))
            reachedTime = tSol(find(all(isfinite(xSol),2),1,'last'));
        else
            reachedTime = tOutput(1);
        end
        error('pemfc_calculate_ice:IntegrationFailed', ...
            ['ode15s 未能完成积分：只推进到 t = %.6g s（请求终点 %.6g s）。\n' ...
            '求解器信息：%s\n' ...
            '这是低温工况的数值刚性问题（冻结速率过强 / 相变 0-1 硬开关 /\n' ...
            '孔隙被冰填满造成的非光滑约束），不是物理模型发散，也不是物性\n' ...
            '参数非法。-25 degC 比 -20 degC 更容易触发，因为该温度下 cCL 开始\n' ...
            '出现液态水，题目式(9) 的液水冻结项才真正激活。可依次尝试：\n' ...
            '  1) 减小 kFreeze、directFreezeFraction 或 kFreezeDirect；\n' ...
            '  2) 放宽 simOpt.AbsTol（如 1e-6）或 RelTol；\n' ...
            '  3) 减小 simOpt.MaxStep（如 0.002）以提高刚性区分辨率；\n' ...
            '  4) 若问题出现在接近冰点的相变区，启用 phaseTransitionWidthK>0。'], ...
            reachedTime,tOutput(end),lastWarningMessage);
    end


    %% 6. 计算代数输出和冰量指标

    Nt = numel(tSol);
    Vcell = zeros(Nt,1);
    Tavg = zeros(Nt,1);
    Tmin = zeros(Nt,1);
    Tmax = zeros(Nt,1);
    lambdaMean = zeros(Nt,1);
    lambdaCCL = zeros(Nt,1);
    lambdaCCLEquilibrium = zeros(Nt,1);
    lambdaSaturationCCL = zeros(Nt,1);
    cCLIonomerWaterMassPerArea = zeros(Nt,1);
    cCLLiquidWaterMassPerArea = zeros(Nt,1);
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
    % [E-5] 适用范围诊断：记录有多少时刻的电压损失被有限化上限截断
    voltageLossClampedCount = 0;
    transportInfeasibleCount = 0;
    maxCurrentOverLimitRatio = 0;

    for k = 1:Nt
        [~,out] = model_rhs_ice( ...
            tSol(k),xSol(k,:).',u,ice,hydration,activation,phaseOpt,clOpt,hydOpt,g,s);
        Vcell(k) = out.Vcell;
        Tavg(k) = out.Tavg;
        Tmin(k) = out.Tmin;
        Tmax(k) = out.Tmax;
        lambdaMean(k) = out.lambdaMean;
        lambdaCCL(k) = out.lambdaCCL;
        lambdaCCLEquilibrium(k) = out.lambdaCCLEquilibrium;
        lambdaSaturationCCL(k) = out.water.lambdaSaturationCCL;
        cCLIonomerWaterMassPerArea(k) = sum( ...
            out.water.mIonomer(g.idx_cCL).*g.dx(g.idx_cCL));
        cCLLiquidWaterMassPerArea(k) = sum( ...
            out.water.ml(g.idx_cCL).*g.dx(g.idx_cCL));
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
        % [E-5] 适用范围诊断
        if out.voltage.voltageLossClamped
            voltageLossClampedCount = voltageLossClampedCount+1;
        end
        if out.voltage.transportInfeasible
            transportInfeasibleCount = transportInfeasibleCount+1;
        end
        if out.jLim > 0
            maxCurrentOverLimitRatio = max(maxCurrentOverLimitRatio, ...
                out.j/out.jLim);
        end
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
    result.lambdaCCL = lambdaCCL;
    result.lambdaCCLState = xSol(:,s.idx.lambdaCCL);
    result.lambdaCCLEquilibrium = lambdaCCLEquilibrium;
    result.lambdaSaturationCCL = lambdaSaturationCCL;
    result.cCLIonomerWaterMassPerArea = cCLIonomerWaterMassPerArea;
    result.cCLLiquidWaterMassPerArea = cCLLiquidWaterMassPerArea;
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
    result.hydrationParameters = hydration;

    mwPorous = xSol(:,s.idx.mw(g.idx_porous));
    iceState = xSol(:,s.idx.mi);
    poreWaterPorous = mwPorous;
    poreWaterPorous(:,s.local.mi_cCL) = ...
        poreWaterPorous(:,s.local.mi_cCL)- ...
        init.cCLIonomerWaterCoefficient.*result.lambdaCCLState;
    result.health.allFinite = all(isfinite(xSol(:))) && all(isfinite(Vcell));
    result.health.nonNegativeGasWaterIce = all(minH2 >= -1e-10) && ...
        all(minO2 >= -1e-10) && all(minMw >= -1e-10) && ...
        all(minIce >= -1e-10);
    result.health.iceNotExceedPoreWater = ...
        all(iceState(:) <= poreWaterPorous(:)+1e-8);
    % 保留旧字段名，便于已有脚本继续读取。
    result.health.iceNotExceedTotalWater = ...
        result.health.iceNotExceedPoreWater;
    result.health.iceSaturationValid = ...
        all(maxIceSaturation >= -1e-10 & maxIceSaturation <= 1+1e-8);
    result.health.currentBelowLimit = all(currentDensity < jLim);
    result.health.lambdaCCLValid = all(result.lambdaCCLState >= -1e-10 & ...
        result.lambdaCCLState <= 22+1e-8);
    % [E-5] 适用范围：false 表示结果已不可用（电流超过极限电流/离聚物干涸）
    result.health.withinModelRange = voltageLossClampedCount == 0;
    result.health.voltageLossClampedCount = voltageLossClampedCount;
    result.health.transportInfeasibleCount = transportInfeasibleCount;
    result.health.maxCurrentOverLimitRatio = maxCurrentOverLimitRatio;

    if ~result.health.withinModelRange && simOpt.verbose
        warning('pemfc_calculate_ice:OutsideModelRange', ...
            ['%d degC 工况有 %d/%d 个时刻超出模型适用范围' ...
            '（最大 j/jLim = %.4f，%d 个时刻 j>=jLim）。\n' ...
            '这些时刻的电压损失已被 [E-5] 的 %.3g V 上限截断以保证' ...
            '积分不产生 NaN，因此该段电压结果不可用于定量结论。\n' ...
            '典型原因：低温下冰堵塞孔隙使极限电流密度降到需求电流以下。'], ...
            temperatureC,voltageLossClampedCount,Nt, ...
            maxCurrentOverLimitRatio,transportInfeasibleCount,2.0);
    end

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
    model.hydrationParameters = hydration;
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
% 局部函数：求解进度输出 [Q2F-12]
% ========================================================================

function status = ice_progress_output( ...
    t,~,flag,intervalS,tEnd,temperatureC)

    % 只按**仿真物理时间**间隔打印，不在每次 Newton/RHS 调用时刷屏；
    % 同时给出实际等待时间。结构参考 Code_Q2_PhaseFix 的
    % ice_progress_output。
    persistent nextPrintTime wallClock lastSimulationTime
    status = 0;

    if strcmp(flag,'init')
        wallClock = tic;
        lastSimulationTime = t(1);
        nextPrintTime = t(1)+intervalS;
        fprintf(['[%d degC] 求解进度：%8.3f / %8.3f s，', ...
            '实际等待 %8.1f s\n'],temperatureC,t(1),tEnd,0);
    elseif isempty(flag)
        currentTime = t(end);
        lastSimulationTime = currentTime;
        if isempty(nextPrintTime)
            nextPrintTime = currentTime+intervalS;
        end
        if currentTime+10*eps(max(abs(currentTime),1)) >= nextPrintTime
            fprintf(['[%d degC] 求解进度：%8.3f / %8.3f s，', ...
                '实际等待 %8.1f s\n'],temperatureC,currentTime,tEnd, ...
                toc(wallClock));
            nextPrintTime = nextPrintTime+intervalS;
            while nextPrintTime <= ...
                    currentTime+10*eps(max(abs(currentTime),1))
                nextPrintTime = nextPrintTime+intervalS;
            end
        end
    elseif strcmp(flag,'done')
        if ~isempty(wallClock)
            fprintf(['[%d degC] 积分结束：%8.3f / %8.3f s，', ...
                '实际等待 %8.1f s\n'],temperatureC,lastSimulationTime, ...
                tEnd,toc(wallClock));
        end
        nextPrintTime = [];
        wallClock = [];
        lastSimulationTime = [];
    end

end


%% ========================================================================
% 局部函数：第一版含冰模型RHS
% RHS只调用水-冰、热、气体三个物理模块。
% ========================================================================

function [dx,out] = model_rhs_ice(t,x,u,ice,hydration,activation,phaseOpt,clOpt,hydOpt,g,s)

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
    lambdaCCLState = x(s.idx.lambdaCCL);

    cH2 = max(cH2State,0);
    cO2 = max(cO2State,0);
    mw = max(mwState,0);
    mi = max(miState,0);
    lambdaCCL = max(lambdaCCLState,0);
    trialProjectionApplied = any(cH2State < 0) || any(cO2State < 0) || ...
        any(mwState < 0) || any(miState < 0) || lambdaCCLState < 0;


    %% 3. 依次计算水-冰、热和气体状态

    [dmw,dmi,dlambdaCCL,water] = water_ice_state_ice( ...
        T,mw,mi,lambdaCCL,j,ice.kFreeze,ice.kMelt, ...
        hydration.tauHyd,phaseOpt,g,s);
    [dT,thermal] = thermal_temperature_state_ice( ...
        T,cH2,cO2,water,j,Tamb,qaux,g,s,ice.gammaIce, ...
        activation.j0Ref,activation.Ea,clOpt,hydOpt);
    [dcH2,dcO2,gas] = gas_transport_state_ice( ...
        T,cH2,cO2,mw,dmw,dT,water,thermal,j,g,s);


    %% 4. 打包状态导数

    dx = zeros(s.Nx,1);
    dx(s.idx.T) = dT;
    dx(s.idx.H2) = dcH2;
    dx(s.idx.O2) = dcO2;
    dx(s.idx.mw) = dmw;
    dx(s.idx.mi) = dmi;
    dx(s.idx.lambdaCCL) = dlambdaCCL;

    % [E-3] 状态导数必须是有限值。若某个分量算出 NaN/Inf，ode15s 会把它
    % 带进下一次试探，随后在 water_ice_state_ice 里以"水扩散系数非法"
    % 的形式报错，掩盖真正原因。这里先定位到具体分量和当时的电压分解。
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
        error('pemfc_calculate_ice:NonFiniteDerivative', ...
            ['t=%.6g s 处状态导数出现 NaN/Inf，分量：%s。\n' ...
            'j=%.6g A/m^2, Vcell=%.6g, Erev=%.6g, etaAct=%.6g, ' ...
            'etaOhm=%.6g, etaCon=%.6g, jLim=%.6g, Tmin=%.6g K。\n' ...
            '常见来源：孔隙被冰填满、etaOhm/etaCon 变为 Inf 后与 0 相乘。'], ...
            t,badComponents,j,thermal.Vcell,thermal.voltage.Erev, ...
            thermal.voltage.etaAct,thermal.voltage.etaOhm, ...
            thermal.voltage.etaCon,thermal.jLim,thermal.Tmin);
    end


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
        out.lambdaCCL = water.lambdaCCL;
        out.lambdaCCLEquilibrium = water.lambdaCCLEquilibrium;
        out.dlambdaCCL = dlambdaCCL;
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











