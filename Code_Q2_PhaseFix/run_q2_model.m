%-------------------------------------------------------------------------------
% 问题2：5片电堆自冷启动调用脚本
%
% 本脚本先用一组“待优化初值”检查三种加载策略。
% 它不把这三组参数声称为最优结果；等热耦合假设和
% 首尾浓差修正经灵敏度分析确认后，再接优化器。
%-------------------------------------------------------------------------------

clearvars;
clc;

codeDir = fileparts(mfilename('fullpath'));
addpath(codeDir);


%% 1. 待优化的三类电流加载策略

strategies(1).name = '恒流';
strategies(1).parameter = struct('type','constant', ...
    'currentAcm2',0.25);

strategies(2).name = '线性升载';
strategies(2).parameter = struct('type','linear', ...
    'rampRateAcm2s',0.005,'plateauAcm2',0.30);

strategies(3).name = '三段阶梯';
strategies(3).parameter = struct('type','step', ...
    'levelsAcm2',[0.10 0.25 0.40], ...
    'switchTimesS',[8 22]);

% 默认只跑线性升载，用于快速检查。改为1:3即比较全部策略。
strategyCasesToRun = 2;


%% 2. 题目工况和求解设置

initialTemperatureC = -10;

simOpt.tMax = 180;
simOpt.outputStep = 0.2;
simOpt.MaxStep = 0.02;
simOpt.RelTol = 1e-6;
simOpt.AbsTol = 1e-8;
simOpt.plot = true;
simOpt.verbose = true;

% 问题2约束
simOpt.qMaxCcm2 = 20;
simOpt.jMaxAcm2 = 0.5;
simOpt.minimumVoltageV = 0.30;
simOpt.maximumIceVolumeFraction = 0.99;


%% 3. 当前Code_ICE单电池参数

simOpt.j0Ref = 0.527582014454;
simOpt.Ea = 8000;
simOpt.fHydDry = 0.00705033754574;
simOpt.lambdaHydOn = 3.5;
simOpt.lambdaHydWet = 8.53122413758;
simOpt.nHyd = 3.0;
simOpt.ionomerBruggemanExponent = 1.77972166611;
simOpt.tauHyd = 10;
simOpt.hVaporCathode = 0.01;
simOpt.hVaporAnode = 0.01;
simOpt.interfaceRateFactor = 0.01;
simOpt.init.lambdaCCL0 = 2.0;
simOpt.endPlateModel = 'separate';
simOpt.bipolarPlateCounting = 'shared_stack';
simOpt.kFreeze = 0.05;
simOpt.kMelt = 1;
% [NUM-PHASE] 相变开关的窄平滑层，设0可恢复题目严格分段式。
simOpt.phaseTransitionWidthK = 0.2;
simOpt.useJacobianPattern = true;
simOpt.directFreezeFraction = 0.01;
simOpt.freezeNucleationLambdaFraction = 0.5;
simOpt.gammaIce = 0.1;

% [Q2-4] 附件1给定的位置浓差极化系数。
% 建议后续将首尾系数在1~10内做灵敏度分析。
simOpt.concentrationFactorEnd = 10;
simOpt.concentrationFactorMiddle = 1;


%% 4. 计算所选策略

results = cell(numel(strategyCasesToRun),1);
models = cell(numel(strategyCasesToRun),1);

for n = 1:numel(strategyCasesToRun)
    caseIndex = strategyCasesToRun(n);
    fprintf('\n========== %s策略 ==========\n',strategies(caseIndex).name);
    [results{n},models{n}] = pemfc_stack5_simulate( ...
        strategies(caseIndex).parameter,initialTemperatureC,simOpt);
end


%% 5. 生成题目表3需要的核心指标

nCase = numel(results);
strategyName = strings(nCase,1);
success = false(nCase,1);
startTimeS = nan(nCase,1);
chargeUsedCcm2 = zeros(nCase,1);
maximumCurrentDensityAcm2 = zeros(nCase,1);
minimumCellVoltageV = zeros(nCase,1);
maximumIceVolumeFraction = zeros(nCase,1);
stopReason = strings(nCase,1);

for n = 1:nCase
    caseIndex = strategyCasesToRun(n);
    strategyName(n) = strategies(caseIndex).name;
    success(n) = results{n}.success;
    startTimeS(n) = results{n}.successTimeS;
    chargeUsedCcm2(n) = results{n}.chargeUsedCcm2;
    maximumCurrentDensityAcm2(n) = ...
        results{n}.maximumCurrentDensityAcm2;
    minimumCellVoltageV(n) = results{n}.minimumCellVoltageV;
    maximumIceVolumeFraction(n) = ...
        results{n}.maximumIceVolumeFraction;
    stopReason(n) = results{n}.stopReason;
end


summaryTable = table(strategyName,success,startTimeS,chargeUsedCcm2, ...
    maximumCurrentDensityAcm2,minimumCellVoltageV, ...
    maximumIceVolumeFraction,stopReason);

fprintf('\n========== 问题2初始策略结果 ==========\n');
disp(summaryTable);
