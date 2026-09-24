%-------------------------------------------------------------------------------
% 问题3：五片电堆外部电热丝辅助冷启动入口
%
% 本脚本只给出两组“模型连通性”候选参数，不代表优化结果。
% 纯预加热模式不加载电流；协同模式严格采用题设
% j=min(0.005*t,0.3) A/cm^2。修改每片功率和加热时间即可比较方案。
%-------------------------------------------------------------------------------

clearvars;
clc;

codeDir = fileparts(mfilename('fullpath'));
addpath(codeDir);


%% 1. 两种辅助冷启动候选方案（尚未优化）

cases(1).name = '纯预加热';
cases(1).mode = 'preheat';
cases(1).powerDensityWcm2 = [0.90 0.70 0.65 0.70 0.90];
cases(1).heatingDurationS = 120;

cases(2).name = '恒功率协同启动';
cases(2).mode = 'coheat';
cases(2).powerDensityWcm2 = [0.55 0.35 0.30 0.35 0.55];
cases(2).heatingDurationS = 80;

caseIndicesToRun = 2;


%% 2. 题设工况与求解设置

simOpt.initialTemperatureC = -30;
simOpt.ambientTemperatureC = -30;
simOpt.tMax = 180;
simOpt.outputStep = 0.5;
simOpt.MaxStep = 0.05;
simOpt.RelTol = 1e-4;
simOpt.AbsTol = 1e-8;
simOpt.nCellLayer = [2 3 4 4 2];
simOpt.plot = true;
simOpt.verbose = true;
simOpt.progressIntervalS = 5;
simOpt.useJacobianPattern = true;

% 问题2继承约束与问题3电热丝功率约束。
simOpt.qMaxCcm2 = 20;
simOpt.jMaxAcm2 = 0.5;
simOpt.minimumVoltageV = 0.30;
simOpt.maximumIceVolumeFraction = 0.99;
simOpt.maximumPoreOccupancy = 0.98;
simOpt.heaterMaxPowerDensityWcm2 = 1;

% 与Code_Q2_PhaseFix保持相同的水-冰和电压参数。
simOpt.init.lambdaCCL0 = 2.0;
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
simOpt.kFreeze = 0.05;
simOpt.kMelt = 1;
simOpt.phaseTransitionWidthK = 0.2;
simOpt.directFreezeFraction = 0.01;
simOpt.freezeNucleationLambdaFraction = 0.5;
simOpt.gammaIce = 0.1;
simOpt.endPlateModel = 'separate';
simOpt.bipolarPlateCounting = 'shared_stack';
simOpt.concentrationFactorEnd = 10;
simOpt.concentrationFactorMiddle = 1;


%% 3. 运行并生成表4所需指标

nCase = numel(caseIndicesToRun);
results = cell(nCase,1);
strategyName = strings(nCase,1);
success = false(nCase,1);
heatingDurationS = zeros(nCase,1);
totalAuxiliaryEnergyJ = zeros(nCase,1);
totalStartupTimeS = zeros(nCase,1);
minimumCellVoltageV = zeros(nCase,1);
maximumIceVolumeFraction = zeros(nCase,1);
energyAllocationJ = cell(nCase,1);
stopReason = strings(nCase,1);

for n = 1:nCase
    k = caseIndicesToRun(n);
    fprintf('\n========== %s ==========' ,cases(k).name);
    fprintf('\n');
    results{n} = pemfc_stack5_simulate_q3( ...
        cases(k).mode,cases(k).powerDensityWcm2, ...
        cases(k).heatingDurationS,simOpt);

    strategyName(n) = cases(k).name;
    success(n) = results{n}.success;
    heatingDurationS(n) = results{n}.actualHeatingTimeS;
    totalAuxiliaryEnergyJ(n) = results{n}.totalAuxiliaryEnergyJ;
    totalStartupTimeS(n) = results{n}.stopTimeS;
    minimumCellVoltageV(n) = results{n}.minimumCellVoltageV;
    maximumIceVolumeFraction(n) = ...
        results{n}.maximumIceVolumeFraction;
    energyAllocationJ{n} = results{n}.auxiliaryEnergyJCell;
    stopReason(n) = results{n}.stopReason;
end

summaryTable = table(strategyName,success,heatingDurationS, ...
    totalAuxiliaryEnergyJ,totalStartupTimeS,minimumCellVoltageV, ...
    maximumIceVolumeFraction,energyAllocationJ,stopReason);

fprintf('\n========== 问题3候选方案结果（非优化结果） ==========\n');
disp(summaryTable);

