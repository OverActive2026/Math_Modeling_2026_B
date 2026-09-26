%-------------------------------------------------------------------------------
% 问题4：20 min和40 min预冷温度场计算入口
%-------------------------------------------------------------------------------

clearvars;
clc;

codeDir = fileparts(mfilename('fullpath'));
addpath(codeDir);
outputDir = fullfile(codeDir,'results_precooling');
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

coolingDurationMin = [20 40];

simOpt.initialTemperatureC = 25;
simOpt.ambientTemperatureC = -30;
simOpt.nCellLayer = [2 3 4 4 2]; % 提交版只使用测试网格
simOpt.outputStepS = 5;
simOpt.MaxStep = 2;
simOpt.RelTol = 1e-7;
simOpt.AbsTol = 1e-9;
simOpt.endPlateModel = 'separate';
simOpt.bipolarPlateCounting = 'shared_stack';
simOpt.plot = false;
simOpt.verbose = true;

nCase = numel(coolingDurationMin);
results = cell(nCase,1);
durationMin = zeros(nCase,1);
finalCellTemperatureC = zeros(nCase,5);
finalEndPlateTemperatureC = zeros(nCase,2);
finalCellSpreadC = zeros(nCase,1);
energyBalanceRelativeError = zeros(nCase,1);

for n = 1:nCase
    durationMin(n) = coolingDurationMin(n);
    simOpt.saveFile = fullfile(outputDir,sprintf( ...
        'q4_precooling_%gmin.mat',durationMin(n)));
    results{n} = q4_precooling_temperature_field( ...
        durationMin(n),simOpt);
    finalCellTemperatureC(n,:) = ...
        results{n}.cellAverageTemperatureC(end,:);
    finalEndPlateTemperatureC(n,:) = ...
        results{n}.endPlateTemperatureC(end,:);
    finalCellSpreadC(n) = results{n}.finalCellAverageSpreadC;
    energyBalanceRelativeError(n) = ...
        results{n}.energyBalanceRelativeError;
end

summaryTable = table(durationMin,finalCellTemperatureC, ...
    finalEndPlateTemperatureC,finalCellSpreadC, ...
    energyBalanceRelativeError);
fprintf('\n========== 问题4预冷温度场汇总 ==========\n');
disp(summaryTable);

% 后续冷启动/优化调用示例：
% data = load(fullfile(outputDir,'q4_precooling_20min.mat'),'initialState');
% startupOpt.initialState = data.initialState;
% startupOpt.ambientTemperatureC = -30;
% [startupResult,startupModel] = pemfc_stack5_simulate_q4( ...
%     'coheat',[0.55 0.35 0.30 0.35 0.55],80,startupOpt);
