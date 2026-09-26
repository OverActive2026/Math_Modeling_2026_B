% Check the Q4 online rule on genuine 20/40-minute precooling states.
% Uses a tuning grid for practical runtime, NOT the official spatial grid.
clearvars;
clc;
codeDir = fileparts(mfilename('fullpath'));
addpath(codeDir);
startClock = tic;
showFigures = usejava('desktop');      % GUI shows plots; batch saves PNG

source = fullfile(codeDir,'results_control', ...
    'q4_online_rule_case1_melt.mat');
if ~isfile(source)
    error('q4_cases23:MissingController', ...
        'Run run_q4_online_rule_case1 before cases 2/3.');
end
saved = load(source,'control');
control = saved.control;
outputDir = fullfile(codeDir,'results_control');
if ~exist(outputDir,'dir'), mkdir(outputDir); end

%% Shared assumptions: only the precooling duration changes by case.
grid = [2 3 4 4 2];                  % Tuning grid, not official grid
durationsMin = [20 40];
horizonS = 80;
preOpt.nCellLayer = grid;
preOpt.initialTemperatureC = 25;
preOpt.ambientTemperatureC = -30;
preOpt.outputStepS = 5;
preOpt.MaxStep = 2;
preOpt.RelTol = 1e-7;
preOpt.AbsTol = 1e-9;
preOpt.plot = false;
preOpt.verbose = false;

simOpt.nCellLayer = grid;
simOpt.kFreeze = 0.4;
simOpt.gammaIce = 3.5;
simOpt.kMelt = 1;                      % Melt rate [1/s], adjustable
simOpt.ambientTemperatureC = -30;
simOpt.MaxStep = 0.05;
simOpt.RelTol = 1e-4;
simOpt.AbsTol = 1e-8;
simOpt.successMarginC = 0.01;
simOpt.recordStepS = 0.5;
simOpt.progressEveryS = 5;
fprintf('Physics: kFreeze=%.3g 1/s, kMelt=%.3g 1/s, gammaIce=%.3g\n', ...
    simOpt.kFreeze,simOpt.kMelt,simOpt.gammaIce);

onlineSchedule.breaksS = [0 horizonS];
onlineSchedule.powerWcm2 = ones(1,5); % Replaced by feedback each sample
onlineSchedule.sampleS = control.sampleS;
onlineSchedule.feedbackFcn = @(~,T,V,dTdt,dVdt,~) ...
    q4_online_rule_controller(T,V,dTdt,dVdt,control);
referenceSchedule.breaksS = [0 horizonS];
referenceSchedule.powerWcm2 = ones(1,5);

for k = 1:numel(durationsMin)
    caseNumber = k+1;
    durationMin = durationsMin(k);
    fprintf('\nCase %d: precooling for %g min begins.\n', ...
        caseNumber,durationMin);
    precoolingResult = q4_precooling_temperature_field( ...
        durationMin,preOpt);
    initialState = precoolingResult.initialState;
    simOpt.initialState = initialState;
    simOpt.initialStateFile = '';
    fprintf('Initial cell mean temperature / C: %s\n', ...
        mat2str(initialState.cellAverageTemperatureC(:).',4));

    fprintf('Case %d ONLINE solve begins.\n',caseNumber);
    [online,model] = q4_simulate_power_schedule_case1( ...
        onlineSchedule,simOpt);
    fprintf('Case %d full-power reference begins.\n',caseNumber);
    reference = q4_simulate_power_schedule_case1( ...
        referenceSchedule,simOpt,model);

    fprintf(['Case %d ONLINE: success=%d, t=%.3f s, E=%.3f J, ', ...
        'spread=%.3f C, Vmin=%.4f V, iceMax=%.4f\n'], ...
        caseNumber,online.success,online.stopTimeS, ...
        online.totalAuxiliaryEnergyJ, ...
        online.maximumCellTemperatureSpreadC, ...
        online.minimumCellVoltageV,online.maximumIceVolumeFraction);
    fprintf(['Case %d REFERENCE: success=%d, t=%.3f s, E=%.3f J, ', ...
        'spread=%.3f C, Vmin=%.4f V, iceMax=%.4f\n'], ...
        caseNumber,reference.success,reference.stopTimeS, ...
        reference.totalAuxiliaryEnergyJ, ...
        reference.maximumCellTemperatureSpreadC, ...
        reference.minimumCellVoltageV,reference.maximumIceVolumeFraction);
    fprintf('Online stop reason: %s\n',online.stopReason);
    fprintf('Online ice peak by cell: %s; final: %s\n', ...
        mat2str(max(online.iceVolumeFraction,[],1),4), ...
        mat2str(online.iceVolumeFraction(end,:),4));
    outputFile = fullfile(outputDir,sprintf( ...
        'q4_online_case%d_tuning_melt.mat',caseNumber));
    save(outputFile,'online','reference','initialState', ...
        'durationMin','simOpt','preOpt','control');
    fprintf('Saved case %d: %s\n',caseNumber,outputFile);
    onlinePlotFile = fullfile(outputDir,sprintf( ...
        'q4_online_case%d_tuning_melt.png',caseNumber));
    referencePlotFile = fullfile(outputDir,sprintf( ...
        'q4_full_power_case%d_tuning_melt.png',caseNumber));
    q4_plot_control_trajectory(online,onlinePlotFile, ...
        sprintf('Q4 case %d online control - tuning grid', ...
        caseNumber),showFigures);
    q4_plot_control_trajectory(reference,referencePlotFile, ...
        sprintf('Q4 case %d full-power reference - tuning grid', ...
        caseNumber),showFigures);
    fprintf('Saved trajectory plots: %s and %s\n', ...
        onlinePlotFile,referencePlotFile);
end
fprintf('Cases 2/3 total wall time: %.1f s\n',toc(startClock));
