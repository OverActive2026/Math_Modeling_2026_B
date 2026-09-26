% Q4 case 1: genuinely state-driven, sampled online heater control.
% Power is recomputed from measured T, V, dT/dt, dV/dt every sample;
% no precomputed switch times are used by this controller.
clearvars;
clc;
codeDir = fileparts(mfilename('fullpath'));
addpath(codeDir);
wallClock = tic;

%% Editable online control parameters
control.sampleS = 0.25;                % Best screened rule; 0.5 s differs by 0.02 J
control.centerEndOffC = -1.10;         % End-cell temperature for center shutoff
control.centerTaperBandC = 0.05;       % Center heater's continuous taper band
control.centerWarmReserveC = 14.0;    % Center must retain thermal reserve
control.adjacentEndOffC = -0.22;       % End-cell temperature for 2/4 shutoff
control.adjacentTaperBandC = 0.05;
control.adjacentWarmReserveC = 11.0;
control.minimumEndRiseCps = 0.35;     % Slow-end reheating trigger
control.voltageWarningV = 0.45;
control.voltageLookaheadS = 1.0;
control.iceRiskVoltageV = 0.55;
control.iceRiskDropVps = 0.02;

%% Physical and numerical settings
simOpt.nCellLayer = [2 3 4 4 2];       % Tuning grid; not final-report grid
simOpt.kFreeze = 0.4;
simOpt.gammaIce = 3.5;
simOpt.kMelt = 1;                       % Melt rate [1/s], adjustable
simOpt.MaxStep = 0.05;
simOpt.RelTol = 1e-4;
simOpt.AbsTol = 1e-8;
simOpt.successMarginC = 0.01;
simOpt.recordStepS = 0.5;
simOpt.progressEveryS = 1;             % Prints simulated, not wall-clock time
horizonS = 80;                         % Numerical search limit, not a deadline

schedule.breaksS = [0 horizonS];
schedule.powerWcm2 = ones(1,5);        % Only a placeholder; feedback replaces it
schedule.sampleS = control.sampleS;
schedule.feedbackFcn = @(~,T,V,dTdt,dVdt,~) ...
    q4_online_rule_controller(T,V,dTdt,dVdt,control);

fprintf('Physics: kFreeze=%.3g 1/s, kMelt=%.3g 1/s, gammaIce=%.3g\n', ...
    simOpt.kFreeze,simOpt.kMelt,simOpt.gammaIce);
fprintf('Starting state-driven Q4 case-1 control simulation...\n');
online = q4_simulate_power_schedule_case1(schedule,simOpt);
power = online.powerDensityWcm2;
partialMask = power > 1e-10 & power < 1-1e-10;
fprintf(['Online result: success=%d, t=%.3f s, E=%.3f J, ', ...
    'Vmin=%.4f V, iceMax=%.4f, spreadMax=%.3f C\n'], ...
    online.success,online.stopTimeS,online.totalAuxiliaryEnergyJ, ...
    online.minimumCellVoltageV,online.maximumIceVolumeFraction, ...
    online.maximumCellTemperatureSpreadC);
fprintf('Intervals with a 0<q<1 power: %d / %d\n', ...
    sum(any(partialMask,2)),size(power,1));
fprintf('Ice peak by cell: %s; final: %s\n', ...
    mat2str(max(online.iceVolumeFraction,[],1),4), ...
    mat2str(online.iceVolumeFraction(end,:),4));
fprintf('First center power reduction: %.3f s\n', ...
    first_time(online.intervalStartS,power(:,3)<1-1e-10));
fprintf('First cells-2/4 power reduction: %.3f s\n', ...
    first_time(online.intervalStartS,any(power(:,[2,4])<1-1e-10,2)));
fprintf('Actual computation wall time: %.2f s\n',toc(wallClock));

%% Save separately from the offline optimizer and old tracking controller
outputDir = fullfile(codeDir,'results_control');
if ~exist(outputDir,'dir'), mkdir(outputDir); end
save(fullfile(outputDir,'q4_online_rule_case1_melt.mat'), ...
    'online','control','simOpt','schedule');
q4_plot_control_trajectory(online,fullfile(outputDir, ...
    'q4_online_rule_case1_melt.png'), ...
    'Q4 case 1 online control - kMelt=1',usejava('desktop'));

function t = first_time(timeS,condition)
    index = find(condition,1);
    t = NaN;
    if ~isempty(index), t = timeS(index); end
end
