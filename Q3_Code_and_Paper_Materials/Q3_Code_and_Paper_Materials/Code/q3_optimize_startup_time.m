% Time-priority Q3 search using the existing physical model unchanged.
% With independent upper bounds q_k<=1, try the maximum heating vector.
% The output is a verified candidate, not a global monotonicity proof.
codeDir = fileparts(mfilename('fullpath'));
addpath(codeDir);
cfg = q3_optimization_config('smoke');
opt = cfg.simOpt;
opt.plot = false;
opt.verbose = false;
opt.enforceChargeCap = true;
q = ones(1,5);

% Locate the first five-cell temperature crossing while preheating.
preProbe = pemfc_stack5_simulate_q3('preheat',q,35,opt);
preEvents = preProbe.event.time(preProbe.event.index == 1);
assert(~isempty(preEvents),'No preheat temperature crossing by 35 s.');
preThresholdTimeS = preEvents(1);
preHeatingDurationS = preThresholdTimeS+1e-5;
pre = pemfc_stack5_simulate_q3( ...
    'preheat',q,preHeatingDurationS,opt);
if ~pre.success && ~pre.solverTerminatedUnexpectedly && ...
        min(pre.TavgC(end,:))>0 && ...
        min(pre.TavgC(end,:))<=opt.temperatureEventMarginK
    preHeatingDurationS=preHeatingDurationS+1e-3;
    pre=pemfc_stack5_simulate_q3('preheat',q,preHeatingDurationS,opt);
end
assert(pre.success && all(pre.currentDensityAcm2 == 0));
assert(all(pre.TavgC(end,:) > 0));

% In coheating, locate success and then recheck a nearly equal planned
% heating duration. The prescribed current starts at t=0 in both runs.
coProbe = pemfc_stack5_simulate_q3('coheat',q,60,opt);
assert(coProbe.success);
coHeatingDurationS = coProbe.successTimeS+0.01;
co = pemfc_stack5_simulate_q3( ...
    'coheat',q,coHeatingDurationS,opt);
assert(co.success && co.chargeUsedCcm2 <= 20);
assert(abs(co.currentDensityAcm2(end)- ...
    min(0.005*co.stopTimeS,0.3)) < 1e-8);

mode = ["preheat";"coheat"];
q1 = ones(2,1); q2 = q1; q3 = q1; q4 = q1; q5 = q1;
requestedHeatingDurationS = ...
    [preHeatingDurationS;coHeatingDurationS];
heatingDurationS = [pre.actualHeatingTimeS;co.actualHeatingTimeS];
startupTimeS = [pre.successTimeS;co.successTimeS];
E1J = [pre.auxiliaryEnergyJCell(1);co.auxiliaryEnergyJCell(1)];
E2J = [pre.auxiliaryEnergyJCell(2);co.auxiliaryEnergyJCell(2)];
E3J = [pre.auxiliaryEnergyJCell(3);co.auxiliaryEnergyJCell(3)];
E4J = [pre.auxiliaryEnergyJCell(4);co.auxiliaryEnergyJCell(4)];
E5J = [pre.auxiliaryEnergyJCell(5);co.auxiliaryEnergyJCell(5)];
ETotalJ = [pre.totalAuxiliaryEnergyJ;co.totalAuxiliaryEnergyJ];
maxIce = [pre.maximumIceVolumeFraction;co.maximumIceVolumeFraction];
minimumVoltageV = [pre.minimumCellVoltageV;co.minimumCellVoltageV];
chargeCcm2 = [pre.chargeUsedCcm2;co.chargeUsedCcm2];
success = [pre.success;co.success];
assessmentScope = ["preheat_threshold_only";"loaded_startup"];
summary = table(mode,q1,q2,q3,q4,q5,requestedHeatingDurationS, ...
    heatingDurationS, ...
    startupTimeS,E1J,E2J,E3J,E4J,E5J,ETotalJ,maxIce, ...
    minimumVoltageV,chargeCcm2,success,assessmentScope);
outDir = fullfile(codeDir,'results');
if ~exist(outDir,'dir'), mkdir(outDir); end
writetable(summary,fullfile(outDir,'time_priority_candidates.csv'), ...
    'Encoding','UTF-8');
save(fullfile(outDir,'time_priority_candidates.mat'), ...
    'summary','preProbe','pre','coProbe','co','opt', ...
    'preThresholdTimeS','preHeatingDurationS', ...
    'coHeatingDurationS','-v7.3');
disp(summary);
fprintf('Preheat crossing %.6f s, verified heater stop %.6f s.\n', ...
    preThresholdTimeS,preHeatingDurationS);
