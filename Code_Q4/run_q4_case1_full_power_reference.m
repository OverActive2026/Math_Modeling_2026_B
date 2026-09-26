% Q4(1) case 1: test-grid full-power feasibility reference.
% This is NOT the optimized constant-power cooperative strategy of Q3.
clearvars;
clc;
codeDir = fileparts(mfilename('fullpath'));
addpath(codeDir);
sourceFile = fullfile(codeDir,'results_control', ...
    'q4_online_rule_case1_melt.mat');
if ~isfile(sourceFile)
    error('q4_case1_ref:MissingOnlineResult', ...
        'Run run_q4_online_rule_case1 first.');
end
saved = load(sourceFile,'simOpt');
simOpt = saved.simOpt;
if ~isequal(simOpt.nCellLayer,[2 3 4 4 2])
    error('q4_case1_ref:WrongGrid', ...
        'This submission uses only the [2 3 4 4 2] test grid.');
end
simOpt.progressEveryS = 1;
schedule.breaksS = [0 80];
schedule.powerWcm2 = ones(1,5);
reference = q4_simulate_power_schedule_case1(schedule,simOpt);
fprintf(['Case 1 full-power reference: success=%d, t=%.3f s, ', ...
    'E=%.3f J, Vmin=%.4f V, iceMax=%.4f\n'], ...
    reference.success,reference.stopTimeS, ...
    reference.totalAuxiliaryEnergyJ,reference.minimumCellVoltageV, ...
    reference.maximumIceVolumeFraction);
outputDir = fullfile(codeDir,'results_control');
if ~exist(outputDir,'dir'), mkdir(outputDir); end
save(fullfile(outputDir,'q4_full_power_case1_test_grid_melt.mat'), ...
    'reference','simOpt');
q4_plot_control_trajectory(reference,fullfile(outputDir, ...
    'q4_full_power_case1_test_grid_melt.png'), ...
    'Q4 case 1 full-power reference - test grid',usejava('desktop'));
