% Recheck the time boundary and one failed perturbation on the same grid.
codeDir = fileparts(mfilename('fullpath'));
addpath(codeDir);
outDir = fullfile(codeDir,'results');
a = load(fullfile(outDir,'time_priority_candidates.mat'), ...
    'pre','co','opt','preHeatingDurationS','coHeatingDurationS');
opt = a.opt;
opt.RelTol = 1e-5;
opt.AbsTol = 1e-9;
opt.MaxStep = 0.025;
opt.plot = false;
opt.verbose = false;
opt.enforceChargeCap = true;

preTight = pemfc_stack5_simulate_q3( ...
    'preheat',ones(1,5),a.preHeatingDurationS,opt);
fprintf('preheat tight: success=%d, minT=%.7f C, E=%.6f J\n', ...
    preTight.success,min(preTight.TavgC(end,:)), ...
    preTight.totalAuxiliaryEnergyJ);

coTight = pemfc_stack5_simulate_q3( ...
    'coheat',ones(1,5),a.coHeatingDurationS,opt);
fprintf('coheat tight: success=%d, ts=%.7f s, minV=%.7f V\n', ...
    coTight.success,coTight.successTimeS, ...
    coTight.minimumCellVoltageV);

q4 = [1 1 1 0.95 1];
cutoffS = a.co.successTimeS-0.001;
opt.tMax = cutoffS;
q4Tight = [];
q4Status = "numerical_failure";
q4Message = "";
try
    q4Tight = pemfc_stack5_simulate_q3( ...
        'coheat',q4,cutoffS,opt);
    if q4Tight.success
        q4Status = "faster_success";
    elseif q4Tight.solverTerminatedUnexpectedly
        q4Status = "numerical_failure";
    else
        q4Status = "not_faster";
    end
    fprintf('q4 tight: %s, minT=%.7f C\n', ...
        q4Status,min(q4Tight.TavgC(end,:)));
catch ME
    q4Message = string(ME.identifier)+": "+string(ME.message);
    fprintf('q4 tight: %s\n',q4Message);
end
save(fullfile(outDir,'time_numerical_check.mat'), ...
    'preTight','coTight','q4Tight','q4Status','q4Message', ...
    'opt','cutoffS','-v7.3');
