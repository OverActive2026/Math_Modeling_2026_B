function summary = q3_probe_25s()
% Fixed-time feasibility probes for the 25 s target; physical model unchanged.
codeDir = fileparts(mfilename('fullpath'));
addpath(codeDir);
outDir = fullfile(codeDir,'results');
saved = load(fullfile(outDir,'time_priority_candidates.mat'),'coProbe');
cfg = q3_optimization_config('smoke');
opt = cfg.simOpt;
opt.tMax = 25;
opt.plot = false;
opt.verbose = false;
opt.enforceChargeCap = true;
q = [1 1 1 1 1; ...
     1 1 0.8 1 1; ...
     1 0.8 1 0.8 1; ...
     1 0.8 0.8 0.8 1];
caseName = ["all_max";"center_reduced"; ...
    "adjacent_reduced";"end_focused"];
minCellTemperatureC = nan(4,1);
endCellTemperatureC = nan(4,2);
minimumVoltageV = nan(4,1);
maximumIceVolumeFraction = nan(4,1);
chargeCcm2 = nan(4,1);
observedStopTimeS = nan(4,1);
successBy25s = false(4,1);
status = strings(4,1);
for i = 1:4
    if i == 1
        r = saved.coProbe;
        idx = find(abs(r.t-25)<1e-9,1);
        assert(~isempty(idx));
    else
        try
            r = pemfc_stack5_simulate_q3('coheat',q(i,:),25,opt);
            idx = numel(r.t);
        catch ME
            status(i) = string(ME.identifier);
            fprintf('%s: %s\n',caseName(i),status(i));
            write_checkpoint();
            continue;
        end
    end
    observedStopTimeS(i) = r.stopTimeS;
    if r.solverTerminatedUnexpectedly
        status(i) = "numerical_failure_before_25s";
        fprintf('%s: numerical solver stopped at %.5f s\n', ...
            caseName(i),r.stopTimeS);
        write_checkpoint();
        continue;
    end
    if r.stopTimeS < 25-1e-7 && ~r.success
        status(i) = "physical_failure_before_25s";
        fprintf('%s: physical failure at %.5f s\n', ...
            caseName(i),r.stopTimeS);
        write_checkpoint();
        continue;
    end
    minCellTemperatureC(i) = min(r.TavgC(idx,:));
    endCellTemperatureC(i,:) = r.TavgC(idx,[1 5]);
    minimumVoltageV(i) = min(r.Vcell(1:idx,:),[],'all');
    maximumIceVolumeFraction(i) = max( ...
        r.maxIceVolumeFractionCell(1:idx,:),[],'all');
    chargeCcm2(i) = r.chargeCcm2(idx);
    successBy25s(i) = r.success && r.successTimeS <= 25;
    if i == 1
        successBy25s(i) = any(r.event.index == 1 & r.event.time <= 25);
    end
    status(i) = "evaluated";
    fprintf('%s: minT=%.5f C, ends=%s, success=%d\n', ...
        caseName(i),minCellTemperatureC(i), ...
        mat2str(endCellTemperatureC(i,:),5),successBy25s(i));
    write_checkpoint();
end
summary = make_table();
save(fullfile(outDir,'q3_25s_probes.mat'), ...
    'summary','q','opt','-v7.3');

function write_checkpoint()
    partial = make_table();
    writetable(partial,fullfile(outDir,'q3_25s_probes.csv'), ...
        'Encoding','UTF-8');
end

function tableOut = make_table()
    tableOut = table(caseName,q(:,1),q(:,2),q(:,3),q(:,4),q(:,5), ...
        minCellTemperatureC,endCellTemperatureC(:,1), ...
        endCellTemperatureC(:,2),minimumVoltageV, ...
        maximumIceVolumeFraction,chargeCcm2,observedStopTimeS, ...
        successBy25s,status, ...
        'VariableNames',{'caseName','q1','q2','q3','q4','q5', ...
        'minCellTemperatureC','T1C','T5C','minimumVoltageV', ...
        'maximumIceVolumeFraction','chargeCcm2','observedStopTimeS', ...
        'successBy25s','status'});
end
end
