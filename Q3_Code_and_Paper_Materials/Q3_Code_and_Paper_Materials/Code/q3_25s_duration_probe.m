function summary = q3_25s_duration_probe()
% Check the sixth decision variable: heater shutoff before 25 s.
    codeDir = fileparts(mfilename('fullpath'));
    addpath(codeDir);
    outDir = fullfile(codeDir,'results');
    cfg = q3_optimization_config('smoke');
    opt = cfg.simOpt;
    opt.tMax = 25;
    opt.plot = false;
    opt.verbose = false;
    opt.enforceChargeCap = true;
    plannedHeatingS = [20;23];
    minTemperatureAt25C = nan(2,1);
    observedStopTimeS = nan(2,1);
    successBy25s = false(2,1);
    status = strings(2,1);
    for i = 1:2
        try
            r = pemfc_stack5_simulate_q3( ...
                'coheat',ones(1,5),plannedHeatingS(i),opt);
            observedStopTimeS(i) = r.stopTimeS;
            if r.solverTerminatedUnexpectedly
                status(i) = "numerical_failure";
            elseif r.stopTimeS < 25-1e-7 && ~r.success
                status(i) = "physical_failure_before_25s";
            else
                status(i) = "evaluated";
                minTemperatureAt25C(i) = min(r.TavgC(end,:));
                successBy25s(i) = r.success && ...
                    r.successTimeS <= 25;
            end
        catch ME
            status(i) = string(ME.identifier);
        end
        fprintf('heater off %.1f s: %s, minT25=%.5f C\n', ...
            plannedHeatingS(i),status(i),minTemperatureAt25C(i));
    end
    summary = table(plannedHeatingS,minTemperatureAt25C, ...
        observedStopTimeS,successBy25s,status);
    writetable(summary,fullfile(outDir,'q3_25s_duration_probe.csv'), ...
        'Encoding','UTF-8');
end
