function summary = q3_25s_gradient_probe()
% One-sided feasible-direction check at the all-max power boundary.
    codeDir = fileparts(mfilename('fullpath'));
    addpath(codeDir);
    outDir = fullfile(codeDir,'results');
    a = load(fullfile(outDir,'time_priority_candidates.mat'),'coProbe');
    cfg = q3_optimization_config('smoke');
    opt = cfg.simOpt;
    opt.tMax = 25;
    opt.plot = false;
    opt.verbose = false;
    opt.enforceChargeCap = true;
    idx = find(abs(a.coProbe.t-25)<1e-9,1);
    referenceMinT = min(a.coProbe.TavgC(idx,:));
    q = ones(5,5)-0.05*eye(5);
    minTemperatureC = nan(5,1);
    deltaFromAllMaxC = nan(5,1);
    observedStopTimeS = nan(5,1);
    status = strings(5,1);
    for k = 1:5
        try
            r = pemfc_stack5_simulate_q3('coheat',q(k,:),25,opt);
            observedStopTimeS(k) = r.stopTimeS;
            if r.solverTerminatedUnexpectedly
                status(k) = "numerical_failure";
            elseif r.stopTimeS < 25-1e-7 && ~r.success
                status(k) = "physical_failure_before_25s";
            else
                status(k) = "evaluated";
                minTemperatureC(k) = min(r.TavgC(end,:));
                deltaFromAllMaxC(k) = ...
                    minTemperatureC(k)-referenceMinT;
            end
        catch ME
            status(k) = string(ME.identifier);
            if strcmp(ME.identifier, ...
                    'pemfc_stack5_simulate_q3:NonfiniteDerivative')
                tightOpt = opt;
                tightOpt.RelTol = 1e-5;
                tightOpt.AbsTol = 1e-9;
                tightOpt.MaxStep = 0.025;
                try
                    r = pemfc_stack5_simulate_q3( ...
                        'coheat',q(k,:),25,tightOpt);
                    observedStopTimeS(k) = r.stopTimeS;
                    if ~r.solverTerminatedUnexpectedly && ...
                            (r.stopTimeS >= 25-1e-7 || r.success)
                        status(k) = "evaluated_tighter";
                        minTemperatureC(k) = min(r.TavgC(end,:));
                        deltaFromAllMaxC(k) = ...
                            minTemperatureC(k)-referenceMinT;
                    else
                        status(k) = "numerical_failure_tighter";
                    end
                catch ME2
                    status(k) = string(ME2.identifier);
                end
            end
        end
        fprintf('cell%d reduction: %s, minT=%.5f C, delta=%.5f C\n', ...
            k,status(k),minTemperatureC(k),deltaFromAllMaxC(k));
        save_checkpoint();
    end
    summary = current_table();

    function save_checkpoint()
        t = current_table();
        writetable(t,fullfile(outDir,'q3_25s_gradient_probe.csv'), ...
            'Encoding','UTF-8');
    end
    function t = current_table()
        t = table((1:5).',q(:,1),q(:,2),q(:,3),q(:,4),q(:,5), ...
            minTemperatureC,deltaFromAllMaxC,observedStopTimeS,status, ...
            'VariableNames',{'changedCell','q1','q2','q3','q4','q5', ...
            'minTemperatureC','deltaFromAllMaxC','observedStopTimeS', ...
            'status'});
    end
end
