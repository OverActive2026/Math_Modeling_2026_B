function retry = q3_retry_energy_failures(mode)
% Re-evaluate solver failures with tighter settings on the same spatial grid.
    mode = char(mode);
    outDir = fullfile(fileparts(mfilename('fullpath')),'results');
    path = fullfile(outDir,['q3_min_energy_' mode '.mat']);
    stored = load(path,'run');
    run = stored.run;
    failed = find(strcmp({run.history.status},'numerical_failure') | ...
        (strcmp({run.history.status},'infeasible') & ...
        isinf([run.history.energyJ])));
    retryPath = fullfile(outDir,['q3_min_energy_' mode '_retry.mat']);
    if exist(retryPath,'file')
        prior = load(retryPath,'retry');
        rows = prior.retry;
        already = [rows.sourceEvaluation];
        failed = setdiff(failed,already);
    else
        rows = struct('sourceEvaluation',{},'q',{},'status',{}, ...
            'energyJ',{},'startupTimeS',{},'reason',{});
    end
    cfg = q3_optimization_config('smoke');
    opt = cfg.simOpt;
    opt.RelTol = 1e-5;
    opt.AbsTol = 1e-9;
    opt.MaxStep = 0.025;
    opt.plot = false;
    opt.verbose = false;
    for index = failed
        q = run.history(index).q;
        horizon = run.history(index).horizonS;
        opt.tMax = horizon;
        status = 'numerical_failure';
        E = NaN; ts = NaN; reason = '';
        try
            r = pemfc_stack5_simulate_q3(mode,q,horizon,opt);
            if strcmp(mode,'preheat')
                crossings = r.event.time(r.event.index == 1);
                if isempty(crossings)
                    status = 'infeasible';
                else
                    th = crossings(1)+0.01;
                    opt.tMax = th;
                    r = pemfc_stack5_simulate_q3(mode,q,th,opt);
                    if r.success, status = 'success';
                    else, status = 'infeasible';
                    end
                end
            elseif r.solverTerminatedUnexpectedly
                status = 'numerical_failure';
            elseif r.success
                status = 'success';
                th = horizon;
            else
                status = 'infeasible';
            end
            E = r.totalAuxiliaryEnergyJ;
            ts = r.successTimeS;
            reason = char(r.stopReason);
            if r.success && E < run.best.energyJ-1e-6
                run.best = struct('q',q,'thRequestedS',th, ...
                    'energyJ',E,'startupTimeS',ts,'result',r);
                save(path,'run','-v7.3');
            end
        catch ME
            reason = [ME.identifier ': ' ME.message];
        end
        rows(end+1) = struct('sourceEvaluation',index,'q',q, ...
            'status',status,'energyJ',E,'startupTimeS',ts, ...
            'reason',reason);
        fprintf('%s retry %d: %s E=%.3f J\n',mode,index,status,E);
        retry = rows;
        save(retryPath,'retry');
    end
    retry = rows;
    save(retryPath,'retry');
end
