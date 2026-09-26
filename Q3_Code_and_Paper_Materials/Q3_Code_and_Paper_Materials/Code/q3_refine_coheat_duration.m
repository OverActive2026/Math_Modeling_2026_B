function logTable = q3_refine_coheat_duration()
% Test whether switching the constant coheater off before first success
% lowers heater energy while retaining the original Q3 startup constraints.
    outDir = fullfile(fileparts(mfilename('fullpath')),'results');
    archivePath = fullfile(outDir,'q3_min_energy_coheat.mat');
    stored = load(archivePath,'run');
    run = stored.run;
    best = run.best;
    cfg = q3_optimization_config('smoke');
    opt = cfg.simOpt;
    opt.tMax = 65;
    opt.plot = false;
    opt.verbose = false;
    offsets = [0.02];
    q = best.q;
    thBase = best.startupTimeS;
    rows = struct('thRequestedS',{},'q',{},'status',{}, ...
        'energyJ',{},'startupTimeS',{},'minimumVoltageV',{}, ...
        'maximumIceVolumeFraction',{},'reason',{},'elapsedS',{});
    for i = 1:numel(offsets)
        th = thBase-offsets(i);
        tick = tic;
        status = 'numerical_failure';
        E = NaN; ts = NaN; v = NaN; ice = NaN; reason = '';
        try
            r = pemfc_stack5_simulate_q3('coheat',q,th,opt);
            E = r.totalAuxiliaryEnergyJ;
            ts = r.successTimeS;
            v = r.minimumCellVoltageV;
            ice = r.maximumIceVolumeFraction;
            reason = char(r.stopReason);
            if r.success
                status = 'success';
                if E < run.best.energyJ-1e-6
                    run.best = struct('q',q,'thRequestedS',th, ...
                        'energyJ',E,'startupTimeS',ts,'result',r);
                    save(archivePath,'run','-v7.3');
                end
            else
                status = 'infeasible';
            end
        catch ME
            reason = [ME.identifier ': ' ME.message];
        end
        rows(end+1) = struct('thRequestedS',th,'q',q, ...
            'status',status,'energyJ',E,'startupTimeS',ts, ...
            'minimumVoltageV',v,'maximumIceVolumeFraction',ice, ...
            'reason',reason,'elapsedS',toc(tick));
        fprintf('coheat early stop %.3f s: %s E=%.3f J t=%.3f s\n', ...
            th,status,E,ts);
    end
    qMatrix = vertcat(rows.q);
    logTable = table([rows.thRequestedS].',qMatrix(:,1), ...
        qMatrix(:,2),qMatrix(:,3),qMatrix(:,4),qMatrix(:,5), ...
        string({rows.status}).',[rows.energyJ].', ...
        [rows.startupTimeS].',[rows.minimumVoltageV].', ...
        [rows.maximumIceVolumeFraction].', ...
        string({rows.reason}).',[rows.elapsedS].', ...
        'VariableNames',{'thRequestedS','q1','q2','q3','q4','q5', ...
        'status','energyJ','startupTimeS','minimumVoltageV', ...
        'maximumIceVolumeFraction','reason','elapsedS'});
    writetable(logTable,fullfile(outDir,'q3_coheat_duration_refine.csv'), ...
        'Encoding','UTF-8');
end
