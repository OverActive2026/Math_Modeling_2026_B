function run = q3_refine_min_energy(mode,maxNewEvaluations)
% Refine minimum heater energy on the existing Q3 physical model.
% Search first in the symmetric 3-D subspace and then in all five powers.
% Every accepted incumbent is independently replayed at its heater stop.
    if nargin < 2, maxNewEvaluations = 30; end
    mode = char(mode);
    assert(ismember(mode,{'preheat','coheat'}));
    cfg = q3_optimization_config('smoke');
    cfg.simOpt.plot = false;
    cfg.simOpt.verbose = false;
    cfg.simOpt.enforceChargeCap = true;
    outDir = fullfile(fileparts(mfilename('fullpath')),'results');
    if ~exist(outDir,'dir'), mkdir(outDir); end
    archivePath = fullfile(outDir,['q3_min_energy_' mode '.mat']);
    if exist(archivePath,'file')
        old = load(archivePath,'run');
        run = old.run;
        assert(strcmp(run.mode,mode));
    else
        run = struct('mode',mode,'best',[],'history', ...
            struct('q',{},'status',{},'energyJ',{},'startupTimeS',{}, ...
            'horizonS',{},'elapsedS',{},'reason',{}));
    end
    baselinePath = fullfile(outDir,'time_priority_candidates.mat');
    if exist(baselinePath,'file')
        if strcmp(mode,'preheat')
            baseline = load(baselinePath,'pre','preHeatingDurationS');
            baseResult = baseline.pre;
            baselineRequestedS = baseline.preHeatingDurationS;
        else
            baseline = load(baselinePath,'co','coHeatingDurationS');
            baseResult = baseline.co;
            baselineRequestedS = baseline.coHeatingDurationS;
        end
        if baseResult.success && ...
                (isempty(run.best) || ...
                baseResult.totalAuxiliaryEnergyJ < run.best.energyJ)
            run.best = struct('q',ones(1,5), ...
                'thRequestedS',baselineRequestedS, ...
                'energyJ',baseResult.totalAuxiliaryEnergyJ, ...
                'startupTimeS',baseResult.successTimeS, ...
                'result',baseResult);
        elseif ~isempty(run.best) && isequal(run.best.q,ones(1,5)) && ...
                abs(baseResult.totalAuxiliaryEnergyJ-run.best.energyJ) < 1e-5
            run.best.thRequestedS = baselineRequestedS;
        end
    end
    cache = containers.Map('KeyType','char','ValueType','logical');
    for i = 1:numel(run.history)
        cache(mat2str(run.history(i).q,12)) = true;
    end
    initialCount = numel(run.history);
    evaluate(ones(1,5));
    % Physically symmetric stack: use 3-D search to locate promising basins.
    for step = [0.20 0.10 0.05 0.02]
        improved = true;
        while improved && remaining() > 0
            improved = false;
            anchor = incumbent_q();
            for group = 1:3
                for direction = [-1 1]
                    q = anchor;
                    if group == 1, cells = [1 5];
                    elseif group == 2, cells = [2 4];
                    else, cells = 3;
                    end
                    q(cells) = min(1,max(0,q(cells)+direction*step));
                    priorEnergy = incumbent_energy();
                    evaluate(q);
                    if incumbent_energy() < priorEnergy-1e-6
                        improved = true;
                    end
                    if remaining() <= 0, break; end
                end
                if remaining() <= 0, break; end
            end
        end
    end
    % Lift the symmetry restriction: each heater may vary independently.
    for step = [0.05 0.02 0.01]
        improved = true;
        while improved && remaining() > 0
            improved = false;
            anchor = incumbent_q();
            for k = 1:5
                for direction = [-1 1]
                    q = anchor;
                    q(k) = min(1,max(0,q(k)+direction*step));
                    priorEnergy = incumbent_energy();
                    evaluate(q);
                    if incumbent_energy() < priorEnergy-1e-6
                        improved = true;
                    end
                    if remaining() <= 0, break; end
                end
                if remaining() <= 0, break; end
            end
        end
    end
    checkpoint();

    function n = remaining()
        n = maxNewEvaluations-(numel(run.history)-initialCount);
    end
    function q = incumbent_q()
        if isempty(run.best), q = ones(1,5);
        else, q = run.best.q;
        end
    end
    function e = incumbent_energy()
        if isempty(run.best), e = inf;
        else, e = run.best.energyJ;
        end
    end
    function evaluate(q)
        q = double(q(:).');
        key = mat2str(q,12);
        if isKey(cache,key) || remaining() <= 0, return; end
        cache(key) = true;
        ticEval = tic;
        energyBound = incumbent_energy();
        if isfinite(energyBound)
            horizon = min(100,max(30,energyBound/(25*sum(q))+2));
        else
            horizon = 70;
        end
        status = 'infeasible';
        reason = '';
        E = inf;
        tStart = NaN;
        try
            opt = cfg.simOpt;
            opt.tMax = horizon;
            if strcmp(mode,'preheat')
                probe = pemfc_stack5_simulate_q3(mode,q,horizon,opt);
                crossings = probe.event.time(probe.event.index == 1);
                if ~isempty(crossings) && ...
                        ~probe.solverTerminatedUnexpectedly
                    th = crossings(1)+1e-5;
                    Eestimate = 25*sum(q)*th;
                    if Eestimate < incumbent_energy()-1e-6
                        opt.tMax = th;
                        r = pemfc_stack5_simulate_q3(mode,q,th,opt);
                        if ~r.success && ~r.solverTerminatedUnexpectedly && ...
                                min(r.TavgC(end,:))>0 && ...
                                min(r.TavgC(end,:))<=opt.temperatureEventMarginK
                            th=th+1e-3;
                            opt.tMax=th;
                            r=pemfc_stack5_simulate_q3(mode,q,th,opt);
                        end
                        if r.success
                            E = r.totalAuxiliaryEnergyJ;
                            tStart = r.successTimeS;
                            status = 'success';
                            consider(r,q,th);
                        else
                            reason = char(r.stopReason);
                        end
                    else
                        status = 'bounded';
                        E = Eestimate;
                        tStart = th;
                    end
                else
                    reason = char(probe.stopReason);
                end
            else
                r = pemfc_stack5_simulate_q3(mode,q,horizon,opt);
                if r.solverTerminatedUnexpectedly
                    status = 'numerical_failure';
                    reason = char(r.stopReason);
                elseif r.success
                    E = r.totalAuxiliaryEnergyJ;
                    tStart = r.successTimeS;
                    status = 'success';
                    consider(r,q,horizon);
                else
                    reason = char(r.stopReason);
                end
            end
        catch ME
            status = 'numerical_failure';
            reason = [ME.identifier ': ' ME.message];
        end
        run.history(end+1) = struct('q',q,'status',status, ...
            'energyJ',E,'startupTimeS',tStart,'horizonS',horizon, ...
            'elapsedS',toc(ticEval),'reason',reason);
        fprintf('%s eval %d: q=%s E=%.3f J t=%.4f s %s\n', ...
            mode,numel(run.history),mat2str(q,4),E,tStart,status);
        checkpoint();
    end
    function consider(r,q,th)
        E = r.totalAuxiliaryEnergyJ;
        if E < incumbent_energy()-1e-6
            run.best = struct('q',q,'thRequestedS',th, ...
                'energyJ',E,'startupTimeS',r.successTimeS, ...
                'result',r);
        end
    end
    function checkpoint()
        save(archivePath,'run','-v7.3');
        if ~isempty(run.history)
            qMatrix = vertcat(run.history.q);
            h = table(qMatrix(:,1),qMatrix(:,2),qMatrix(:,3), ...
                qMatrix(:,4),qMatrix(:,5), ...
                string({run.history.status}).', ...
                [run.history.energyJ].', ...
                [run.history.startupTimeS].', ...
                [run.history.horizonS].', ...
                [run.history.elapsedS].', ...
                string({run.history.reason}).', ...
                'VariableNames',{'q1','q2','q3','q4','q5', ...
                'status','energyJ','startupTimeS','horizonS', ...
                'elapsedS','reason'});
            writetable(h,fullfile(outDir,['q3_min_energy_' mode '_history.csv']));
        end
    end
end
