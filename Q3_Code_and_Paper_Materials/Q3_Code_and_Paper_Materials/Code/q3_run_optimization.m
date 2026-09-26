function run = q3_run_optimization(mode,cfg)
% Six-variable, feasibility-first search on the existing physical model.
% Preheat: multi-start coordinate pattern search. Coheat: DE/rand/1/bin.
    if nargin < 2, cfg = q3_optimization_config('smoke'); end
    mode = char(mode);
    if ~ismember(mode,{'preheat','coheat'})
        error('q3_run_optimization:Mode','Mode must be preheat or coheat.');
    end
    rng(cfg.seed,'twister');
    cache = containers.Map('KeyType','char','ValueType','any');
    count = 0;
    guidanceCount = 0;
    history = struct('evaluation',{},'mode',{},'z',{},'status',{}, ...
        'energyJ',{},'violation',{},'successTimeS',{},'elapsedS',{});
    best = [];
    wallTimer = tic;
    if strcmp(mode,'preheat')
        search_preheat();
    else
        search_coheat();
    end
    run = struct('mode',mode,'cfg',cfg,'best',best, ...
        'history',history,'evaluations',count, ...
        'guidanceEvaluations',guidanceCount,'wallTimeS',toc(wallTimer), ...
        'rngState',rng);
    outdir = fullfile(fileparts(mfilename('fullpath')),'results');
    if ~exist(outdir,'dir'), mkdir(outdir); end
    stem = sprintf('q3_%s_%s_%s_guide%d_seed%d', ...
        mode,cfg.profile,cfg.acceptanceMode,logical(cfg.guidance),cfg.seed);
    save(fullfile(outdir,[stem '.mat']),'run','-v7.3');
    if ~isempty(history)
        writetable(struct2table(history),fullfile(outdir,[stem '.csv']));
    end
    if isempty(best) || ~strcmp(best.status,'success')
        fprintf('%s: no feasible solution in %d evaluations.\n',mode,count);
    else
        fprintf('%s: best exploratory E=%.3f J, q=%s, th=%.4f s, t_s=%.4f s (%d evaluations).\n', ...
            mode,best.energyJ,mat2str(best.q,5),best.thRequestedS, ...
            best.successTimeS,count);
    end

    function [candidate,wasNew] = evaluate(z)
        z = min(max(double(z(:).'),0),1);
        key = [mode ':' reshape(num2hex(z).',1,[])];
        if isKey(cache,key)
            candidate = cache(key);
            wasNew = false;
            return;
        end
        if count >= cfg.maxEvaluations
            candidate = [];
            wasNew = false;
            return;
        end
        timer = tic;
        candidate = q3_evaluate_candidate(mode,z,cfg);
        elapsed = toc(timer);
        candidate.result = [];
        count = count+1;
        wasNew = true;
        cache(key) = candidate;
        history(count) = struct('evaluation',count,'mode',mode, ...
            'z',z,'status',candidate.status,'energyJ',candidate.energyJ, ...
            'violation',candidate.violation, ...
            'successTimeS',candidate.successTimeS,'elapsedS',elapsed);
        if isempty(best) || q3_compare_candidates(candidate,best)
            best = candidate;
        end
        fprintf('%s eval %d/%d: %s E=%.3f J th=%.2f s (%.1f s)\n', ...
            mode,count,cfg.maxEvaluations,candidate.status, ...
            candidate.energyJ,candidate.thRequestedS,elapsed);
    end

    function [metric,cellIndex] = probe(z,tref,cellIndex)
        metric = NaN;
        if nargin < 3, cellIndex = []; end
        if count >= cfg.maxEvaluations, return; end
        probeCfg = cfg;
        probeCfg.simOpt.stopOnSuccess = false;
        if strcmp(mode,'coheat')
            probeCfg.coheatObservationS = tref;
        end
        timer = tic;
        trial = q3_evaluate_candidate(mode,z,probeCfg);
        elapsed = toc(timer);
        count = count+1;
        guidanceCount = guidanceCount+1;
        history(count) = struct('evaluation',count,'mode',mode, ...
            'z',z,'status','probe','energyJ',NaN, ...
            'violation',NaN,'successTimeS',NaN,'elapsedS',elapsed);
        if isempty(trial.result) || ...
                trial.result.stopTimeS < tref-1e-6
            return;
        end
        T = trial.result.TavgC(end,:);
        if isempty(cellIndex)
            [~,cellIndex] = min(T);
        end
        metric = T(cellIndex)/30;
    end

    function propose_transfer(anchor)
        if ~cfg.guidance || guidanceCount+8 > ...
                floor(0.2*cfg.maxEvaluations) || ...
                count+8 > cfg.maxEvaluations
            return;
        end
        if strcmp(mode,'preheat')
            tref = anchor.thRequestedS;
        else
            tref = min(anchor.thRequestedS,anchor.successTimeS);
            if ~isfinite(tref), return; end
        end
        if tref <= 0, return; end
        [base,k] = probe(anchor.z,tref);
        if ~isfinite(base), return; end
        sensitivity = nan(1,5);
        for j = 1:5
            zPert = anchor.z;
            if zPert(j) <= 0.99
                zPert(j) = zPert(j)+0.01;
                sign = 1;
            else
                zPert(j) = zPert(j)-0.01;
                sign = -1;
            end
            [value,~] = probe(zPert,tref,k);
            if ~isfinite(value), return; end
            sensitivity(j) = (value-base)/(0.01*sign);
        end
        receiver = find(anchor.q < 0.98);
        donor = find(anchor.q > 0.02);
        if isempty(receiver) || isempty(donor), return; end
        [~,i] = max(sensitivity(receiver));
        gain = receiver(i);
        donor(donor == gain) = [];
        if isempty(donor), return; end
        [~,i] = min(sensitivity(donor));
        give = donor(i);
        if sensitivity(gain) <= sensitivity(give), return; end
        proposal = anchor.z;
        proposal(gain) = proposal(gain)+0.02;
        proposal(give) = proposal(give)-0.02;
        [candidate,wasNew] = evaluate(proposal);
        guidanceCount = guidanceCount+double(wasNew);
        if ~isempty(candidate) && q3_compare_candidates(candidate,anchor)
            return;
        end
        if count >= cfg.maxEvaluations, return; end
        proposal = anchor.z;
        proposal(gain) = proposal(gain)+0.01;
        proposal(give) = proposal(give)-0.01;
        [~,wasNew] = evaluate(proposal);
        guidanceCount = guidanceCount+double(wasNew);
    end

    function search_preheat()
        qSeeds = [0.9 0.7 0.65 0.7 0.9; 1 1 1 1 1; ...
            0.8 0.8 0.8 0.8 0.8; ...
            1 0.7 0.5 0.7 1; 0.5 0.7 1 0.7 0.5; ...
            0.2+0.8*rand(1,5)];
        startRecords = cell(size(qSeeds,1),1);
        for i = 1:size(qSeeds,1)
            if count >= cfg.maxEvaluations, break; end
            initialTh = min(cfg.timeUpperS, ...
                max(35,40/mean(qSeeds(i,:))));
            if i == 1, initialTh = min(50.3,cfg.timeUpperS); end
            if i == 2, initialTh = min(35,cfg.timeUpperS); end
            [startRecords{i},~] = evaluate([qSeeds(i,:), ...
                initialTh/cfg.timeUpperS]);
        end
        delta = 0.10;
        while count < cfg.maxEvaluations && delta >= 1e-4
            if isempty(best), break; end
            current = best;
            propose_transfer(current);
            improved = false;
            order = [6 randperm(5)];
            for ii = 1:6
                d = order(ii);
                for sign = [-1 1]
                    if count >= cfg.maxEvaluations, break; end
                    z = current.z;
                    step = delta;
                    if d == 6
                        step = delta*30/cfg.timeUpperS;
                    end
                    z(d) = min(max(z(d)+sign*step,0),1);
                    if z(d) == current.z(d), continue; end
                    [candidate,~] = evaluate(z);
                    if q3_compare_candidates(candidate,current)
                        current = candidate;
                        improved = true;
                    end
                end
            end
            if improved
                delta = min(0.20,1.2*delta);
            else
                delta = delta/2;
            end
        end
    end

    function search_coheat()
        np = min(48,max(4,floor(cfg.maxEvaluations/2)));
        P = [0.6+0.4*rand(np,5), ...
            (30+90*rand(np,1))/cfg.timeUpperS];
        P(1,:) = [ones(1,5),min(60,cfg.timeUpperS)/cfg.timeUpperS];
        P(2,:) = [[1 0.7 0.5 0.7 1], ...
            min(80,cfg.timeUpperS)/cfg.timeUpperS];
        if np >= 3
            P(3,:) = [[0.95 0.9 0.85 0.9 0.95], ...
                min(60,cfg.timeUpperS)/cfg.timeUpperS];
        end
        if np >= 4
            P(4,:) = [[1 0.9 0.9 0.9 1], ...
                min(60,cfg.timeUpperS)/cfg.timeUpperS];
        end
        if np >= 5
            P(5,:) = [[0.9 1 1 1 0.9], ...
                min(60,cfg.timeUpperS)/cfg.timeUpperS];
        end
        R = cell(np,1);
        for i = 1:np
            if count >= cfg.maxEvaluations, break; end
            [R{i},~] = evaluate(P(i,:));
        end
        generation = 0;
        while count < cfg.maxEvaluations
            generation = generation+1;
            Pold = P;
            Rold = R;
            for i = 1:np
                if count >= cfg.maxEvaluations, break; end
                others = setdiff(1:np,i);
                picked = others(randperm(numel(others),3));
                v = min(max(Pold(picked(1),:) + ...
                    0.65*(Pold(picked(2),:)-Pold(picked(3),:)),0),1);
                mask = rand(1,6) < 0.85;
                mask(randi(6)) = true;
                u = Pold(i,:);
                u(mask) = v(mask);
                [Ru,~] = evaluate(u);
                if ~isempty(Ru) && q3_compare_candidates(Ru,Rold{i})
                    P(i,:) = u;
                    R{i} = Ru;
                end
            end
            if mod(generation,5) == 0 && ~isempty(best)
                propose_transfer(best);
            end
        end
    end
end
