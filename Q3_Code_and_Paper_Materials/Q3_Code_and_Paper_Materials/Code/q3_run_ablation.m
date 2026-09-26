function summary = q3_run_ablation(cfg,seeds)
% Equal-budget baseline versus physics-guided search for both modes.
    if nargin < 1, cfg = q3_optimization_config('smoke'); end
    if nargin < 2, seeds = [20260925 20260926 20260927]; end
    if ~isnumeric(seeds) || any(~isfinite(seeds)) || ...
            any(seeds < 0) || any(mod(seeds,1) ~= 0)
        error('q3_run_ablation:Seeds','Seeds must be nonnegative integers.');
    end
    n = 4*numel(seeds);
    mode = strings(n,1);
    seed = zeros(n,1);
    guidance = false(n,1);
    evaluations = zeros(n,1);
    guidanceEvaluations = zeros(n,1);
    feasible = false(n,1);
    bestEnergyJ = nan(n,1);
    bestTimeS = nan(n,1);
    wallTimeS = nan(n,1);
    row = 0;
    for m = {'preheat','coheat'}
        for g = [false true]
            for s = seeds(:).'
                row = row+1;
                localCfg = cfg;
                localCfg.seed = s;
                localCfg.guidance = g;
                run = q3_run_optimization(m{1},localCfg);
                mode(row) = string(m{1});
                seed(row) = s;
                guidance(row) = g;
                evaluations(row) = run.evaluations;
                guidanceEvaluations(row) = run.guidanceEvaluations;
                wallTimeS(row) = run.wallTimeS;
                feasible(row) = ~isempty(run.best) && ...
                    strcmp(run.best.status,'success');
                if feasible(row)
                    bestEnergyJ(row) = run.best.energyJ;
                    bestTimeS(row) = run.best.successTimeS;
                end
            end
        end
    end
    summary = table(mode,seed,guidance,evaluations, ...
        guidanceEvaluations,feasible,bestEnergyJ,bestTimeS,wallTimeS);
    outdir = fullfile(fileparts(mfilename('fullpath')),'results');
    if ~exist(outdir,'dir'), mkdir(outdir); end
    writetable(summary,fullfile(outdir, ...
        sprintf('ablation_%s_budget%d.csv',cfg.profile,cfg.maxEvaluations)));
end
