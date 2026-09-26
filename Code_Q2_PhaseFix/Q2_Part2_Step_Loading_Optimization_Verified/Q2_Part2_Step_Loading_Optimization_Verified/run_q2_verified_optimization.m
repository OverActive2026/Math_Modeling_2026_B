function report=run_q2_verified_optimization(newEvaluations,stage,coldTemperatureC)
% No legacy cache import. Cold feasibility is part of every accepted point.
if nargin<1, newEvaluations=6; end
if nargin<2, stage=7; end
if nargin<3
    if stage==7, coldTemperatureC=-11.2; else, coldTemperatureC=-11; end
end
assert(ismember(stage,[6 7]) && isscalar(coldTemperatureC) && ...
    isfinite(coldTemperatureC) && coldTemperatureC<-10);
assert(newEvaluations>=3 && newEvaluations==floor(newEvaluations));
folder=fileparts(mfilename('fullpath')); addpath(folder,'-begin');
simOpt=q2_verified_options(); temps=[-10 coldTemperatureC];
timeCapS=600; % Computational budget, not a claimed physical startup bound.
sources={'q2_evaluate_schedule.m','q2_step_strategy.m','q2_step_encode.m', ...
    'q2_verified_options.m','q2_optimization_options.m', ...
    'pemfc_stack5_simulate.m','pemfc_setup_ice.m', ...
    'gas_transport_state_ice.m','thermal_temperature_state_ice.m', ...
    'water_ice_state_ice.m','q2_current_strategy.m'};
code=cellfun(@(s)fileread(fullfile(folder,s)),sources,'UniformOutput',false);
meta=struct('version',1,'simOpt',simOpt,'temps',temps,'stage',stage, ...
    'sourceContents',{code});
cache=struct('x',{},'result',{},'wallSeconds',{});
cacheFile=fullfile(folder,sprintf('corrected_cache_%d_T%g.mat',stage,coldTemperatureC));
if isfile(cacheFile)
    old=load(cacheFile,'cache','meta');
    assert(isequaln(meta,old.meta),'Q2:CacheMismatch', ...
        'Settings/code changed; archive the old cache before restarting.');
    cache=old.cache;
end
known=arrayfun(@(c)~c.result.solver_failed && ~c.result.time_censored,cache);
if any(known)
    initial=vertcat(cache(known).x);
    if size(initial,1)==1
        second=q2_step_encode(q2_reference_strategy(stage));
        second(stage)=max(0,second(stage)-.01);
        initial=[initial;second];
    end
else
    seed=q2_step_encode(q2_reference_strategy(stage));
    second=seed; second(stage)=max(0,second(stage)-.01);
    initial=[seed;second];
end
opt=struct('maxEvaluations',size(initial,1)+newEvaluations, ...
    'initialSamples',size(initial,1),'initialPoints',initial, ...
    'simOpt',simOpt,'seed',20260926+numel(cache),'evaluator',@evaluate, ...
    'initialRadius',.01,'minRadius',.002, ...
    'constraintLearning','event_aware','globalCandidateFraction',.15, ...
    'localRiskWeight',.4,'maxDistanceFromSuccess',.15,'globalScoutPeriod',6, ...
    'nInspectors',512,'nCandidates',512,'acquisitionRestarts',3, ...
    'useIncumbentCutoff',true,'cutoffSlackS',.5);
% First call evaluates the two seeds plus newEvaluations further candidates.
report=q2_optimize_trbo(sprintf('step%d_full',stage),opt);
report.simOpt=simOpt; report.temperaturesC=temps; report.stage=stage;
report.timeCapS=timeCapS; report.actualSimulationCalls=numel(cache);
report.cumulativeEvaluationWallSeconds=sum([cache.wallSeconds]);
report.note='Observed feasible best only. Quick grid; FuRBO-inspired/LogEI hybrid, not exact FuRBO or a minimum-temperature proof.';
save(fullfile(folder,sprintf('corrected_report_%d_T%g.mat',stage,coldTemperatureC)), ...
    'report','-v7.3');
fprintf('OBSERVED BEST %.8f simulated seconds; feasible=%d; censored=%d\n', ...
    report.bestTimeS,report.nSuccess,report.nTimeCensored);
    function r=evaluate(x,override)
        if nargin<2, override=simOpt; end
        % Overrides may alter only computational horizons/printing, never physics.
        allowed={'tMax','coldTMaxS','optimizationCutoffActive'};
        assert(isequaln(rm_existing(override,allowed),rm_existing(simOpt,allowed)), ...
            'Q2:PhysicsOverride','Cutoff overlay changed physics or integration settings.');
        for j=1:numel(cache)
            if isequaln(x,cache(j).x) && ~cache(j).result.solver_failed && ...
                    ~cache(j).result.time_censored
                r=cache(j).result; return
            end
        end
        clock=tic;
        r=q2_evaluate_schedule(q2_step_strategy(x),temps,override,timeCapS);
        r.x=x; wallSeconds=toc(clock);
        cache(end+1)=struct('x',x,'result',r,'wallSeconds',wallSeconds);
        save(cacheFile,'cache','meta','-v7.3');
    end
end
function s=rm_existing(s,names)
names=names(isfield(s,names)); if ~isempty(names), s=rmfield(s,names); end
end
