function report = q2_optimize_furbo_real(opt)
%Q2_OPTIMIZE_FURBO_REAL Run the inspector-TR + analytic LogEI hybrid on Q2.
% FuRBO's inspector-defined region is paired with Ament et al.'s LogEI,
% rather than FuRBO's original joint Thompson sampling.
% Prior real trajectories are reused as an initial design, while a small
% LHS extension and every subsequent infill point use the real stack model.
if nargin<1, opt=struct(); end
codeDir=fileparts(mfilename('fullpath'));
opt=with_default(opt,'initialTemperatureC',-10);
opt=with_default(opt,'simOpt',q2_optimization_options());
opt=with_default(opt,'seed',20260925);
opt=with_default(opt,'maxWarmStart',200);
opt=with_default(opt,'nLHSNew',2);
opt=with_default(opt,'maxNewEvaluations',8);
opt=with_default(opt,'initialRadius',1);
opt=with_default(opt,'minRadius',5e-8);
opt=with_default(opt,'minFeasibilityProbability',0);
% No hand-picked new candidates: subsequent points come from the
% literature-based surrogate acquisition. Existing real evaluations are
% retained as warm-start observations under the same physics/grid.
opt=with_default(opt,'seedPoints',zeros(0,3));
% Allow slower successful starts to teach the time GP. A 0.5 s window
% censored almost every non-improving run and left the GP under-informed.
opt=with_default(opt,'cutoffSlackS',12);
opt=with_default(opt,'warmStartCacheFile',fullfile(codeDir, ...
    'q2_search_continuation_checkpoint.mat'));
opt=with_default(opt,'checkpointFile',fullfile(codeDir, ...
    'q2_furbo_real_checkpoint.mat'));
assert(opt.maxWarmStart>=2 && opt.maxWarmStart==floor(opt.maxWarmStart));
assert(opt.nLHSNew>=0 && opt.nLHSNew==floor(opt.nLHSNew));
assert(opt.maxNewEvaluations>=opt.nLHSNew && ...
    opt.maxNewEvaluations==floor(opt.maxNewEvaluations));
assert(opt.cutoffSlackS>=0 && isfinite(opt.cutoffSlackS));
assert(opt.minRadius>0 && opt.minRadius<=opt.initialRadius && ...
    isfinite(opt.minRadius));
assert(opt.minFeasibilityProbability>=0 && ...
    opt.minFeasibilityProbability<1);
assert(size(opt.seedPoints,2)==3 && all(isfinite(opt.seedPoints(:))) && ...
    all(opt.seedPoints(:)>=0) && all(opt.seedPoints(:)<=1), ...
    'seedPoints must be rows in [0,1]^3.');

sourceFiles={'pemfc_stack5_simulate.m','pemfc_setup_ice.m', ...
    'gas_transport_state_ice.m','thermal_temperature_state_ice.m', ...
    'water_ice_state_ice.m','q2_current_strategy.m', ...
    'q2_decode_strategy.m','q2_evaluate_real.m'};
stamp=zeros(numel(sourceFiles),2);
for k=1:numel(sourceFiles)
    info=dir(fullfile(codeDir,sourceFiles{k}));
    if numel(info)~=1
        error('Q2:MissingSource','Missing model file: %s',sourceFiles{k});
    end
    stamp(k,:)=[info.bytes,info.datenum];
end

cache=empty_cache();
cache=append_cache(cache,opt.warmStartCacheFile);
cache=append_cache(cache,opt.checkpointFile);
[cache,nRemapped]=remap_legacy_linear_cache(cache,stamp,sourceFiles);
cache=deduplicate_cache(cache);
if nRemapped>0
    fprintf('Remapped %d legacy linear evaluations to the expanded j0 domain.\n',nRemapped);
end
resumeState=struct('radius',opt.initialRadius, ...
    'goodStreak',0,'badStreak',0);
if isfile(opt.checkpointFile)
    saved=load(opt.checkpointFile);
    if isfield(saved,'cache') && isfield(saved.cache,'trustState')
        prior=saved.cache.trustState;
        if isfield(prior,'sourceStamp') && ...
                isequaln(prior.sourceStamp,stamp) && ...
                isequaln(prior.temperatureC,opt.initialTemperatureC) && ...
                isequaln(prior.options,core_options(opt.simOpt))
            resumeState=struct('radius',prior.radius, ...
                'goodStreak',prior.goodStreak, ...
                'badStreak',prior.badStreak);
        end
    end
end
cache.trustState=struct('radius',resumeState.radius, ...
    'goodStreak',resumeState.goodStreak, ...
    'badStreak',resumeState.badStreak,'sourceStamp',stamp, ...
    'temperatureC',opt.initialTemperatureC, ...
    'options',core_options(opt.simOpt));
eligible=[];
seenX=zeros(0,3);
for k=1:numel(cache.entries)
    e=cache.entries(k);
    if numel(e.x)~=3 || any(~isfinite(e.x)) || ...
            any(e.x<0) || any(e.x>1) || ...
            ~isequaln(e.temperatureC,opt.initialTemperatureC) || ...
            ~isequaln(e.sourceStamp,stamp) || ...
            ~isequaln(core_options(e.options),core_options(opt.simOpt))
        continue
    end
    r=e.result;
    if ~(isfield(r,'solver_failed') && r.solver_failed) && ...
            ~(isfield(r,'time_censored') && r.time_censored) && ...
            ~any(all(seenX==e.x,2))
        eligible(end+1)=k; %#ok<AGROW>
        seenX(end+1,:)=e.x; %#ok<AGROW>
    end
end
% Preserve the best successes and the nearest observed failure boundary.
successIndex=[]; failureIndex=[];
for k=eligible
    if cache.entries(k).result.success
        successIndex(end+1)=k; %#ok<AGROW>
    else
        failureIndex(end+1)=k; %#ok<AGROW>
    end
end
baselineX=[0.02/0.50,(0.01-0.001)/0.049,(0.4-0.02)/(0.50-0.02)];
bestTime=inf;
bestPoint=baselineX;
if ~isempty(successIndex)
    [~,order]=sort(arrayfun(@(k)cache.entries(k).result.startup_time, ...
        successIndex));
    successIndex=successIndex(order);
    bestTime=cache.entries(successIndex(1)).result.startup_time;
    bestPoint=cache.entries(successIndex(1)).x;
end
if ~isempty(failureIndex)
    distances=arrayfun(@(k)norm(cache.entries(k).x-bestPoint), ...
        failureIndex);
    [~,order]=sort(distances);
    failureIndex=failureIndex(order);
end
nSuccessWarm=min(numel(successIndex),ceil(0.6*opt.maxWarmStart));
nFailureWarm=min(numel(failureIndex),opt.maxWarmStart-nSuccessWarm);
% Use spare slots for additional successes instead of dropping useful data.
nSuccessWarm=min(numel(successIndex),opt.maxWarmStart-nFailureWarm);
selected=[successIndex(1:nSuccessWarm),failureIndex(1:nFailureWarm)];
warmX=zeros(numel(selected),3);
for k=1:numel(selected)
    warmX(k,:)=cache.entries(selected(k)).x;
end
% Duplicate strategy coordinates add no information to the GP.
warmX=unique(warmX,'rows','stable');
nWarm=size(warmX,1);
if isempty(successIndex)
    % A changed physical model invalidates old successes. Evaluate the
    % unoptimized reference curve first, then explore with new LHS points.
    warmX=[baselineX;warmX];
end
initialPoints=unique([warmX;opt.seedPoints],'rows','stable');
nInitial=max(size(initialPoints,1)+opt.nLHSNew,2);
nTotal=nWarm+opt.maxNewEvaluations;
if nInitial>nTotal
    error('Q2:BudgetTooSmall','Increase maxNewEvaluations for initial design.');
end
fprintf('FuRBO启发式：复用%d条真实记录；最多%d次新评价。\n', ...
    nWarm,opt.maxNewEvaluations);
newCalls=0; cacheHits=0;
boOpt=struct('maxEvaluations',nTotal,'initialSamples',nInitial, ...
    'initialPoints',initialPoints,'seed',opt.seed, ...
    'evaluator',@cached_real_evaluate,'useIncumbentCutoff',false, ...
    'initialRadius',max(resumeState.radius,opt.minRadius), ...
    'minRadius',opt.minRadius, ...
    'minFeasibilityProbability',opt.minFeasibilityProbability, ...
    'initialGoodStreak',resumeState.goodStreak, ...
    'initialBadStreak',resumeState.badStreak, ...
    'stateCallback',@persist_state);
report=q2_optimize_trbo('linear',boOpt);
report.method='furbo_inspector_tr_plus_stable_analytic_logei';
report.nWarmStart=nWarm;
report.nNewModelCalls=newCalls;
report.nCacheHits=cacheHits;
report.nLegacyRemapped=nRemapped;
report.resumedTrustRadius=resumeState.radius;
report.initialTemperatureC=opt.initialTemperatureC;
report.simOpt=opt.simOpt;
report.sourceStamp=stamp;
report.sourceFiles=sourceFiles;
report.note=['FuRBO inspector-defined trust region plus stable analytic ', ...
    'LogEI times feasibility, with real-model confirmation on the same ', ...
    'grid as the unoptimized script. This hybrid uses LogEI instead of ', ...
    'FuRBO joint Thompson sampling and is not exact FuRBO.'];
save_checkpoint();

    function r=cached_real_evaluate(x)
        x=double(x(:).');
        sim=opt.simOpt;
        if isfinite(bestTime)
            baseTmax=horizon(sim);
            if bestTime+opt.cutoffSlackS<baseTmax
                sim.tMax=bestTime+opt.cutoffSlackS;
                sim.optimizationCutoffActive=true;
            end
        end
        r=[];
        for entryIndex=1:numel(cache.entries)
            e=cache.entries(entryIndex);
            if ~isequaln(e.x,x) || ...
                    ~isequaln(e.temperatureC,opt.initialTemperatureC) || ...
                    ~isequaln(e.sourceStamp,stamp) || ...
                    ~isequaln(core_options(e.options),core_options(sim))
                continue
            end
            knownSuccess=e.result.success && ...
                isfinite(e.result.startup_time);
            knownFailure=~e.result.success && ...
                ~e.result.solver_failed && ~e.result.time_censored && ...
                e.result.termination_time<=horizon(sim);
            if knownSuccess || knownFailure || isequaln(e.options,sim)
                r=e.result;
                cacheHits=cacheHits+1;
                fprintf('复用真实五片模型记录。\n');
                break
            end
        end
        if isempty(r)
            newCalls=newCalls+1;
            r=q2_evaluate_real(x,'linear',opt.initialTemperatureC,sim);
            if ~r.solver_failed
                e=struct('x',x,'temperatureC',opt.initialTemperatureC, ...
                    'options',sim,'sourceStamp',stamp,'result',r);
                cache.entries(end+1)=e;
                save_checkpoint();
            end
        end
        if r.success && isfinite(r.startup_time)
            bestTime=min(bestTime,r.startup_time);
        end
    end

    function persist_state(state)
        cache.trustState.radius=state.radius;
        cache.trustState.goodStreak=state.goodStreak;
        cache.trustState.badStreak=state.badStreak;
        save_checkpoint();
    end

    function save_checkpoint()
        pending=[opt.checkpointFile '.writing'];
        save(pending,'cache','-v7.3');
        [moved,message]=movefile(pending,opt.checkpointFile,'f');
        if ~moved
            error('Q2:CheckpointSave', ...
                'Could not save BO checkpoint: %s',message);
        end
    end
end

function [cache,nRemapped]=remap_legacy_linear_cache(cache,stamp,sourceFiles)
% Reuse a real trajectory only when the physical model and solver settings
% are unchanged and its recorded strategy proves the old j0=0.10*x(1)
% parameterization.  Do not reinterpret old normalized coordinates.
nRemapped=0;
decoderIndex=find(strcmp(sourceFiles,'q2_decode_strategy.m'),1);
for k=1:numel(cache.entries)
    e=cache.entries(k);
    if numel(e.x)~=3 || any(~isfinite(e.x)) || ...
            any(e.x<0) || any(e.x>1) || ...
            ~isequal(size(e.sourceStamp),size(stamp)) || ...
            isequaln(e.sourceStamp,stamp)
        continue
    end
    sameOtherSources=isequaln(e.sourceStamp([1:decoderIndex-1, ...
        decoderIndex+1:end],:),stamp([1:decoderIndex-1, ...
        decoderIndex+1:end],:));
    if ~sameOtherSources || ~isfield(e.result,'strategy') || ...
            ~isfield(e.result.strategy,'type') || ...
            ~strcmpi(char(e.result.strategy.type),'linear')
        continue
    end
    s=e.result.strategy;
    if ~all(isfield(s,{'initialAcm2','rampRateAcm2s', ...
            'plateauAcm2'})) || ...
            ~all(isfinite([s.initialAcm2,s.rampRateAcm2s,s.plateauAcm2]))
        continue
    end
    oldJ0=0.10*e.x(1);
    oldRate=0.001+0.049*e.x(2);
    oldPlateau=oldJ0+(0.50-oldJ0)*e.x(3);
    if max(abs([s.initialAcm2-oldJ0,s.rampRateAcm2s-oldRate, ...
            s.plateauAcm2-oldPlateau]))>1e-10
        continue
    end
    newX=[s.initialAcm2/0.50, ...
        (s.rampRateAcm2s-0.001)/0.049, ...
        (s.plateauAcm2-s.initialAcm2)/(0.50-s.initialAcm2)];
    if any(~isfinite(newX)) || any(newX<0) || any(newX>1)
        continue
    end
    decoded=q2_decode_strategy(newX,'linear');
    if max(abs([decoded.initialAcm2-s.initialAcm2, ...
            decoded.rampRateAcm2s-s.rampRateAcm2s, ...
            decoded.plateauAcm2-s.plateauAcm2]))>1e-10
        continue
    end
    cache.entries(k).x=newX;
    cache.entries(k).sourceStamp=stamp;
    cache.entries(k).result.x=newX;
    nRemapped=nRemapped+1;
end
end

function cache=empty_cache()
cache=struct('entries',struct('x',{},'temperatureC',{}, ...
    'options',{},'sourceStamp',{},'result',{}));
end

function cache=append_cache(cache,path)
if ~isfile(path), return; end
saved=load(path);
if isfield(saved,'cache')
    source=saved.cache;
elseif isfield(saved,'cacheToSave')
    source=saved.cacheToSave;
else
    return
end
if isstruct(source) && isfield(source,'entries')
    cache.entries=[cache.entries,source.entries];
end
end

function cache=deduplicate_cache(cache)
% Repeated resumptions must not multiply historical physical trajectories.
kept=cache.entries([]);
for k=1:numel(cache.entries)
    e=cache.entries(k);
    duplicate=false;
    for j=1:numel(kept)
        old=kept(j);
        if isequaln(old.x,e.x) && ...
                isequaln(old.temperatureC,e.temperatureC) && ...
                isequaln(old.sourceStamp,e.sourceStamp) && ...
                isequaln(old.options,e.options)
            duplicate=true;
            break
        end
    end
    if ~duplicate
        kept(end+1)=e; %#ok<AGROW>
    end
end
cache.entries=kept;
end

function s=with_default(s,name,value)
if ~isfield(s,name) || isempty(s.(name)), s.(name)=value; end
end

function options=core_options(options)
remove={'tMax','optimizationCutoffActive'};
for k=1:numel(remove)
    if isfield(options,remove{k})
        options=rmfield(options,remove{k});
    end
end
end

function tMax=horizon(options)
tMax=180;
if isfield(options,'tMax') && ~isempty(options.tMax)
    tMax=options.tMax;
end
end
