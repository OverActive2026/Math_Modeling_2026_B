function report=run_q2_step5_full_search(maxNewEvaluations)
%RUN_Q2_STEP5_FULL_SEARCH Real-model, nine-variable five-step optimization.
% Example: report=run_q2_step5_full_search(20);
% All five current levels and all four switch times are optimized. Earlier
% five-step simulations are replayed only if physical sources and options
% match; the eight fresh perturbations move EVERY one of the nine inputs.
if nargin<1, maxNewEvaluations=20; end
assert(isscalar(maxNewEvaluations) && maxNewEvaluations>=2 && ...
    maxNewEvaluations==floor(maxNewEvaluations));
codeDir=fileparts(mfilename('fullpath'));
addpath(codeDir);
simOpt=q2_optimization_options();
simOpt.progressIntervalS=5;
physicalFiles={'pemfc_stack5_simulate.m','pemfc_setup_ice.m', ...
    'gas_transport_state_ice.m','thermal_temperature_state_ice.m', ...
    'water_ice_state_ice.m','q2_current_strategy.m'};
sourceFiles=[physicalFiles,{'q2_step5_full_strategy.m', ...
    'q2_step5_full_encode.m','q2_evaluate_step5_full_real.m'}];
sourceStamp=zeros(numel(sourceFiles),2);
for k=1:numel(sourceFiles)
    info=dir(fullfile(codeDir,sourceFiles{k}));
    assert(numel(info)==1,'Missing model source: %s',sourceFiles{k});
    sourceStamp(k,:)=[info.bytes,info.datenum];
end
cacheMeta=struct('sourceStamp',sourceStamp,'simOpt',simOpt, ...
    'initialTemperatureC',-10);
cacheFile=fullfile(codeDir,'q2_step5_full_cache.mat');
cache=struct('x',{},'result',{});
if isfile(cacheFile)
    saved=load(cacheFile,'cache','cacheMeta');
    if isequaln(saved.cacheMeta,cacheMeta)
        cache=saved.cache;
    else
        fprintf('Source/settings changed: old full-search cache ignored.\n');
    end
end
if isempty(cache)
    assert(maxNewEvaluations>=3, ...
        'Portable cold start needs at least three new evaluations.');
    bridgeFile=fullfile(codeDir,'q2_step5_bridge_report.mat');
    if isfile(bridgeFile)
        prior=load(bridgeFile,'report');
        p=prior.report;
        samePhysics=isfield(p,'sourceStamp') && ...
            size(p.sourceStamp,1)>=numel(physicalFiles) && ...
            isequaln(p.sourceStamp(1:numel(physicalFiles),:), ...
                sourceStamp(1:numel(physicalFiles),:));
        if samePhysics && isequaln(p.simOpt,simOpt)
            for k=1:numel(p.evaluations)
                result=p.evaluations{k};
                if isempty(result) || result.solver_failed || ...
                        result.time_censored, continue; end
                x=q2_step5_full_encode(result.strategy);
                result.x=x;
                cache(end+1)=struct('x',x,'result',result); %#ok<AGROW>
            end
            fprintf('Reused %d prior physical runs; no model calls for them.\n', ...
                numel(cache));
        end
    end
end
if isempty(cache)
    % A portable release may contain the published report but not its
    % source-timestamp-bound cache. Re-evaluate its incumbent locally;
    % never replay an unverified physical label after moving the files.
    reportFile=fullfile(codeDir,'q2_step5_full_report.mat');
    assert(isfile(reportFile), ...
        'No compatible cache or published step5 report is available.');
    published=load(reportFile,'report');
    assert(isfield(published.report,'bestX') && ...
        numel(published.report.bestX)==9 && ...
        all(isfinite(published.report.bestX)), ...
        'Published step5 report has no valid nine-variable incumbent.');
    initial=double(published.report.bestX(:).');
    bestX=initial;
    fprintf('Portable start: published incumbent will be re-simulated.\n');
else
    initial=vertcat(cache.x);
    [~,bestIndex]=min(arrayfun(@(c) feasible_time(c.result),cache));
    bestX=cache(bestIndex).x;
end
newDesign=min(8,maxNewEvaluations-2);
rng(20261001+numel(cache),'twister');
for k=1:newDesign
    % Each new proposal changes every coordinate, including early levels
    % and switch times that the old bridge search held fixed.
    if numel(cache)>24
        stepSize=0.008+0.055*rand(1,9);
    else
        stepSize=0.08+0.18*rand(1,9);
    end
    delta=stepSize.*sign(rand(1,9)-0.5);
    initial(end+1,:)=min(1,max(0,bestX+delta)); %#ok<AGROW>
end
newCalls=0;
boOpt=struct('maxEvaluations',numel(cache)+maxNewEvaluations, ...
    'initialSamples',size(initial,1),'initialPoints',initial, ...
    'seed',20261001+numel(cache),'initialRadius',0.65, ...
    'minRadius',0.10,'minFeasibilityProbability',0.10, ...
    'localRiskWeight',0.5,'constraintLearning','event_aware', ...
    'globalCandidateFraction',0.20,'maxDistanceFromSuccess',0.32, ...
    'globalScoutPeriod',6,'nInspectors',512, ...
    'nCandidates',512,'acquisitionRestarts',3,'evaluator',@evaluate);
report=q2_optimize_trbo('step5_full',boOpt);
report.simOpt=simOpt;
report.sourceStamp=sourceStamp;
report.sourceFiles=sourceFiles;
report.newModelCalls=newCalls;
report.note=['All 9 step parameters vary. Historical physical runs are ', ...
    'replayed only under an identical physical-source stamp. Event-aware ', ...
    'constraint regression, globally mixed candidates and analytic LogEI ', ...
    'are engineering adaptations, not a strict FuRBO reproduction.'];
save(fullfile(codeDir,'q2_step5_full_report.mat'),'report','-v7.3');
fprintf('Nine-variable best startup: %.9f s; new calls: %d\n', ...
    report.bestTimeS,newCalls);
disp(report.bestStrategy);

    function out=evaluate(x)
        x=double(x(:).');
        for i=1:numel(cache)
            if isequaln(cache(i).x,x)
                out=cache(i).result;
                fprintf('Reusing compatible real-model evaluation.\n');
                return
            end
        end
        out=q2_evaluate_step5_full_real(x,-10,simOpt);
        newCalls=newCalls+1;
        if ~out.solver_failed
            cache(end+1)=struct('x',x,'result',out); %#ok<AGROW>
            save(cacheFile,'cache','cacheMeta','-v7.3');
        end
    end
end

function t=feasible_time(result)
if result.success, t=result.startup_time; else, t=inf; end
end
