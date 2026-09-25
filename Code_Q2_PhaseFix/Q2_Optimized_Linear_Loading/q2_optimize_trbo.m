function report = q2_optimize_trbo(kind,opt)
%Q2_OPTIMIZE_TRBO Feasibility-guided inspector TR with analytic LogEI.
% FuRBO supplies the inspector-defined trust region; Ament et al. supply
% stable analytic LogEI for ranking points inside it. This hybrid replaces
% FuRBO's joint Thompson-sampling acquisition and is NOT exact FuRBO.
% Every reported best point has been evaluated by the supplied oracle.
% Usage: report=q2_optimize_trbo('linear',struct('maxEvaluations',20));
if nargin<2, opt=struct(); end
kind = lower(char(kind));
switch kind
    case 'constant', d=3;
    case 'strict_constant', d=1;
    case 'linear', d=3;
    case 'step', d=7;
    otherwise, error('Q2:UnknownStrategy','Unknown strategy: %s',kind);
end
opt = default(opt,'maxEvaluations',24);
opt = default(opt,'initialSamples',min(12,opt.maxEvaluations));
opt = default(opt,'seed',20260924);
opt = default(opt,'initialTemperatureC',-10);
opt = default(opt,'simOpt',struct());
opt = default(opt,'evaluator',[]);
opt = default(opt,'initialPoints',zeros(0,d));
isRealEvaluator = isempty(opt.evaluator);
opt = default(opt,'useIncumbentCutoff',isRealEvaluator);
opt = default(opt,'cutoffSlackS',0.5);
opt = default(opt,'initialRadius',1);
opt = default(opt,'minRadius',5e-8);
opt = default(opt,'minFeasibilityProbability',0);
opt = default(opt,'initialGoodStreak',0);
opt = default(opt,'initialBadStreak',0);
opt = default(opt,'stateCallback',[]);
opt = default(opt,'nInspectors',1024);
opt = default(opt,'nCandidates',1024);
opt = default(opt,'acquisitionRestarts',4);
opt = default(opt,'inspectorTopFraction',0.10);
if isRealEvaluator
    opt.evaluator = @(x) q2_evaluate_real(x,kind, ...
        opt.initialTemperatureC,opt.simOpt);
end
assert(isa(opt.evaluator,'function_handle'),'evaluator must be a function handle.');
assert(opt.maxEvaluations>=2 && opt.initialSamples>=2 && ...
    opt.initialSamples<=opt.maxEvaluations,'Invalid evaluation budget.');
assert(isscalar(opt.cutoffSlackS) && isfinite(opt.cutoffSlackS) && ...
    opt.cutoffSlackS>=0,'cutoffSlackS must be nonnegative.');
assert(opt.initialRadius>0 && isfinite(opt.initialRadius));
assert(opt.minRadius>0 && opt.minRadius<=opt.initialRadius && ...
    isfinite(opt.minRadius),'Invalid minimum trust-region radius.');
assert(opt.minFeasibilityProbability>=0 && ...
    opt.minFeasibilityProbability<1 && ...
    isfinite(opt.minFeasibilityProbability));
assert(opt.initialGoodStreak>=0 && opt.initialGoodStreak==floor(opt.initialGoodStreak));
assert(opt.initialBadStreak>=0 && opt.initialBadStreak==floor(opt.initialBadStreak));
assert(isempty(opt.stateCallback) || isa(opt.stateCallback,'function_handle'));
assert(opt.nInspectors>=20 && opt.nInspectors==floor(opt.nInspectors));
assert(opt.nCandidates>=20 && opt.nCandidates==floor(opt.nCandidates));
assert(opt.acquisitionRestarts>=0 && ...
    opt.acquisitionRestarts==floor(opt.acquisitionRestarts));
assert(opt.inspectorTopFraction>0 && opt.inspectorTopFraction<=1);
assert(size(opt.initialPoints,2)==d && ...
    size(opt.initialPoints,1)<=opt.initialSamples && ...
    all(isfinite(opt.initialPoints(:))) && ...
    all(opt.initialPoints(:)>=0) && all(opt.initialPoints(:)<=1), ...
    'initialPoints must be rows in [0,1]^d.');
assert(exist('lhsdesign','file')==2 && exist('fitrgp','file')==2, ...
    'Statistics and Machine Learning Toolbox is required.');
rng(opt.seed,'twister');
N = opt.maxEvaluations;
X = nan(N,d); Y = nan(N,1); G = nan(N,1);
wallSeconds = nan(N,1);
success = false(N,1); numericalFailure = false(N,1);
timeCensored = false(N,1);
records = cell(N,1); n = 0;
acquisitionTrace=struct('x',{},'radius',{},'predictedTime',{}, ...
    'timeSD',{},'feasibilityProbability',{},'logScore',{}, ...
    'poolLogScore',{},'actualTime',{},'success',{}, ...
    'trLower',{},'trUpper',{},'center',{});
initial = lhsdesign(opt.initialSamples,d,'criterion','maximin','iterations',5);
initial(1:size(opt.initialPoints,1),:)=opt.initialPoints;
for i=1:opt.initialSamples
    evaluate(initial(i,:));
end
radius=max(opt.initialRadius,opt.minRadius);
goodStreak=opt.initialGoodStreak;
badStreak=opt.initialBadStreak;
while n<N
    valid=find(~numericalFailure(1:n) & ~timeCensored(1:n) & ...
        isfinite(G(1:n)));
    % In d dimensions, start constraint learning as soon as d+1 distinct
    % physical labels exist. Waiting longer can repeat a known failure.
    % Event-driven solvers stop precisely at a constraint boundary, so the
    % terminal voltage margin of a failed run is often zero. Learn the
    % binary feasibility label rather than treating that censored margin
    % as a continuous violation of known magnitude.
    gModel=fit_gp(X(valid,:),0.5-double(success(valid)),max(3,d+1));
    feasible=find(success(1:n) & isfinite(Y(1:n)));
    % Two distinct successful starts are already informative about time.
    % Requiring five made the search feasibility-only for many expensive
    % evaluations, contrary to the intended EI x feasibility acquisition.
    yModel=fit_gp(X(feasible,:),Y(feasible),2);
    if isempty(valid)
        center=0.5*ones(1,d);
    elseif isempty(feasible)
        [~,j]=min(G(valid));
        center=X(valid(j),:);
    else
        [~,j]=min(Y(feasible));
        center=X(feasible(j),:);
    end
    % FuRBO Algorithm 1, lines 4-7: inspectors in a ball around the best
    % observed point, then the top predicted points define a hyperrectangle.
    inspectors=sample_inspectors(center,radius,opt.nInspectors,d);
    [iGMu,iGSd]=gp_predict(gModel,inspectors);
    iLogP=log_normal_cdf(-iGMu./max(iGSd,1e-12));
    iFeasible=iLogP>=log(0.5);
    if ~isempty(yModel)
        [iYMu,~]=gp_predict(yModel,inspectors);
        iRank=iYMu;
        iRank(~iFeasible)=iGMu(~iFeasible);
    else
        iRank=-iLogP;
    end
    [~,iOrder]=sortrows([~iFeasible,iRank],[1 2]);
    top=inspectors(iOrder(1:max(2,ceil(opt.inspectorTopFraction* ...
        opt.nInspectors))),:);
    if ~isempty(feasible)
        % Hybrid safeguard: the observed feasible incumbent remains inside
        % the inspector box even when surrogate rankings are inaccurate.
        top=[top;center]; %#ok<AGROW>
    end
    trLower=max(0,min(top,[],1));
    trUpper=min(1,max(top,[],1));
    % Avoid a degenerate box when the inspectors collapse near a boundary.
    width=max(trUpper-trLower,min(0.02,0.1*radius));
    midpoint=(trLower+trUpper)/2;
    trLower=max(0,midpoint-width/2);
    trUpper=min(1,midpoint+width/2);
    pool=trLower+(trUpper-trLower).*lhsdesign(opt.nCandidates,d, ...
        'criterion','maximin','iterations',2);
    pool(end,:)=min(trUpper,max(trLower,center));
    bestObserved=inf;
    if ~isempty(feasible), bestObserved=min(Y(feasible)); end
    [pMu,pSd]=gp_predict(gModel,pool);
    poolP=exp(log_normal_cdf(-pMu./max(pSd,1e-12)));
    minP=opt.minFeasibilityProbability;
    if max(poolP)<minP
        minP=max(0,max(poolP)-0.02);
    end
    score=acquisition_score(pool,gModel,yModel,bestObserved,minP);
    poolBest=max(score);
    % Optimize the cheap acquisition itself, never the physical model.
    % Retain all pool points so refinement cannot lower its best score.
    [~,starts]=sort(score,'descend');
    for localIndex=1:min(opt.acquisitionRestarts,numel(starts))
        u=(pool(starts(localIndex),:)-trLower)./max(trUpper-trLower,eps);
        u=min(1-1e-6,max(1e-6,u));
        z0=log(u./(1-u));
        transform=@(z) trLower+(trUpper-trLower)./(1+exp(-max(-40,min(40,z))));
        objective=@(z) -acquisition_score(transform(z),gModel, ...
            yModel,bestObserved,minP);
        z=fminsearch(objective,z0,optimset('Display','off', ...
            'MaxIter',60,'MaxFunEvals',140,'TolX',1e-4,'TolFun',1e-5));
        pool(end+1,:)=transform(z); %#ok<AGROW>
    end
    score=acquisition_score(pool,gModel,yModel,bestObserved,minP);
    [~,order]=sort(score,'descend');
    chosen=[];
    requiredSeparation=min(0.002,0.02*radius);
    for pass=1:4
        for k=1:numel(order)
            c=pool(order(k),:);
            if all(vecnorm(X(1:n,:)-c,2,2)>requiredSeparation)
                chosen=c; break
            end
        end
        if ~isempty(chosen), break; end
        requiredSeparation=requiredSeparation/2;
    end
    if isempty(chosen)
        % Preserve the trust-region semantics even if all proposed points
        % are close to previous evaluations.
        chosen=pool(order(1),:);
    end
    oldBest=inf;
    if ~isempty(feasible), oldBest=min(Y(feasible)); end
    [chosenMu,chosenSd]=gp_predict(yModel,chosen);
    [chosenG,chosenGs]=gp_predict(gModel,chosen);
    chosenP=exp(log_normal_cdf(-chosenG./max(chosenGs,1e-12)));
    chosenScore=acquisition_score(chosen,gModel,yModel,bestObserved,minP);
    fprintf('Acquisition: R=%.5g predicted=%.4f s SD=%.4f P(feasible)=%.3f logCEI=%.4g\n', ...
        radius,chosenMu,chosenSd,chosenP,chosenScore);
    idx=evaluate(chosen);
    acquisitionTrace(end+1)=struct('x',chosen,'radius',radius, ...
        'predictedTime',chosenMu,'timeSD',chosenSd, ...
        'feasibilityProbability',chosenP,'logScore',chosenScore, ...
        'poolLogScore',poolBest,'actualTime',Y(idx),'success',success(idx), ...
        'trLower',trLower,'trUpper',trUpper,'center',center); %#ok<AGROW>
    if numericalFailure(idx)
        % Numerical failures consume budget but are not physical labels.
        if ~isempty(opt.stateCallback)
            opt.stateCallback(struct('radius',radius, ...
                'goodStreak',goodStreak,'badStreak',badStreak));
        end
        continue
    elseif success(idx) && Y(idx)<oldBest-1e-8
        goodStreak=goodStreak+1; badStreak=0;
    else
        badStreak=badStreak+1; goodStreak=0;
    end
    if goodStreak>=2
        radius=min(2,2*radius); goodStreak=0; badStreak=0;
    elseif badStreak>=3
        radius=max(opt.minRadius,0.5*radius); badStreak=0; goodStreak=0;
    end
    if ~isempty(opt.stateCallback)
        opt.stateCallback(struct('radius',radius, ...
            'goodStreak',goodStreak,'badStreak',badStreak));
    end
end
feasible=find(success & isfinite(Y));
report.kind=kind;
report.seed=opt.seed;
report.nEvaluations=n;
report.nSuccess=numel(feasible);
report.nNumericalFailures=sum(numericalFailure);
report.nTimeCensored=sum(timeCensored);
report.X=X; report.startupTimeS=Y;
report.signedConstraint=G; report.success=success;
report.numericalFailure=numericalFailure;
report.timeCensored=timeCensored;
report.wallSeconds=wallSeconds;
report.trustRadiusFinal=radius;
report.goodStreakFinal=goodStreak;
report.badStreakFinal=badStreak;
report.acquisition='analytic_logei_times_feasibility_in_inspector_TR';
report.acquisitionTrace=acquisitionTrace;
report.algorithmSettings=struct('nInspectors',opt.nInspectors, ...
    'nCandidates',opt.nCandidates, ...
    'acquisitionRestarts',opt.acquisitionRestarts, ...
    'minRadius',opt.minRadius, ...
    'minFeasibilityProbability',opt.minFeasibilityProbability, ...
    'inspectorTopFraction',opt.inspectorTopFraction, ...
    'anchorObservedIncumbent',true, ...
    'nominalCandidateSeparation',0.002, ...
    'uncertainty','latent_function_standard_deviation');
report.evaluations=records;
report.best=[];
report.bestX=[];
report.bestStrategy=[];
report.bestTimeS=NaN;
if ~isempty(feasible)
    [report.bestTimeS,j]=min(Y(feasible));
    idx=feasible(j);
    report.best=records{idx};
    report.bestX=X(idx,:);
    if isfield(report.best,'strategy')
        report.bestStrategy=report.best.strategy;
    end
end
report.note=['The incumbent is an observed feasible evaluation. ', ...
    'A mock evaluator does not validate the physical model.'];

    function idx=evaluate(x)
        fprintf('\n[%s %d/%d] 开始评价，归一化参数 = [', ...
            kind,n+1,N);
        fprintf(' %.4f',x);
        fprintf(' ]\n');
        evaluationClock=tic;
        if isRealEvaluator && opt.useIncumbentCutoff
            incumbent=Y(success(1:n) & isfinite(Y(1:n)));
            thisSimOpt=opt.simOpt;
            baseTmax=180;
            if isfield(thisSimOpt,'tMax') && ~isempty(thisSimOpt.tMax)
                baseTmax=thisSimOpt.tMax;
            end
            if ~isempty(incumbent)
                candidateTmax=min(baseTmax,min(incumbent)+opt.cutoffSlackS);
                if candidateTmax<baseTmax
                    thisSimOpt.tMax=candidateTmax;
                    thisSimOpt.optimizationCutoffActive=true;
                    fprintf('已知最快 %.3f s；本方案若 %.3f s 前未启动则提前结束。\n', ...
                        min(incumbent),candidateTmax);
                end
            end
            result=q2_evaluate_real(x,kind, ...
                opt.initialTemperatureC,thisSimOpt);
        else
            result=opt.evaluator(x);
        end
        evaluationWallSeconds=toc(evaluationClock);
        assert(isstruct(result) && isfield(result,'success') && ...
            isfield(result,'startup_time') && isfield(result,'V_min') && ...
            isfield(result,'q_use') && isfield(result,'ice_peak') && ...
            isfield(result,'T_min_end'), ...
            'Evaluator is missing required fields.');
        n=n+1; idx=n; X(idx,:)=x;
        wallSeconds(idx)=evaluationWallSeconds;
        success(idx)=logical(result.success) && ...
            isfinite(result.startup_time);
        if success(idx)
            Y(idx)=result.startup_time;
        else
            assert(isnan(result.startup_time), ...
                'A failed run must have startup_time=NaN.');
        end
        if isfield(result,'solver_failed') && logical(result.solver_failed)
            assert(~success(idx),'A solver failure cannot be successful.');
            numericalFailure(idx)=true;
            records{idx}=result;
            fprintf('%s eval %d/%d NUMERICAL_FAILURE (excluded from GP), wall=%.1f s\n', ...
                kind,n,N,wallSeconds(idx));
            return
        end
        if isfield(result,'time_censored') && logical(result.time_censored)
            assert(~success(idx),'A time-censored run cannot be successful.');
            timeCensored(idx)=true;
            records{idx}=result;
            fprintf('%s eval %d/%d NO_IMPROVEMENT_BY_CUTOFF (excluded from GP), wall=%.1f s\n', ...
                kind,n,N,wallSeconds(idx));
            return
        end
        iceForScore=result.ice_peak;
        if success(idx) && isfield(result,'ice_start')
            iceForScore=result.ice_start;
        end
        margins=[(0.30-result.V_min)/0.30, ...
            (result.q_use-20)/20, ...
            (iceForScore-0.99)/0.99];
        if isfield(result,'raw') && isfield(result.raw, ...
                'maximumPoreOccupancy')
            margins(end+1)=(result.raw.maximumPoreOccupancy-0.98)/0.98;
        end
        hasPhysicalEvent=isfield(result,'raw') && ...
            isfield(result.raw,'event') && ...
            isfield(result.raw.event,'index') && ...
            ~isempty(result.raw.event.index);
        if ~success(idx) && ~hasPhysicalEvent
            % Only a timeout has a meaningful distance-to-warmup signal.
            % A voltage-triggered run ends early by design: its remaining
            % temperature deficit is NOT an extra voltage violation.
            margins(end+1)=(273.15-result.T_min_end)/10;
        end
        assert(all(isfinite(margins)),'Non-finite diagnostic returned.');
        G(idx)=max(margins);
        % Success is authoritative; unsuccessful warm-up/solver termination
        % must not be learned as feasible just because limits were not hit.
        if success(idx), G(idx)=min(G(idx),-1e-5);
        else, G(idx)=max(G(idx),1e-5); end
        records{idx}=result;
        fprintf('%s eval %d/%d success=%d time=%.4g s G=%.3g wall=%.1f s\n', ...
            kind,n,N,success(idx),Y(idx),G(idx),wallSeconds(idx));
    end
end

function s=default(s,name,value)
if ~isfield(s,name) || isempty(s.(name)), s.(name)=value; end
end

function model=fit_gp(X,y,minPoints)
if size(X,1)<minPoints || numel(unique(y))<2
    model=[];
else
    model=fitrgp(X,y,'KernelFunction','ardsquaredexponential', ...
        'Standardize',true,'FitMethod','exact','PredictMethod','exact');
end
end

function [mu,sd]=gp_predict(model,X)
if isempty(model)
    mu=zeros(size(X,1),1); sd=ones(size(X,1),1);
else
    [mu,sd]=predict(model,X);
    % MATLAB predict returns response variance = latent variance + noise.
    % EI concerns the unknown objective, not another noisy observation.
    mu=double(mu(:));
    sd=sqrt(max(0,double(sd(:)).^2-double(model.Sigma)^2));
end
end

function score=acquisition_score(X,gModel,yModel,bestObserved,minP)
[gMu,gSd]=gp_predict(gModel,X);
logP=log_normal_cdf(-gMu./max(gSd,1e-12));
if isempty(yModel) || ~isfinite(bestObserved)
    score=logP+0.05*log(max(gSd,1e-12));
else
    [mu,sd]=gp_predict(yModel,X);
    score=q2_logei(bestObserved,mu,sd)+logP;
end
score=max(-1e100,score);
score(exp(logP)<minP)=-1e100;
end

function X=sample_inspectors(center,radius,n,d)
% Rejection sampling is uniform on ball intersected with [0,1]^d.
% Clipping sphere samples would put spurious probability mass on faces.
lo=max(0,center-radius); hi=min(1,center+radius);
X=zeros(n,d); filled=0;
while filled<n
    trial=lo+(hi-lo).*rand(max(256,4*(n-filled)),d);
    trial=trial(sum((trial-center).^2,2)<=radius^2,:);
    count=min(n-filled,size(trial,1));
    X(filled+(1:count),:)=trial(1:count,:);
    filled=filled+count;
end
end

function y=log_normal_cdf(x)
y=zeros(size(x));
negative=x<0;
y(~negative)=log(0.5*erfc(-x(~negative)/sqrt(2)));
xn=x(negative);
y(negative)=-0.5*xn.^2+log(0.5*erfcx(-xn/sqrt(2)));
end
