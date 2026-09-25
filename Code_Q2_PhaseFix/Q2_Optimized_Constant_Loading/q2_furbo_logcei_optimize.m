function out = q2_furbo_logcei_optimize(oracle,d,budget,nInit,seed, ...
    initialPoint,checkpointFile,resumeFile,initialDesignBounds)
% FuRBO检查点信赖域 + Ament解析LogCEI的混合创新算法。
% FuRBO部分：检查点排序、前10%最小包围超矩形、R倍增/减半。
% LogCEI部分：稳定解析LogEI + 稳定log(P(feasible))。
% 目标时间只用真实成功样本拟合；最终最优只从真实可行评价中选择。

rng(seed,'twister');
if nargin<6 || isempty(initialPoint), initialPoint=0.5*ones(1,d); end
if nargin<7, checkpointFile=''; end
if nargin<8, resumeFile=''; end
assert(nInit>=1 && nInit<=budget && fix(nInit)==nInit && fix(budget)==budget);
assert(isrow(initialPoint) && numel(initialPoint)==d && ...
    all(isfinite(initialPoint)) && all(initialPoint>=0 & initialPoint<=1));
assert(exist('lhsdesign','file')==2 && exist('fitrgp','file')==2, ...
    'Statistics and Machine Learning Toolbox is required.');

X=nan(budget,d); Y=nan(budget,1); G=nan(budget,1);
success=false(budget,1); records=cell(budget,1); n=0;
radiusHistory=nan(budget,1); trustRegionLower=nan(budget,d);
trustRegionUpper=nan(budget,d);

X0=lhsdesign(nInit,d,'criterion','maximin','iterations',10);
% Optional focused Latin-hypercube start around a known feasible reference.
% The adaptive trust region and candidate domain remain the full [0,1]^d.
if nargin>=9 && ~isempty(initialDesignBounds)
    assert(isequal(size(initialDesignBounds),[2 d]) && ...
        all(isfinite(initialDesignBounds(:))) && ...
        all(initialDesignBounds(1,:)>=0) && ...
        all(initialDesignBounds(2,:)<=1) && ...
        all(initialDesignBounds(2,:)>initialDesignBounds(1,:)));
    X0=initialDesignBounds(1,:)+X0.* ...
        (initialDesignBounds(2,:)-initialDesignBounds(1,:));
end
X0(1,:)=initialPoint;
sourceFile=resumeFile;
if ~isempty(checkpointFile) && isfile(checkpointFile)
    sourceFile=checkpointFile;
end
if ~isempty(sourceFile)
    saved=load(sourceFile);
    if isfield(saved,'report')
        cp=saved.report.search;
    else
        assert(isfield(saved,'checkpoint'),'Resume file has no search data.');
        cp=saved.checkpoint;
    end
    n=cp.nfe;
    assert(n>=1 && n<=budget && size(cp.X,1)==n && ...
        size(cp.X,2)==d && numel(cp.evaluations)==n, ...
        'Resume file dimensions or budget do not match.');
    assert(max(abs(cp.X(1:min(n,nInit),:)- ...
        X0(1:min(n,nInit),:)),[],'all')<1e-12, ...
        'Resume file does not match the seeded initial design.');
    X(1:n,:)=cp.X; Y(1:n)=cp.Y; G(1:n)=cp.G;
    success(1:n)=cp.success; records(1:n)=cp.evaluations;
    if isfield(cp,'radius_history')
        radiusHistory(1:n)=cp.radius_history(1:n);
        trustRegionLower(1:n,:)=cp.trust_region_lower(1:n,:);
        trustRegionUpper(1:n,:)=cp.trust_region_upper(1:n,:);
    end
    % Old result files do not contain the RNG state.  New checkpoints do.
    if isfield(cp,'rng_state')
        rng(cp.rng_state);
    else
        rng(seed+n,'twister');
    end
    fprintf('Continuing from %d observed evaluations in %s\n',n,sourceFile);
end
for i=n+1:nInit, evaluate(X0(i,:)); end

% FuRBO论文实验设置：R0=1, tau_s=2, tau_f=3, top P=10%。
radius=1;
successThreshold=2;
failureThreshold=3;
successCounter=0;
failureCounter=0;
inspectorFraction=0.10;
nInspectors=1000;
nCandidates=512;
minimumRadius=5e-8;
if n>nInit
    [radius,successCounter,failureCounter]=replay_trust_updates( ...
        Y(1:n),G(1:n),success(1:n),nInit, ...
        successThreshold,failureThreshold,minimumRadius);
    fprintf('Recovered FuRBO radius %.9g; success/failure counters %d/%d\n', ...
        radius,successCounter,failureCounter);
end

while n<budget
    % Numerical noncompletion is not evidence about the physical constraint.
    usable=true(n,1);
    for k=1:n
        if isfield(records{k},'solver_failed') && records{k}.solver_failed
            usable(k)=false;
        end
    end
    gModel=fit_gp(X(find(usable),:),G(usable)); %#ok<FNDSB>
    good=find(success(1:n) & isfinite(Y(1:n)));
    yModel=fit_gp(X(good,:),Y(good));
    if isempty(good)
        [previousMetric,centerIndex]=min(G(1:n));
        previousBest=inf;
    else
        [previousBest,j]=min(Y(good));
        centerIndex=good(j);
        previousMetric=previousBest;
    end
    center=X(centerIndex,:);

    inspectors=uniform_ball_samples(center,radius,nInspectors,d);
    [gMean,~]=gp_predict(gModel,inspectors);
    if isempty(yModel)
        yMean=zeros(nInspectors,1);
    else
        [yMean,~]=gp_predict(yModel,inspectors);
    end
    order=rank_predictions(yMean,gMean,~isempty(good) && ~isempty(yModel));
    nTop=max(2,ceil(inspectorFraction*nInspectors));
    bestInspectors=inspectors(order(1:nTop),:);
    lower=min(bestInspectors,[],1);
    upper=max(bestInspectors,[],1);
    [lower,upper]=enforce_minimum_width(lower,upper,center,0.01);
    trustRegionLower(n+1,:)=lower;
    trustRegionUpper(n+1,:)=upper;
    radiusHistory(n+1)=radius;

    unitPool=lhsdesign(nCandidates,d,'criterion','maximin','iterations',2);
    pool=lower+(upper-lower).*unitPool;
    pool(end,:)=center;
    [gMean,gSd]=gp_predict(gModel,pool);
    logP=log_normal_feasibility(gMean,gSd);
    if isempty(good) || isempty(yModel)
        % 尚无足够成功时间样本时，先按可行概率并辅以约束不确定度。
        score=logP+0.05*min(gSd,1);
    else
        [yMean,ySd]=gp_predict(yModel,pool);
        score=analytic_log_expected_improvement(previousBest,yMean,ySd)+logP;
    end
    [~,candidateOrder]=sort(score,'descend');
    chosen=candidateOrder(1);
    for j=1:numel(candidateOrder)
        candidate=candidateOrder(j);
        if all(vecnorm(X(1:n,:)-pool(candidate,:),2,2)>1e-6)
            chosen=candidate;
            break;
        end
    end

    newIndex=evaluate(pool(chosen,:));
    if isempty(good)
        improved=G(newIndex)<previousMetric-1e-8;
    else
        improved=success(newIndex) && Y(newIndex)<previousBest-1e-8;
    end
    if improved
        successCounter=successCounter+1;
        failureCounter=0;
    else
        failureCounter=failureCounter+1;
        successCounter=0;
    end
    if successCounter>=successThreshold
        radius=2*radius;
        successCounter=0;
        failureCounter=0;
    elseif failureCounter>=failureThreshold
        radius=0.5*radius;
        successCounter=0;
        failureCounter=0;
    end
    if radius<=minimumRadius
        radius=1; % 论文允许停止或重启；有限预算下选择重启。
        successCounter=0;
        failureCounter=0;
    end
end

good=find(success & isfinite(Y));
best=NaN; bestIndex=[];
if ~isempty(good)
    [best,j]=min(Y(good));
    bestIndex=good(j);
end
out.nfe=n;
out.best_time=best;
out.best_x=X(bestIndex,:);
out.success_count=numel(good);
out.best_is_observed=~isempty(bestIndex) && success(bestIndex) && isfinite(best);
out.X=X; out.Y=Y; out.G=G; out.success=success;
out.evaluations=records;
out.radius_history=radiusHistory;
out.trust_region_lower=trustRegionLower;
out.trust_region_upper=trustRegionUpper;
out.rng_state=rng;
out.method=['FuRBO inspector-defined trust region with analytic LogCEI; ', ...
    'hybrid method, not FuRBO Thompson sampling reproduction.'];

    function idx=evaluate(x)
        assert(all(isfinite(x)) && all(x>=0) && all(x<=1) && n<budget);
        fprintf('\n========== 优化评估 %d / %d 开始 ==========\n',n+1,budget);
        fprintf('归一化参数：'); fprintf(' %.6f',x); fprintf('\n');
        evaluationClock=tic;
        r=oracle(x);
        wallSeconds=toc(evaluationClock);
        n=n+1; idx=n; X(n,:)=x;
        success(n)=logical(r.success) && isfinite(r.startup_time);
        if success(n), Y(n)=r.startup_time; else, assert(isnan(r.startup_time)); end
        G(n)=r.constraint_violation;
        records{n}=r;
        fprintf(['========== 优化评估 %d / %d 完成：success=%d，', ...
            '启动时间=%g s，实际耗时=%.2f s，约束指标=%g ==========\n'], ...
            n,budget,success(n),Y(n),wallSeconds,G(n));
        assert(isfinite(G(n)));
        if ~isempty(checkpointFile)
            checkpoint=struct('nfe',n,'budget',budget,'X',X(1:n,:), ...
                'Y',Y(1:n),'G',G(1:n),'success',success(1:n), ...
                'evaluations',{records(1:n)},'method','FuRBO-LogCEI hybrid', ...
                'radius_history',radiusHistory(1:n), ...
                'trust_region_lower',trustRegionLower(1:n,:), ...
                'trust_region_upper',trustRegionUpper(1:n,:), ...
                'rng_state',rng);
            save(checkpointFile,'checkpoint');
        end
    end
end

function [radius,successCounter,failureCounter]=replay_trust_updates( ...
    Y,G,success,nInit,successThreshold,failureThreshold,minimumRadius)
radius=1; successCounter=0; failureCounter=0;
for k=nInit+1:numel(Y)
    oldGood=find(success(1:k-1) & isfinite(Y(1:k-1)));
    if isempty(oldGood)
        improved=G(k)<min(G(1:k-1))-1e-8;
    else
        improved=success(k) && Y(k)<min(Y(oldGood))-1e-8;
    end
    if improved
        successCounter=successCounter+1; failureCounter=0;
    else
        failureCounter=failureCounter+1; successCounter=0;
    end
    if successCounter>=successThreshold
        radius=2*radius; successCounter=0; failureCounter=0;
    elseif failureCounter>=failureThreshold
        radius=0.5*radius; successCounter=0; failureCounter=0;
    end
    if radius<=minimumRadius
        radius=1; successCounter=0; failureCounter=0;
    end
end
end

function samples=uniform_ball_samples(center,radius,nSamples,d)
directions=randn(nSamples,d);
norms=vecnorm(directions,2,2);
norms(norms==0)=1;
directions=directions./norms;
% 按FuRBO正文：方向单位化后乘以[0,R]均匀随机尺度。
distances=radius*rand(nSamples,1);
samples=center+directions.*distances;
samples=min(1,max(0,samples));
end

function order=rank_predictions(yMean,gMean,useObjective)
feasible=find(gMean<=0);
infeasible=find(gMean>0);
if useObjective && ~isempty(feasible)
    [~,j]=sort(yMean(feasible),'ascend');
    feasible=feasible(j);
else
    [~,j]=sort(gMean(feasible),'ascend');
    feasible=feasible(j);
end
if ~isempty(infeasible)
    scale=max(abs(gMean(infeasible)));
    if scale<=0, scale=1; end
    violation=gMean(infeasible)/scale;
    [~,j]=sort(violation,'ascend');
    infeasible=infeasible(j);
end
order=[feasible;infeasible];
end

function [lower,upper]=enforce_minimum_width(lower,upper,center,minWidth)
tooNarrow=(upper-lower)<minWidth;
lower(tooNarrow)=center(tooNarrow)-minWidth/2;
upper(tooNarrow)=center(tooNarrow)+minWidth/2;
lower=max(0,lower); upper=min(1,upper);
end

function model=fit_gp(X,y)
if size(X,1)<5 || numel(unique(y))<2
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
    mu=double(mu(:)); sd=max(0,double(sd(:)));
end
end

function logP=log_normal_feasibility(mu,sd)
logP=-inf(size(mu));
zero=sd<=1e-12;
logP(zero & mu<=0)=0;
positive=~zero;
z=-mu(positive)./sd(positive);
values=zeros(size(z));
regular=z>=-1;
values(regular)=log(normcdf(z(regular)));
tail=~regular;
values(tail)=log(0.5)+log(erfcx(-z(tail)/sqrt(2)))-0.5*z(tail).^2;
logP(positive)=min(0,values);
end

function logEI=analytic_log_expected_improvement(best,mu,sd)
delta=best-mu;
logEI=-inf(size(mu));
zero=sd<=1e-12;
positiveImprovement=zero & delta>0;
logEI(positiveImprovement)=log(delta(positiveImprovement));
positive=~zero;
z=delta(positive)./sd(positive);
logEI(positive)=log_h(z)+log(sd(positive));
end

function value=log_h(z)
value=zeros(size(z));
c1=0.5*log(2*pi);
c2=0.5*log(pi/2);
threshold=-1/sqrt(eps);
regular=z>-1;
value(regular)=log(normpdf(z(regular))+z(regular).*normcdf(z(regular)));
middle=z<=-1 & z>threshold;
if any(middle)
    zm=z(middle);
    a=log(erfcx(-zm/sqrt(2)).*abs(zm))+c2;
    value(middle)=-0.5*zm.^2-c1+log1mexp(a);
end
tail=z<=threshold;
if any(tail)
    zt=z(tail);
    value(tail)=-0.5*zt.^2-c1-2*log(abs(zt));
end
end

function y=log1mexp(x)
% 稳定计算log(1-exp(x))，定义域x<=0。
x=min(x,0);
y=zeros(size(x));
near=x>-log(2);
y(near)=log(-expm1(x(near)));
y(~near)=log1p(-exp(x(~near)));
end
