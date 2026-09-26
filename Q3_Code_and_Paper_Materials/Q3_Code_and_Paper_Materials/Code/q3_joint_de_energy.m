function run = q3_joint_de_energy(maxNewEvaluations,settings)
% Six-variable DE/rand/1/bin, with exact resume and optional parallel batches.
% The physical evaluator observes startup independently of heater duration.
    if nargin<1, maxNewEvaluations=24; end
    if nargin<2, settings=struct(); end
    cfg=q3_optimization_config('smoke');
    cfg.simOpt.RelTol=1e-5; cfg.simOpt.AbsTol=1e-9;
    defaults=struct('populationSize',12,'seed',20260926, ...
        'energyTieToleranceJ',0.1,'batchSize',1,'parallel',false, ...
        'evaluator',@q3_cached_energy6,'simOpt',cfg.simOpt, ...
        'archivePath',fullfile(fileparts(mfilename('fullpath')), ...
        'results','q3_joint_de_checkpoint.mat'));
    fields=fieldnames(defaults);
    for k=1:numel(fields)
        if ~isfield(settings,fields{k}), settings.(fields{k})=defaults.(fields{k}); end
    end
    if ~isfield(settings,'eliteCandidate'), settings.eliteCandidate=[]; end
    if ~isfield(settings,'initialPopulation'), settings.initialPopulation=[]; end
    if ~isfield(settings,'initialRecords'), settings.initialRecords={}; end
    assert(maxNewEvaluations>=0 && maxNewEvaluations==floor(maxNewEvaluations));
    np=settings.populationSize;
    assert(np>=4 && np==floor(np));
    cap=60+(20-9)/0.3;
    fingerprint=struct('version',1,'populationSize',np,'seed',settings.seed, ...
        'energyTieToleranceJ',settings.energyTieToleranceJ, ...
        'evaluator',func2str(settings.evaluator),'simOpt',settings.simOpt);
    if ~isempty(settings.initialPopulation)
        P0=settings.initialPopulation;
        assert(isequal(size(P0),[np 6]) && all(isfinite(P0),'all') && ...
            all(P0>=0 & P0<=1,'all'),'Invalid custom population.');
        fingerprint.initialPopulation=P0;
        if ~isempty(settings.initialRecords)
            assert(numel(settings.initialRecords)==np);
            seedKeys=cell(np,1);
            for k=1:np
                c=settings.initialRecords{k};
                assert(isequal(c.q,P0(k,1:5)) && abs(c.th-P0(k,6)*cap)<1e-8);
                assert(ismember(c.status,{'success','physical_infeasible'}));
                seedKeys{k}={c.status,c.energyJ,c.startupTimeS,c.violation};
            end
            fingerprint.initialRecordKeys=seedKeys;
        end
    else
        assert(isempty(settings.initialRecords),'Seed records require an explicit population.');
    end
    if exist(settings.archivePath,'file')
        stored=load(settings.archivePath,'run'); run=stored.run;
        if ~isequaln(run.fingerprint,fingerprint)
            error('q3_joint_de_energy:ArchiveMismatch','Checkpoint settings differ.');
        end
    else
        rng(settings.seed,'twister');
        P=rand(np,6);
        seeds=[1 1 1 1 1 24;1 1 1 1 1 30.2537663; ...
            1 .8 .5 .8 1 28;1 .5 .2 .5 1 34;.8 .8 .8 .8 .8 28; ...
            1 .85 .7 .85 1 3770/(25*4.4); ...
            1 .6 .2 .6 1 3770/(25*3.4); ...
            1 1 0 1 1 3770/(25*4); ...
            .8 .8 .8 .8 .8 3770/(25*4); ...
            1 .7 .3 .7 1 45; .6 .8 1 .8 .6 45; ...
            1 .9 .8 .9 1 40];
        n=min(size(seeds,1),np); P(1:n,:)=[seeds(1:n,1:5),seeds(1:n,6)/cap];
        if ~isempty(settings.initialPopulation), P=settings.initialPopulation; end
        run=struct('fingerprint',fingerprint,'population',P, ...
            'records',{cell(np,1)},'parentPopulation',[], ...
            'parentRecords',{cell(0,1)},'pendingIndex',1,'generation',0, ...
            'evaluations',0,'rngState',rng,'best',[], ...
            'history',struct([]),'pendingTrials',struct([]));
        if ~isempty(settings.initialRecords), run.seedRecords=settings.initialRecords; end
    end
    % Import an independently verified improvement only between generations.
    elite=settings.eliteCandidate;
    if ~isempty(elite) && ~isempty(elite.result)
        fixed={'initialTemperatureC','ambientTemperatureC','gammaIce','kFreeze', ...
            'nCellLayer','iceScope','temperatureEventMarginK'};
        for f=1:numel(fixed)
            field=fixed{f};
            assert(isequaln(elite.result.options.(field),settings.simOpt.(field)), ...
                'Elite candidate uses different physical assumptions: %s',field);
        end
    end
    if ~isempty(elite) && strcmp(elite.status,'success') && ...
            isempty(run.pendingTrials) && run.pendingIndex>np && ...
            better(elite,run.best,0)
        slot=find(~cellfun(@(c) strcmp(c.status,'success'),run.records),1);
        if isempty(slot)
            [~,slot]=max(cellfun(@(c) c.energyJ,run.records));
        end
        assert(all(elite.q>=0 & elite.q<=1) && elite.th>=0 && elite.th<=cap);
        run.population(slot,:)=[elite.q,elite.th/cap];
        small=elite; small.result=[]; run.records{slot}=small;
        run.best=elite;
        if ~isfield(run,'eliteImports'), run.eliteImports={}; end
        run.eliteImports{end+1}=struct('slot',slot,'afterEvaluations',run.evaluations, ...
            'q',elite.q,'th',elite.th,'energyJ',elite.energyJ);
    end
    targetCount=run.evaluations+maxNewEvaluations;
    evaluator=settings.evaluator; opt=settings.simOpt;
    if settings.parallel && isempty(gcp('nocreate'))
        parpool('Processes',settings.batchSize);
    end
    while run.evaluations<targetCount
        if isempty(run.pendingTrials)
            if run.pendingIndex>np
                run.generation=run.generation+1;
                run.parentPopulation=run.population;
                run.parentRecords=run.records;
                run.pendingIndex=1;
            end
            n=min([settings.batchSize,np-run.pendingIndex+1, ...
                targetCount-run.evaluations]);
            rng(run.rngState);
            for b=1:n
                i=run.pendingIndex+b-1;
                parents=[]; donors=[]; target=[]; mask=false(1,6);
                if run.generation==0
                    z=run.population(i,:);
                else
                    others=setdiff(1:np,i);
                    parents=others(randperm(numel(others),3));
                    donors=run.parentPopulation(parents,:);
                    mutant=min(1,max(0,donors(1,:)+0.65*(donors(2,:)-donors(3,:))));
                    mask=rand(1,6)<0.85; mask(randi(6))=true;
                    target=run.parentPopulation(i,:);
                    z=target; z(mask)=mutant(mask);
                end
                trial=struct('z',z,'populationIndex',i, ...
                    'generation',run.generation,'parentIndices',parents, ...
                    'donorZ',donors,'targetZ',target,'crossoverMask',mask, ...
                    'rngAfter',rng);
                if isempty(run.pendingTrials), run.pendingTrials=trial;
                else, run.pendingTrials(end+1)=trial;
                end
            end
            checkpoint();
        end
        n=min(numel(run.pendingTrials),targetCount-run.evaluations);
        trials=run.pendingTrials(1:n);
        results=cell(n,1);
        seedRecords={};
        if isfield(run,'seedRecords'), seedRecords=run.seedRecords; end
        if settings.parallel
            parfor b=1:n
                if trials(b).generation==0 && ~isempty(seedRecords)
                    results{b}=seedRecords{trials(b).populationIndex};
                else
                    results{b}=evaluate_one(evaluator,trials(b).z,cap,opt);
                end
            end
        else
            for b=1:n
                if trials(b).generation==0 && ~isempty(seedRecords)
                    results{b}=seedRecords{trials(b).populationIndex};
                else
                    results{b}=evaluate_one(evaluator,trials(b).z,cap,opt);
                end
            end
        end
        for b=1:n
            h=trials(b); candidate=results{b}; i=h.populationIndex;
            if h.generation==0
                accepted=true;
            else
                accepted=better(candidate,run.parentRecords{i}, ...
                    settings.energyTieToleranceJ);
            end
            if isempty(run.best) || better(candidate,run.best,0)
                run.best=candidate;
            end
            small=candidate;
            if isfield(small,'result'), small.result=[]; end
            if accepted
                run.population(i,:)=h.z;
                run.records{i}=small;
            end
            h=rmfield(h,'rngAfter');
            h.th=candidate.th; h.status=candidate.status;
            h.energyJ=candidate.energyJ; h.startupTimeS=candidate.startupTimeS;
            h.violation=candidate.violation; h.accepted=accepted;
            if isempty(run.history), run.history=h;
            else, run.history(end+1)=h;
            end
            run.evaluations=run.evaluations+1;
            run.pendingIndex=i+1;
            run.rngState=trials(b).rngAfter;
            run.pendingTrials(1)=[];
            checkpoint();
            fprintf('Joint DE eval %d gen %d q=%s th=%.4f: %s E=%.3f t=%.3f\n', ...
                run.evaluations,h.generation,mat2str(candidate.q,4),candidate.th, ...
                candidate.status,candidate.energyJ,candidate.startupTimeS);
        end
    end
    checkpoint();

    function checkpoint()
        save([settings.archivePath '.tmp'],'run','-v7.3');
        movefile([settings.archivePath '.tmp'],settings.archivePath,'f');
    end
end

function r=evaluate_one(evaluator,z,cap,opt)
    q=z(1:5); th=z(6)*cap;
    try
        r=evaluator('coheat',q,th,opt);
    catch ME
        r=struct('mode','coheat','q',q,'th',th,'status','numerical_unknown', ...
            'energyJ',Inf,'startupTimeS',NaN,'violation',Inf, ...
            'reason',[ME.identifier ': ' ME.message],'result',[]);
    end
end

function yes=better(a,b,tolerance)
    if isempty(b), yes=true; return; end
    rank=@(r) 2*strcmp(r.status,'success')+strcmp(r.status,'physical_infeasible');
    ra=rank(a); rb=rank(b);
    if ra~=rb, yes=ra>rb; return; end
    if ra==2
        if abs(a.energyJ-b.energyJ)>tolerance
            yes=a.energyJ<b.energyJ;
        else
            yes=a.startupTimeS<b.startupTimeS;
        end
    elseif ra==1
        yes=a.violation<b.violation;
    else
        yes=false;
    end
end
