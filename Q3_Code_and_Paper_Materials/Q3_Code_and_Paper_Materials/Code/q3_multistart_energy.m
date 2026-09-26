function result=q3_multistart_energy(generations)
% Two genuinely different DE populations, initialized from saved physics.
    if nargin<1, generations=2; end
    assert(generations>=1 && generations==floor(generations));
    codeDir=fileparts(mfilename('fullpath')); out=fullfile(codeDir,'results');
    cfg=q3_optimization_config('smoke'); opt=cfg.simOpt;
    opt.RelTol=1e-5; opt.AbsTol=1e-9;
    physics={'pemfc_stack5_simulate_q3.m','thermal_temperature_state_ice.m', ...
        'water_ice_state_ice.m','gas_transport_state_ice.m','pemfc_setup_ice.m', ...
        'q3_auxiliary_heating.m','q3_evaluate_energy6.m','q2_current_strategy.m'};
    dates=zeros(size(physics));
    for k=1:numel(physics), d=dir(fullfile(codeDir,physics{k})); dates(k)=d.datenum; end
    manifestPath=fullfile(out,'q3_multistart_manifest.mat');
    if exist(manifestPath,'file')
        a=load(manifestPath,'manifest'); manifest=a.manifest;
        assert(isequal(manifest.physicsDates,dates) && isequaln(manifest.options,opt));
    else
        a=load(fullfile(out,'q3_final_energy_table.mat'),'verified'); pool={a.verified{2}};
        keys={mat2str([pool{1}.q,pool{1}.th],12)};
        files=dir(fullfile(out,'energy6_coheat_*.mat'));
        for k=1:numel(files)
            a=load(fullfile(out,files(k).name),'c','signature'); c=a.c;
            if ~isequal(a.signature.modified,dates) || ~valid_record(c,opt), continue; end
            key=mat2str([c.q,c.th],12);
            if any(strcmp(keys,key)), continue; end
            pool{end+1}=c; keys{end+1}=key; %#ok<AGROW>
        end
        good=find(cellfun(@(c) strcmp(c.status,'success'),pool));
        bad=find(cellfun(@(c) strcmp(c.status,'physical_infeasible'),pool));
        assert(numel(good)>=8 && numel(bad)>=4,'Insufficient physical seed records.');
        [~,j]=min(cellfun(@(c) c.energyJ,pool(good))); first=good(j);
        features=zeros(numel(pool),6);
        for k=1:numel(pool)
            c=pool{k}; duration=c.th;
            if strcmp(c.status,'success'), duration=min(duration,c.startupTimeS); end
            features(k,:)=[c.q,2*duration/(60+(20-9)/.3)];
        end
        selected1=first;
        for k=2:8, selected1(end+1)=diverse(good,selected1,features); end
        selected2=first;
        for k=2:4, selected2(end+1)=diverse(good,selected2,features); end
        [~,order]=sort(cellfun(@(c) c.violation,pool(bad)));
        nearBad=bad(order(1:min(16,numel(bad))));
        for k=5:8, selected2(end+1)=diverse(nearBad,selected2,features); end
        manifest=struct('options',opt,'physicsDates',dates, ...
            'records',{{pool(selected1),pool(selected2)}},'seeds',[7311 8429], ...
            'poolCount',numel(pool),'feasibleCount',numel(good));
        save(manifestPath,'manifest','-v7.3');
    end
    fprintf('Multi-start: %d saved records, %d feasible; two populations of 8.\n', ...
        manifest.poolCount,manifest.feasibleCount);
    result=struct('runs',{cell(2,1)},'best',[],'newTrialCount',0);
    for start=1:2
        records=manifest.records{start}; P=zeros(8,6);
        for k=1:8, P(k,:)=[records{k}.q,records{k}.th/(60+(20-9)/.3)]; end
        archive=fullfile(out,sprintf('q3_multistart_%02d.mat',start));
        existing=0;
        if exist(archive,'file'), a=load(archive,'run'); existing=a.run.evaluations; end
        settings=struct('populationSize',8,'seed',manifest.seeds(start), ...
            'initialPopulation',P,'initialRecords',{records}, ...
            'simOpt',opt,'parallel',true,'batchSize',2,'archivePath',archive);
        fprintf('Starting population %d, target %d DE generations.\n',start,generations);
        r=q3_joint_de_energy(max(0,8*(1+generations)-existing),settings);
        result.runs{start}=r;
        result.newTrialCount=result.newTrialCount+max(0,r.evaluations-8);
        if isempty(result.best) || r.best.energyJ<result.best.energyJ, result.best=r.best; end
        save(fullfile(out,'q3_multistart_result.mat'),'result','-v7.3');
    end
    c=result.best; save(fullfile(out,'q3_multistart_best.mat'),'c','opt','-v7.3');
    rows=table;
    for start=1:2
        h=result.runs{start}.history; Z=vertcat(h.z);
        t=table(repmat(start,numel(h),1),(1:numel(h)).', ...
            [h.generation].',Z(:,1),Z(:,2),Z(:,3),Z(:,4),Z(:,5), ...
            [h.th].',string({h.status}).',[h.energyJ].',[h.startupTimeS].', ...
            [h.generation].'==0,'VariableNames',{'start','evaluation','generation', ...
            'q1','q2','q3','q4','q5','th','status','energyJ','startupTimeS','savedSeed'});
        rows=[rows;t]; %#ok<AGROW>
    end
    writetable(rows,fullfile(out,'q3_multistart_history.csv'),'Encoding','UTF-8');
end
function yes=valid_record(c,opt)
    yes=~isempty(c.result) && ismember(c.status,{'success','physical_infeasible'}) && ...
        c.th>=0 && c.th<=60+(20-9)/.3 && ~c.result.solverTerminatedUnexpectedly;
    if ~yes, return; end
    fields={'initialTemperatureC','ambientTemperatureC','gammaIce','kFreeze', ...
        'nCellLayer','iceScope','temperatureEventMarginK'};
    for k=1:numel(fields)
        f=fields{k}; yes=yes && isequaln(c.result.options.(f),opt.(f));
    end
    yes=yes && c.result.options.RelTol<=1e-5 && c.result.options.AbsTol<=1e-9;
end
function index=diverse(available,selected,F)
    available=setdiff(available,selected,'stable');
    assert(~isempty(available)); distances=inf(size(available));
    for k=selected
        distances=min(distances,sum((F(available,:)-F(k,:)).^2,2).');
    end
    [~,j]=max(distances); index=available(j);
end
