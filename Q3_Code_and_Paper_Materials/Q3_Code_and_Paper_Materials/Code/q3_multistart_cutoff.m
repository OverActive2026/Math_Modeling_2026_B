function run=q3_multistart_cutoff()
% Probe shared heater shutoff at lower budgets for distinct feasible profiles.
    out=fullfile(fileparts(mfilename('fullpath')),'results');
    source=load(fullfile(out,'q3_multistart_result.mat'),'result');
    m=load(fullfile(out,'q3_multistart_manifest.mat'),'manifest');
    opt=m.manifest.options;
    opt.RelTol=1e-6; opt.AbsTol=1e-10; opt.MaxStep=.0125;
    archive=fullfile(out,'q3_multistart_cutoff.mat');
    if exist(archive,'file')
        a=load(archive,'run'); run=a.run;
        assert(isequaln(run.options,opt));
    else
        final=load(fullfile(out,'q3_final_energy_table.mat'),'verified');
        pool={final.verified{2},source.result.best};
        for s=1:numel(source.result.runs)
            r=source.result.runs{s};
            pool=[pool,r.seedRecords(:).']; %#ok<AGROW>
            for k=1:numel(r.history)
                h=r.history(k);
                if h.generation>0 && strcmp(h.status,'success')
                    pool{end+1}=q3_cached_energy6('coheat',h.z(1:5), ...
                        h.th,m.manifest.options); %#ok<AGROW>
                end
            end
        end
        pool=pool(cellfun(@(c) strcmp(c.status,'success'),pool));
        [~,order]=sort(cellfun(@(c) c.energyJ,pool)); pool=pool(order);
        selected=1;
        for k=2:numel(pool)
            distance=cellfun(@(c) norm(c.q-pool{k}.q),pool(selected));
            if min(distance)>.05, selected(end+1)=k; end %#ok<AGROW>
            if numel(selected)==3, break; end
        end
        Q=cell2mat(cellfun(@(c) c.q,pool(selected).','UniformOutput',false));
        baseline=min(final.verified{2}.energyJ,source.result.best.energyJ);
        run=struct('options',opt,'Q',Q,'budgetJ',baseline-[.5 10], ...
            'evaluations',0,'history',struct([]),'best',final.verified{2});
        save(archive,'run','-v7.3');
    end
    total=2*size(run.Q,1);
    if isempty(gcp('nocreate')), parpool('local',2); end
    while run.evaluations<total
        indices=run.evaluations+1:min(run.evaluations+2,total);
        candidates=cell(size(indices));
        parfor j=1:numel(indices)
            i=indices(j); q=run.Q(ceil(i/2),:);
            budget=run.budgetJ(1+mod(i-1,2));
            th=budget/(25*sum(q));
            candidates{j}=q3_cached_energy6('coheat',q,th,opt);
        end
        for j=1:numel(indices)
            i=indices(j); c=candidates{j};
            budget=run.budgetJ(1+mod(i-1,2));
            save(fullfile(out,sprintf('q3_multistart_stop_%02d.mat',i)), ...
                'c','opt','budget','-v7.3');
            accepted=strcmp(c.status,'success') && c.energyJ<run.best.energyJ;
            if accepted, run.best=c; end
            h=struct('index',i,'q',c.q,'th',c.th,'targetEnergyJ',budget, ...
                'status',c.status,'energyJ',c.energyJ, ...
                'startupTimeS',c.startupTimeS,'reason',c.reason,'accepted',accepted);
            if isempty(run.history), run.history=h; else, run.history(end+1)=h; end
            run.evaluations=i;
            save(archive,'run','-v7.3');
            fprintf('Cutoff %d/%d q=%s th=%.9f %s E=%.6f t=%.6f reason=%s\n', ...
                i,total,mat2str(c.q,7),c.th,c.status,c.energyJ,c.startupTimeS,c.reason);
        end
        writetable(struct2table(run.history), ...
            fullfile(out,'q3_multistart_cutoff_history.csv'),'Encoding','UTF-8');
    end
end
