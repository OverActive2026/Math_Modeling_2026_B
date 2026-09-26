function run=q3_energy_budget_search(maxNewEvaluations)
% Search between earlier power levels at a lower common energy budget.
    if nargin<1, maxNewEvaluations=6; end
    out=fullfile(fileparts(mfilename('fullpath')),'results');
    path=fullfile(out,'q3_budget_checkpoint.mat');
    cfg=q3_optimization_config('smoke'); opt=cfg.simOpt;
    opt.RelTol=1e-6; opt.AbsTol=1e-10; opt.MaxStep=.0125;
    if exist(path,'file'), a=load(path,'run'); run=a.run;
    else
        a=load(fullfile(out,'q3_final_energy_table.mat'),'verified');
        run=struct('evaluations',0,'best',a.verified{2},'options',opt,'history',struct([]));
    end
    assert(isequaln(run.options,opt));
    designs=[.7 .7 .7 .7 .7;.6 .6 .6 .6 .6;.9 .9 .9 .9 .9;.85 1 1 1 .85];
    for k=1:maxNewEvaluations
        i=run.evaluations+1;
        if i<=4
            q=designs(i,:); budget=run.best.energyJ-5;
            th=budget/(25*sum(q)); label='lower_energy_budget';
        elseif i==7 || i==8
            levels=[.5 .4]; q=levels(i-6)*ones(1,5);
            budget=run.best.energyJ-5;
            th=budget/(25*sum(q)); label='long_duration_budget_gap';
        else
            q=run.best.q; q(5)=min(1,max(0,q(5)+ ...
                (-1)^i*.00025/(2^floor((i-5)/2))));
            th=40; budget=NaN; label='fine_end5';
        end
        fprintf('Budget search %d %s q=%s th=%.8f\n',i,label,mat2str(q,9),th);
        c=q3_cached_energy6('coheat',q,th,opt);
        accepted=strcmp(c.status,'success') && (c.energyJ<run.best.energyJ || ...
            (c.energyJ==run.best.energyJ && c.startupTimeS<run.best.startupTimeS));
        if accepted, run.best=c; end
        save(fullfile(out,sprintf('q3_budget_case_%03d.mat',i)), ...
            'c','opt','budget','label','-v7.3');
        h=struct('index',i,'q',q,'th',th,'status',c.status, ...
            'energyJ',c.energyJ,'startupTimeS',c.startupTimeS, ...
            'accepted',accepted,'bestEnergyJ',run.best.energyJ);
        if isempty(run.history), run.history=h; else, run.history(end+1)=h; end
        run.evaluations=i;
        save([path '.tmp'],'run','-v7.3'); movefile([path '.tmp'],path,'f');
        writetable(struct2table(run.history),fullfile(out,'q3_budget_history.csv'),'Encoding','UTF-8');
        fprintf('Budget search %d: %s E=%.6f t=%.6f\n',i,c.status,c.energyJ,c.startupTimeS);
    end
end
