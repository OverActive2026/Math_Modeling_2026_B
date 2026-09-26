function run=q3_continue_energy_search(maxNewEvaluations)
% Resume targeted full-model searches around the incumbent and omitted faces.
% All heaters share one duration; no extra physical constraints are imposed.
    if nargin<1, maxNewEvaluations=11; end
    assert(maxNewEvaluations>=0 && maxNewEvaluations==floor(maxNewEvaluations));
    out=fullfile(fileparts(mfilename('fullpath')),'results');
    path=fullfile(out,'q3_targeted_checkpoint.mat');
    cfg=q3_optimization_config('smoke'); opt=cfg.simOpt;
    opt.RelTol=1e-6; opt.AbsTol=1e-10; opt.MaxStep=.0125;
    if exist(path,'file')
        a=load(path,'run'); run=a.run;
        assert(isequaln(run.options,opt),'Search options changed.');
        latest=load(fullfile(out,'q3_final_energy_table.mat'),'verified');
        if latest.verified{1}.energyJ<run.bestPre.energyJ
            run.bestPre=latest.verified{1};
        end
        if latest.verified{2}.energyJ<run.bestCo.energyJ
            run.bestCo=latest.verified{2};
        end
    else
        a=load(fullfile(out,'q3_final_energy_table.mat'),'verified');
        run=struct('options',opt,'bestPre',a.verified{1}, ...
            'bestCo',a.verified{2},'anchor',a.verified{2}, ...
            'evaluations',0,'history',struct([]));
    end
    if ~isfield(run,'proposals'), run.proposals=run.evaluations; end
    if ~isfield(run,'skippedNoChange'), run.skippedNoChange=0; end
    target=run.evaluations+maxNewEvaluations;
    while run.evaluations<target
        i=run.proposals+1; run.proposals=i;
        mode='coheat'; q=run.bestCo.q; th=40;
        if i==1
            a=load(fullfile(out,'q3_end_margin_01.mat'),'c');
            q=a.c.q; label='retry_unresolved_end5';
        elseif i==2
            a=load(fullfile(out,'q3_end_margin_01.mat'),'c');
            q=(run.anchor.q+a.c.q)/2; label='joint_end_margin';
        elseif i==3
            th=run.bestCo.startupTimeS-.01; label='incumbent_early_shutoff';
        elseif i==4
            q=[1 .3 0 .3 1]; th=55; label='sparse_long_heating';
        elseif i==5
            q=.5*ones(1,5); th=60; label='uniform_low_power_long_heating';
        elseif i==6
            mode='preheat'; q=[1 1 .9965 1 1]; th=45;
            label='preheat_center_reduction';
        else
            % Poll each independent power around the best point found so far.
            offset=i-7; cellIndex=mod(offset,5)+1;
            step=.002/(2^floor(offset/10));
            direction=2*mod(floor(offset/5),2)-1;
            q(cellIndex)=min(1,max(0,q(cellIndex)+direction*step));
            label=sprintf('adaptive_cell%d_direction%d',cellIndex,direction);
            if isequal(q,run.bestCo.q)
                run.skippedNoChange=run.skippedNoChange+1;
                continue;
            end
        end
        fprintf('Targeted %d %s q=%s th=%.8f\n',i,label,mat2str(q,9),th);
        c=q3_cached_energy6(mode,q,th,opt);
        accepted=false;
        if strcmp(c.status,'success')
            if strcmp(mode,'preheat'), incumbent=run.bestPre;
            else, incumbent=run.bestCo;
            end
            accepted=c.energyJ<incumbent.energyJ || ...
                (c.energyJ==incumbent.energyJ && ...
                c.startupTimeS<incumbent.startupTimeS);
            if accepted
                if strcmp(mode,'preheat'), run.bestPre=c;
                else, run.bestCo=c;
                end
            end
        end
        save(fullfile(out,sprintf('q3_targeted_%03d.mat',i)), ...
            'c','label','accepted','opt','-v7.3');
        h=struct('index',i,'label',label,'mode',mode,'q',q,'th',th, ...
            'status',c.status,'energyJ',c.energyJ, ...
            'startupTimeS',c.startupTimeS,'accepted',accepted, ...
            'bestCoEnergyJ',run.bestCo.energyJ,'bestPreEnergyJ',run.bestPre.energyJ);
        if isempty(run.history), run.history=h; else, run.history(end+1)=h; end
        run.evaluations=run.evaluations+1;
        save([path '.tmp'],'run','-v7.3'); movefile([path '.tmp'],path,'f');
        writetable(struct2table(run.history),fullfile(out,'q3_targeted_history.csv'), ...
            'Encoding','UTF-8');
        fprintf('Targeted %d: %s E=%.6f t=%.6f accepted=%d\n', ...
            i,c.status,c.energyJ,c.startupTimeS,accepted);
    end
end
