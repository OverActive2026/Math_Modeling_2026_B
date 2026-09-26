function candidates=q3_refine_end_margin()
% Use measured end-temperature slack to propose energy-saving allocations.
% Every sensitivity-guided proposal is evaluated with the full physical model.
    outDir=fullfile(fileparts(mfilename('fullpath')),'results');
    reference=load(fullfile(outDir,'q3_final_energy_table.mat'));
    if isfield(reference,'audits'), r0=reference.audits{2};
    elseif numel(reference.verified)>=3, r0=reference.verified{3}.result;
    else, r0=reference.verified{2}.result;
    end
    D0=r0.TavgC(end,5)-r0.TavgC(end,1);
    cfg=q3_optimization_config('smoke'); opt=cfg.simOpt;
    opt.RelTol=1e-5; opt.AbsTol=1e-9; opt.MaxStep=.025;
    candidates=cell(2,1);
    for i=1:2
        sample=load(fullfile(outDir,sprintf('q3_small_allocation_%02d.mat',i)),'c');
        assert(strcmp(sample.c.status,'success'));
        D1=sample.c.result.TavgC(end,5)-sample.c.result.TavgC(end,1);
        changedCell=6-i;
        sensitivity=(D1-D0)/(sample.c.q(changedCell)-1);
        if abs(sensitivity)<1e-9, continue; end
        q=ones(1,5); q(changedCell)=min(1,max(.98,1-D0/sensitivity));
        c=q3_evaluate_energy6('coheat',q,40,opt);
        candidates{i}=c;
        save(fullfile(outDir,sprintf('q3_end_margin_%02d.mat',i)), ...
            'c','D0','D1','sensitivity','-v7.3');
        fprintf('End-margin allocation q=%s: %s E=%.6f J t=%.6f s\n', ...
            mat2str(q,8),c.status,c.energyJ,c.startupTimeS);
    end
end
