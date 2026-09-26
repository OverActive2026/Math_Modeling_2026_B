function c=q3_small_allocation_probe(index)
% Resolve narrow feasible directions when only one end cell is active.
    Q=[1 1 1 1 .998;1 1 1 .998 1;1 1 1 .999 .999];
    cfg=q3_optimization_config('smoke');
    opt=cfg.simOpt; opt.RelTol=1e-5; opt.AbsTol=1e-9; opt.MaxStep=.025;
    c=q3_evaluate_energy6('coheat',Q(index,:),40,opt);
    outDir=fullfile(fileparts(mfilename('fullpath')),'results');
    save(fullfile(outDir,sprintf('q3_small_allocation_%02d.mat',index)),'c','-v7.3');
    fprintf('Small allocation %d q=%s: %s E=%.6f J t=%.6f s\n', ...
        index,mat2str(c.q,6),c.status,c.energyJ,c.startupTimeS);
end
