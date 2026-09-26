function reproduce_selected()
% OPTIONAL: independent ODE reruns, potentially slow. Writes new files only.
    root=fileparts(fileparts(mfilename('fullpath')));
    src=fullfile(root,'Code'); addpath(src,'-begin');
    assert(strcmp(which('pemfc_stack5_simulate_q3'),fullfile(src,'pemfc_stack5_simulate_q3.m')));
    out=fullfile(root,'独立复算'); if ~exist(out,'dir'), mkdir(out); end
    for k=1:3
        s=load(fullfile(root,'结果与图表',sprintf('q3_strict_verification_%02d.mat',k)),'checkpoint','opt');
        original=s.checkpoint;
        repeated=q3_evaluate_energy6(original.mode,original.q,original.th,s.opt);
        save(fullfile(out,sprintf('repeat_%02d.mat',k)),'original','repeated','-v7.3');
        fprintf('%s: %s, deltaE=%.9g J, deltaTime=%.9g s\n',original.mode,repeated.status,repeated.energyJ-original.energyJ,repeated.startupTimeS-original.startupTimeS);
    end
end
