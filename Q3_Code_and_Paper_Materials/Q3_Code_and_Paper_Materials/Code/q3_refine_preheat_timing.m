function candidates=q3_refine_preheat_timing()
% Remove the old 0.01-second waiting buffer after the temperature event.
% Retain the base model's temperature event margin; only refine heater time.
    out=fullfile(fileparts(mfilename('fullpath')),'results');
    a=load(fullfile(out,'q3_final_energy_table.mat'),'verified');
    probes={a.verified{1}};
    p=fullfile(out,'q3_targeted_006.mat');
    if exist(p,'file')
        b=load(p,'c'); probes{end+1}=b.c;
    end
    cfg=q3_optimization_config('smoke'); opt=cfg.simOpt;
    opt.RelTol=1e-6; opt.AbsTol=1e-10; opt.MaxStep=.0125;
    candidates=cell(size(probes));
    for k=1:numel(probes)
        probe=probes{k}; r=probe.result;
        if isempty(r) || r.solverTerminatedUnexpectedly, continue; end
        times=r.event.time(r.event.index==1);
        if isempty(times), continue; end
        th=times(1)+1e-5;
        savedPath=fullfile(out,sprintf('q3_preheat_timing_%02d.mat',k));
        old=[];
        if exist(savedPath,'file'), old=load(savedPath,'c','opt'); end
        if ~isempty(old) && isequaln(old.opt,opt) && isequal(old.c.q,probe.q)
            c=old.c;
        else
            c=q3_evaluate_energy6('preheat',probe.q,th,opt);
        end
        attemptHistory={c};
        for retry=1:2
            if ~strcmp(c.status,'physical_infeasible') || isempty(c.result), break; end
            rr=c.result;
            deficit=opt.temperatureEventMarginK-min(rr.TavgC(end,:));
            if deficit<0 || deficit>.01 || numel(rr.t)<2, break; end
            slope=(rr.TavgC(end,:)-rr.TavgC(end-1,:))/(rr.t(end)-rr.t(end-1));
            if any(slope<=0), break; end
            % Correct the measured shortfall; this only proposes a duration.
            % Full simulation, not extrapolation, determines acceptance.
            extra=max(1e-4,2*max(max(0,opt.temperatureEventMarginK- ...
                rr.TavgC(end,:))./slope));
            th=c.th+extra;
            c=q3_evaluate_energy6('preheat',probe.q,th,opt);
            attemptHistory{end+1}=c; %#ok<AGROW>
        end
        candidates{k}=c;
        save(savedPath,'c','opt','attemptHistory','-v7.3');
        fprintf('Preheat event refinement q=%s: %s E=%.6f t=%.8f\n', ...
            mat2str(c.q,7),c.status,c.energyJ,c.startupTimeS);
    end
end
