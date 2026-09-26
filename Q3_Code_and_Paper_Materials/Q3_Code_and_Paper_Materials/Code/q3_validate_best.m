function audit = q3_validate_best(run,fullGrid,tightTolerance)
% Re-evaluate an archived candidate using the task's startup criteria.
    if nargin < 2, fullGrid = false; end
    if nargin < 3, tightTolerance = false; end
    if isempty(run.best) || ~strcmp(run.best.status,'success')
        error('q3_validate_best:NoFeasible','No feasible candidate to validate.');
    end
    opt = run.cfg.simOpt;
    opt.plot = false;
    opt.verbose = false;
    if tightTolerance
        opt.RelTol = min(opt.RelTol/10,1e-6);
        opt.MaxStep = opt.MaxStep/2;
    end
    if fullGrid
        opt.nCellLayer = [6 13 19 19 6];
        label = 'full_grid';
    else
        label = 'same_coarse_grid';
    end
    th = run.best.thRequestedS;
    if strcmp(run.mode,'preheat')
        opt.tMax = th;
    else
        opt.tMax = max(run.cfg.coheatObservationS,th);
    end
    r = pemfc_stack5_simulate_q3(run.mode,run.best.q,th,opt);
    audit = struct('grid',label,'mode',run.mode,'q',run.best.q, ...
        'thRequestedS',th,'searchEnergyJ',run.best.energyJ, ...
        'verifiedEnergyJ',r.totalAuxiliaryEnergyJ, ...
        'searchSuccessTimeS',run.best.successTimeS, ...
        'verifiedSuccessTimeS',r.successTimeS, ...
        'questionSuccess',r.questionSuccess, ...
        'minimumVoltageV',r.minimumCellVoltageV, ...
        'maximumIceVolumeFraction',r.maximumIceVolumeFraction, ...
        'result',r);
    outdir = fullfile(fileparts(mfilename('fullpath')),'results');
    if ~exist(outdir,'dir'), mkdir(outdir); end
    save(fullfile(outdir,sprintf('q3_%s_%s_%s_validation.mat', ...
        run.mode,run.cfg.acceptanceMode,label)),'audit','-v7.3');
    fprintf('%s %s: success=%d, E=%.3f J, t=%.3f s, Vmin=%.4f V\n', ...
        run.mode,label,audit.questionSuccess, ...
        audit.verifiedEnergyJ,audit.verifiedSuccessTimeS, ...
        audit.minimumVoltageV);
end
