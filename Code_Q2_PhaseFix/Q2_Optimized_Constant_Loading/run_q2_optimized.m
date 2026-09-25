function report = run_q2_optimized(cfg)
% report = run_q2_optimized(); 或传入修改后的配置。
codeDir = fileparts(mfilename('fullpath'));
addpath(codeDir,'-begin');
assert(strcmp(fileparts(which('pemfc_stack5_simulate')),codeDir), ...
    'Change MATLAB current folder to Code_Q2_Optimized to avoid old functions.');
if nargin<1, cfg = q2_optimization_config(); end
switch lower(cfg.strategyType)
    case 'constant', d = 1;
    case 'linear', d = 3;
    case 'step'
        stages = 3;
        if isfield(cfg,'stepStageCount'), stages = cfg.stepStageCount; end
        d = 2*stages-1;
    otherwise, error('Unknown strategyType.');
end
if ~exist(cfg.outputDir,'dir'), mkdir(cfg.outputDir); end
filePath = [tempname(cfg.outputDir) '.mat'];
if isfield(cfg,'outputFile') && ~isempty(cfg.outputFile)
    filePath = cfg.outputFile;
end
resumeFile = '';
if isfield(cfg,'resumeFile') && ~isempty(cfg.resumeFile)
    resumeFile = cfg.resumeFile;
    previous = load(resumeFile,'report');
    assert(isfield(previous,'report'), ...
        'Resume source must be a completed optimization report.');
    old = previous.report.config;
    oldStages = 3; newStages = 3;
    if isfield(old,'stepStageCount'), oldStages = old.stepStageCount; end
    if isfield(cfg,'stepStageCount'), newStages = cfg.stepStageCount; end
    assert(oldStages==newStages,'Resume source has a different number of stages.');
    if strcmpi(cfg.strategyType,'constant')
        assert(isequaln(old.constantPlateauBounds,cfg.constantPlateauBounds) && ...
            old.constantInitialCurrentAcm2==cfg.constantInitialCurrentAcm2 && ...
            old.constantRampTimeS==cfg.constantRampTimeS, ...
            'Resume source uses different constant-current search bounds.');
    end
    if isfield(cfg,'experiment')
        assert(isfield(old,'experiment') && ...
            strcmp(old.experiment,cfg.experiment), ...
            'Resume source belongs to another experiment.');
    end
    if isfield(cfg,'initialDesignBounds')
        assert(isfield(old,'initialDesignBounds') && ...
            isequaln(old.initialDesignBounds,cfg.initialDesignBounds), ...
            'Resume source uses a different Latin-hypercube initial range.');
    end
    assert(isfield(old,'integrationScheme') && ...
        strcmp(old.integrationScheme,cfg.integrationScheme), ...
        'q2:IntegrationMismatch', ...
        '积分方式不同：旧分段积分结果不能混入恢复原积分后的搜索。');
    assert(strcmpi(old.strategyType,cfg.strategyType) && ...
        old.initialTemperatureC==cfg.initialTemperatureC && ...
        old.seed==cfg.seed && old.nInit==cfg.nInit && ...
        isequaln(old.simOpt,cfg.simOpt) && ...
        old.currentMax==cfg.currentMax && ...
        old.switchMinS==cfg.switchMinS && ...
        old.switchMaxS==cfg.switchMaxS && ...
        old.switchGapS==cfg.switchGapS && ...
        isequaln(old.initialNormalizedPoint,cfg.initialNormalizedPoint), ...
        'Resume source uses different physical or search settings.');
    assert(cfg.budget>previous.report.search.nfe, ...
        'New budget must exceed completed evaluations.');
end
assert(~strcmp(filePath,resumeFile), ...
    'Output must not overwrite the previous report.');
report.config = cfg;
report.scope = 'Coarse-grid FuRBO trust region plus analytic LogCEI search.';
initialPoint = [];
if isfield(cfg,'initialNormalizedPoint')
    initialPoint = cfg.initialNormalizedPoint;
end
designBounds = [];
if isfield(cfg,'initialDesignBounds')
    designBounds = cfg.initialDesignBounds;
end
report.search = q2_furbo_logcei_optimize(@(x) q2_real_blackbox(x,cfg), ...
    d,cfg.budget,cfg.nInit,cfg.seed,initialPoint,filePath,resumeFile,designBounds);
report.bestStrategy = [];
report.verified = false;
report.verificationPerformed = false;
report.resultFile = filePath;
if report.search.best_is_observed
    report.bestStrategy = q2_decode_parameters(report.search.best_x,cfg);
end
save(filePath,'report'); % 保存粗网格真实评价结果。
if ~report.search.best_is_observed
    fprintf('No feasible strategy observed within this budget.\n');
else
    fprintf('Best observed coarse-grid startup time=%g s after %d evaluations.\n', ...
        report.search.best_time,report.search.nfe);
end
fprintf('Saved: %s\n',filePath);
end
