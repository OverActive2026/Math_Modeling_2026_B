function result = run_q2_constant_best_current(reportFile)
% Reproduce the best verified new-parameter constant-current result.
codeDir = fileparts(mfilename('fullpath'));
addpath(codeDir,'-begin');
assert(strcmp(fileparts(which('pemfc_stack5_simulate')),codeDir));
cfg = q2_constant_fixed_reference_config();
if nargin<1 || isempty(reportFile)
    files = dir(fullfile(cfg.outputDir,'*_report.mat'));
    bestTime = Inf;
    reportFile = '';
    for k = 1:numel(files)
        path = fullfile(files(k).folder,files(k).name);
        data = load(path,'report');
        if isfield(data,'report') && data.report.verified && ...
                strcmp(data.report.config.experiment,cfg.experiment) && ...
                data.report.search.best_time<bestTime
            bestTime = data.report.search.best_time;
            reportFile = path;
        end
    end
    assert(~isempty(reportFile),'No verified report in this experiment.');
end
data = load(reportFile,'report');
report = data.report;
assert(report.verified && strcmp(report.config.experiment,cfg.experiment));
assert(isequaln(report.config.simOpt,cfg.simOpt));
[check,result] = q2_real_blackbox(report.search.best_x,report.config);
assert(check.success && abs(check.startup_time-report.search.best_time)<0.02);
fprintf('Plateau %.12f A/cm2: %.6f s vs reference %.6f s, saved %.3f%%; Vmin %.6f V\n', ...
    result.strategy.plateauAcm2,result.successTimeS,cfg.referenceTimeS, ...
    100*(cfg.referenceTimeS-result.successTimeS)/cfg.referenceTimeS, ...
    result.minimumCellVoltageV);
report.verifiedResult = result;
plot_q2_constant_fixed_comparison(report);
end
