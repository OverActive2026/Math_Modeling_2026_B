function report = run_q2_constant_fixed_reference_optimization(budget,resumeFile)
% Fresh FuRBO+LogCEI search relative to the successful fixed loading.
codeDir = fileparts(mfilename('fullpath'));
addpath(codeDir,'-begin');
cfg = q2_constant_fixed_reference_config();
if nargin>=1 && ~isempty(budget), cfg.budget = budget; end
if nargin>=2 && ~isempty(resumeFile), cfg.resumeFile = resumeFile; end
if ~exist(cfg.outputDir,'dir'), mkdir(cfg.outputDir); end
cfg.outputFile = [tempname(cfg.outputDir) '_report.mat'];
fprintf('FIXED_REFERENCE plateau=%.15g A/cm2 startup=%.12g s.\n', ...
    cfg.reference.strategy.plateauAcm2,cfg.referenceTimeS);
fprintf('GLOBAL_SEARCH plateau=[0.02,0.5]; LHS_START=[0.245,0.270]; gammaIce=3.5 kFreeze=0.3.\n');
report = run_q2_optimized(cfg);
if report.search.best_is_observed
    [check,trajectory] = q2_real_blackbox(report.search.best_x,cfg);
    assert(check.success && abs(check.startup_time-report.search.best_time)<0.02, ...
        'Best constant-current candidate did not reproduce.');
    report.verificationPerformed = true;
    report.verified = true;
    report.verifiedResult = trajectory;
    report.referenceTimeS = cfg.referenceTimeS;
    report.improvementS = cfg.referenceTimeS-check.startup_time;
    report.improvementPercent = 100*report.improvementS/cfg.referenceTimeS;
    fprintf('COMPARISON reference=%.9f optimized=%.9f saved=%.9f s (%.4f%%) Vmin=%.9f\n', ...
        cfg.referenceTimeS,check.startup_time,report.improvementS, ...
        report.improvementPercent,check.V_min);
    fprintf('PLATEAU=%.15g A/cm2, charge=%.9f C/cm2.\n', ...
        report.bestStrategy.plateauAcm2,check.q_use);
end
save(cfg.outputFile,'report');
fprintf('VERIFIED_REPORT=%s\n',cfg.outputFile);
end
