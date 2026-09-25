% Expand the initial-current search to the problem's 0.50 A/cm^2 limit.
% Preserve the old strict checkpoint; compatible real trajectories are
% remapped by their physical current strategy, not their old coordinates.
codeDir=fileparts(mfilename('fullpath'));
addpath(codeDir);
opt=struct();
opt.simOpt=q2_optimization_options();
opt.simOpt.RelTol=5e-5;
opt.simOpt.AbsTol=5e-9;
opt.nLHSNew=0;
opt.maxNewEvaluations=6;
opt.initialRadius=0.12;
opt.minRadius=0.015; % Keep searching beyond a collapsed local boundary.
opt.minFeasibilityProbability=0.25;
opt.warmStartCacheFile=fullfile(codeDir,'q2_furbo_strict_checkpoint.mat');
opt.checkpointFile=fullfile(codeDir,'q2_furbo_expanded_strict_checkpoint.mat');
report=q2_optimize_furbo_real(opt);
save(fullfile(codeDir,'q2_furbo_expanded_strict_latest_report.mat'), ...
    'report','-v7.3');
fprintf('Best expanded-search startup time: %.9f s\n',report.bestTimeS);
fprintf('New real-model calls: %d; successes: %d/%d\n', ...
    report.nNewModelCalls,report.nSuccess,report.nEvaluations);
disp(report.bestStrategy);
