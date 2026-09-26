function result=verify_q2_optimized(stage,coldTemperatureC)
if nargin<1, stage=7; end
if nargin<2
    if stage==7, coldTemperatureC=-11.2; else, coldTemperatureC=-11; end
end
folder=fileparts(mfilename('fullpath')); addpath(folder,'-begin');
a=load(fullfile(folder,sprintf('corrected_report_%d_T%g.mat',stage,coldTemperatureC)),'report');
assert(~isempty(a.report.best),'No feasible observed candidate.');
opt=a.report.simOpt; opt.tMax=600; opt.coldTMaxS=600;
% Independent, uncapped recomputation: no cache lookup or early objective cutoff.
clock=tic;
result=q2_evaluate_schedule(q2_step_strategy(a.report.bestX), ...
    a.report.temperaturesC,opt,600);
wallSeconds=toc(clock);
save(fullfile(folder,sprintf('corrected_verification_%d_T%g.mat',stage,coldTemperatureC)), ...
    'result','wallSeconds','-v7.3');
fprintf('INDEPENDENT joint=%d warm=%.10f s wall=%.1f s\n', ...
    result.success,result.startup_time,wallSeconds);
end
