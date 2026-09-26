function result=verify_q2_reference(stage)
if nargin<1, stage=7; end
folder=fileparts(mfilename('fullpath')); addpath(folder,'-begin');
if stage==7, cold=-11.2; else, cold=-11; end
clock=tic;
result=q2_evaluate_schedule(q2_reference_strategy(stage),[-10 cold], ...
    q2_verified_options(),600);
wallSeconds=toc(clock);
save(fullfile(folder,sprintf('verified_reference_%d.mat',stage)), ...
    'result','wallSeconds','-v7.3');
fprintf('REFERENCE stage=%d joint=%d warm=%.10f s wall=%.1f s\n', ...
    stage,result.success,result.startup_time,wallSeconds);
assert(~result.success || all(result.status==1));
end
