function test_q3_energy_incumbent()
% Repeated near-ties must not walk the reported minimum energy upwards.
    out=fullfile(fileparts(mfilename('fullpath')),'results');
    path=[tempname(out) '.mat']; clean=onCleanup(@() cleanup(path)); %#ok<NASGU>
    calls=0;
    settings=struct('populationSize',8,'evaluator',@drifting_objective, ...
        'archivePath',path);
    run=q3_joint_de_energy(6,settings);
    assert(run.best.energyJ==100,'Reported minimum drifted upward across near-ties.');
    disp('test_q3_energy_incumbent passed');
    function c=drifting_objective(mode,q,th,~)
        calls=calls+1;
        c=struct('mode',mode,'q',q,'th',th,'status','success', ...
            'energyJ',100+.09*(calls-1),'startupTimeS',100-calls, ...
            'violation',0,'result',[],'reason','mock');
    end
end
function cleanup(path)
    if exist(path,'file'), delete(path); end
    if exist([path '.tmp'],'file'), delete([path '.tmp']); end
end
