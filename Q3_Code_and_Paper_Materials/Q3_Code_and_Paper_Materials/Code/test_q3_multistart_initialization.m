function test_q3_multistart_initialization()
% Custom starts must survive initialization; trusted seed records avoid solves.
    out=fullfile(fileparts(mfilename('fullpath')),'results');
    path=[tempname(out) '.mat']; cleaner=onCleanup(@() cleanup(path)); %#ok<NASGU>
    P=[.8 .7 .6 .5 .4 .3;.4 .5 .6 .7 .8 .4; ...
       1 .9 .8 .7 .6 .5;.6 .7 .8 .9 1 .6];
    calls=0; seeds=cell(4,1); cap=60+(20-9)/.3;
    for k=1:4, seeds{k}=record(P(k,1:5),P(k,6)*cap,k); end
    settings=struct('populationSize',4,'initialPopulation',P, ...
        'initialRecords',{seeds},'evaluator',@evaluate,'archivePath',path);
    run=q3_joint_de_energy(4,settings);
    assert(isequal(run.population,P),'Custom initial population was overwritten.');
    assert(calls==0,'A saved seed record was needlessly simulated again.');
    run=q3_joint_de_energy(4,settings);
    assert(calls==4 && run.generation==1,'New DE trials must use the evaluator.');
    disp('test_q3_multistart_initialization passed');
    function c=evaluate(~,q,th,~)
        calls=calls+1;c=record(q,th,10+calls);
    end
end
function c=record(q,th,k)
    c=struct('mode','coheat','q',q,'th',th,'status','success', ...
        'energyJ',100+k,'startupTimeS',th,'violation',0,'result',[],'reason','mock');
end
function cleanup(path)
    if exist(path,'file'), delete(path); end
    if exist([path '.tmp'],'file'), delete([path '.tmp']); end
end
