function summary=run_q3_energy_optimization(additionalEvaluations)
% Continue the archived Q3 joint energy search and regenerate verified tables.
    if nargin<1, additionalEvaluations=24; end
    codeDir=fileparts(mfilename('fullpath'));
    addpath(codeDir);
    outDir=fullfile(codeDir,'results');
    if ~exist(fullfile(outDir,'q3_min_energy_preheat.mat'),'file')
        q3_refine_min_energy('preheat',16);
    end
    settings=struct('parallel',false,'batchSize',1);
    finalPath=fullfile(outDir,'q3_final_energy_table.mat');
    if exist(finalPath,'file')
        verified=load(finalPath,'verified');
        settings.eliteCandidate=verified.verified{2};
    end
    run=q3_joint_de_energy(additionalEvaluations,settings);
    assert(~isempty(run.best) && strcmp(run.best.status,'success'), ...
        'Search has not found a feasible coheat candidate.');
    summary=q3_finalize_joint_results();
end
