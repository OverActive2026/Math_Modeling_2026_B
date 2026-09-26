function test_q3_joint_de_energy()
% Cheap deterministic objectives test optimizer state, not the physical model.
    root = fullfile(fileparts(mfilename('fullpath')),'results');
    if ~exist(root,'dir'), mkdir(root); end
    paths = {[tempname(root) '.mat'],[tempname(root) '.mat'], ...
        [tempname(root) '.mat'],[tempname(root) '.mat']};
    clean = onCleanup(@() delete_archives(paths)); %#ok<NASGU>
    settings = struct('populationSize',8,'seed',7321, ...
        'evaluator',@mock_objective,'archivePath',paths{1});
    uninterrupted = q3_joint_de_energy(29,settings);
    settings.archivePath = paths{2};
    q3_joint_de_energy(7,settings);
    q3_joint_de_energy(12,settings);
    resumed = q3_joint_de_energy(10,settings);
    assert(resumed.evaluations == 29);
    assert(isequal(uninterrupted.population,resumed.population), ...
        'Resuming changed the DE population.');
    assert(isequal(uninterrupted.rngState,resumed.rngState));
    assert(isequal(vertcat(uninterrupted.history.z), ...
        vertcat(resumed.history.z)), ...
        'Resuming changed the sequence of trial vectors.');
    assert(isequal(uninterrupted.best.q,resumed.best.q) && ...
        uninterrupted.best.th == resumed.best.th);
    settings.archivePath=paths{4};
    settings.batchSize=3;
    batched=q3_joint_de_energy(29,settings);
    assert(isequal(uninterrupted.population,batched.population) && ...
        isequal(uninterrupted.rngState,batched.rngState), ...
        'Batching evaluations changed frozen-generation DE results.');

    history = resumed.history;
    population = vertcat(history(1:8).z);
    parentGeneration = 0;
    parentPopulation = [];
    movedDuration = false;
    acceptedTrials = 0;
    for i = 9:numel(history)
        h = history(i);
        if h.generation ~= parentGeneration
            parentGeneration = h.generation;
            parentPopulation = population;
        end
        assert(isequal(h.donorZ,parentPopulation(h.parentIndices,:)), ...
            'DE mutation used an updated member within the same generation.');
        assert(isequal(h.targetZ,parentPopulation(h.populationIndex,:)));
        expectedMutant = min(1,max(0,h.donorZ(1,:) + ...
            0.65*(h.donorZ(2,:)-h.donorZ(3,:))));
        expectedTrial = h.targetZ;
        expectedTrial(h.crossoverMask) = expectedMutant(h.crossoverMask);
        assert(isequal(h.z,expectedTrial));
        assert(all(h.z >= 0 & h.z <= 1));
        assert(abs(h.th-h.z(6)*(60+(20-9)/0.3)) < 1e-12);
        movedDuration = movedDuration || ...
            abs(h.z(6)-h.targetZ(6)) > 1e-12;
        if h.accepted
            acceptedTrials = acceptedTrials+1;
            population(h.populationIndex,:) = h.z;
        end
    end
    assert(movedDuration && acceptedTrials > 0);
    assert(isequal(population,resumed.population));
    stored = load(paths{2},'run');
    assert(stored.run.pendingIndex == resumed.pendingIndex);
    assert(isequal(stored.run.parentPopulation,resumed.parentPopulation));

    settings = struct('populationSize',8,'seed',7321, ...
        'evaluator',@classification_objective,'archivePath',paths{3});
    incomplete = q3_joint_de_energy(2,settings);
    assert(strcmp(incomplete.best.status,'physical_infeasible'), ...
        'A numerical unknown must rank below a physical evaluation.');
    classified = q3_joint_de_energy(3,settings);
    assert(strcmp(classified.best.status,'success'));
    assert(abs(classified.best.energyJ-500) < 1e-10 && ...
        classified.best.startupTimeS == 40, ...
        'The reported energy minimum must not increase for a faster candidate.');
    assert(strcmp(classified.history(1).status,'numerical_unknown'));
    assert(strcmp(classified.history(2).status,'physical_infeasible'));
    changed = settings;
    changed.energyTieToleranceJ = 0.2;
    rejected = false;
    try
        q3_joint_de_energy(0,changed);
    catch ME
        rejected = strcmp(ME.identifier,'q3_joint_de_energy:ArchiveMismatch');
    end
    assert(rejected,'An incompatible checkpoint must not be reused.');
    fprintf('test_q3_joint_de_energy: passed\n');
end

function r = mock_objective(mode,q,th,~)
    z = [q,th/(60+(20-9)/0.3)];
    target = [0.85 0.6 0.25 0.6 0.85 0.24];
    r = result_record(mode,q,th,'success',100+1e3*sum((z-target).^2),th,0);
end

function r = classification_objective(mode,q,th,~)
    if all(q == 1) && th < 25
        r = result_record(mode,q,th,'numerical_unknown',0,1,Inf);
    elseif all(q == 1)
        r = result_record(mode,q,th,'physical_infeasible',0.5,2,1);
    elseif all(abs(q-0.8) < 1e-12)
        r = result_record(mode,q,th,'success',500.05,5,0);
    elseif q(3) < 0.3
        r = result_record(mode,q,th,'success',500,40,0);
    else
        r = result_record(mode,q,th,'success',600,10,0);
    end
end

function r = result_record(mode,q,th,status,energy,time,violation)
    r = struct('mode',mode,'q',q,'th',th,'status',status, ...
        'energyJ',energy,'startupTimeS',time,'violation',violation, ...
        'result',[],'reason','Deterministic optimizer regression objective');
end

function delete_archives(paths)
    for i = 1:numel(paths)
        if exist(paths{i},'file'), delete(paths{i}); end
        if exist([paths{i} '.tmp'],'file'), delete([paths{i} '.tmp']); end
    end
end
