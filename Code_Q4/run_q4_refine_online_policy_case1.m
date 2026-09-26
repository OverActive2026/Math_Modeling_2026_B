function refinement = run_q4_refine_online_policy_case1()
%RUN_Q4_REFINE_ONLINE_POLICY_CASE1 Refine the last ~2 s of online heating.
% Requires the saved broad search; compares every finalist on the same
% melting model/grid and retains the best genuinely state-driven policy.

    codeDir = fileparts(mfilename('fullpath'));
    addpath(codeDir);
    outputDir = fullfile(codeDir,'results_control');
    previous = load(fullfile(outputDir, ...
        'q4_online_policy_search_case1_melt.mat'),'study');
    base = previous.study.bestControl;
    oldBest = previous.study.best;
    % Columns: center off, center band, adjacent off, adjacent band,
    % feedback sampling period. Row 1 is the already verified baseline.
    candidates = [ ...
        -1.10 0.05 -0.22 0.05 0.50; ...
        -1.10 0.05 -0.12 0.30 0.50; ...
        -1.10 0.05 -0.18 0.30 0.50; ...
        -1.10 0.05 -0.22 0.40 0.50; ...
        -1.12 0.05 -0.22 0.05 0.50; ...
        -1.08 0.05 -0.22 0.05 0.50; ...
        -1.10 0.05 -0.22 0.05 0.25];
    n = size(candidates,1);
    coarse = nan(n,5);
    verified = nan(n,5);
    verified(1,:) = summary_row(oldBest);
    trajectories = cell(n,1);
    trajectories{1} = oldBest;
    timer = tic;

    opt = previous.study.simOpt;
    opt.nCellLayer = [1 1 1 1 1];
    opt.MaxStep = 0.2;
    opt.RelTol = 1e-3;
    opt.AbsTol = 1e-7;
    opt.recordStepS = 0;
    opt.progressEveryS = 0;
    model = [];
    fprintf('Local online-policy screening (%d rules)...\n',n);
    for j = 1:n
        c = set_candidate(base,candidates(j,:));
        [r,model] = run_rule(c,opt,model);
        coarse(j,:) = summary_row(r);
        fprintf('%d/%d: dt=%.2f, E=%.3f J, t=%.3f s, success=%d\n', ...
            j,n,c.sampleS,r.totalAuxiliaryEnergyJ,r.stopTimeS,r.success);
    end

    valid = find(coarse(:,1) == 1 & (1:n).' ~= 1);
    [~,rank] = sort(coarse(valid,3));
    % Always verify the faster-sampling candidate as a separate check.
    finalists = unique([valid(rank(1:min(2,numel(rank))));n]);
    opt.nCellLayer = [2 3 4 4 2];
    opt.MaxStep = 0.05;
    opt.RelTol = 1e-4;
    opt.AbsTol = 1e-8;
    opt.recordStepS = 0.5;
    model = [];
    for k = 1:numel(finalists)
        j = finalists(k);
        c = set_candidate(base,candidates(j,:));
        [r,model] = run_rule(c,opt,model);
        verified(j,:) = summary_row(r);
        trajectories{j} = r;
        fprintf(['Verified %d: E=%.3f J, t=%.3f s, success=%d, ', ...
            'Vmin=%.4f V, iceMax=%.4f\n'],j,r.totalAuxiliaryEnergyJ, ...
            r.stopTimeS,r.success,r.minimumCellVoltageV, ...
            r.maximumIceVolumeFraction);
    end
    feasible = find(verified(:,1) == 1);
    [~,rank] = sort(verified(feasible,3));
    bestIndex = feasible(rank(1));
    refinement.candidates = candidates;
    refinement.coarse = coarse;
    refinement.verified = verified;
    refinement.finalists = finalists;
    refinement.bestIndex = bestIndex;
    refinement.bestControl = set_candidate(base,candidates(bestIndex,:));
    refinement.best = trajectories{bestIndex};
    refinement.benchmarkEnergyJ = previous.study.benchmarkEnergyJ;
    refinement.simOpt = opt;
    refinement.wallTimeS = toc(timer);
    save(fullfile(outputDir,'q4_online_policy_refinement_case1_melt.mat'), ...
        'refinement');
    q4_plot_control_trajectory(refinement.best,fullfile(outputDir, ...
        'q4_online_policy_refinement_case1_melt.png'), ...
        'Q4 case 1 - refined online policy',usejava('desktop'));
    fprintf(['Final online candidate %d: %.3f J / %.3f s; ', ...
        'old fixed schedule %.3f J; search took %.1f s.\n'], ...
        bestIndex,refinement.best.totalAuxiliaryEnergyJ, ...
        refinement.best.stopTimeS,refinement.benchmarkEnergyJ, ...
        refinement.wallTimeS);
end

function c = set_candidate(base,row)
    c = base;
    c.centerEndOffC = row(1);
    c.centerTaperBandC = row(2);
    c.adjacentEndOffC = row(3);
    c.adjacentTaperBandC = row(4);
    c.sampleS = row(5);
end

function [r,model] = run_rule(c,opt,model)
    schedule.breaksS = [0 80];
    schedule.powerWcm2 = ones(1,5);
    schedule.sampleS = c.sampleS;
    schedule.feedbackFcn = @(~,T,V,dT,dV,~) ...
        q4_online_rule_controller(T,V,dT,dV,c);
    [r,model] = q4_simulate_power_schedule_case1(schedule,opt,model);
end

function row = summary_row(r)
    row = [r.success,r.stopTimeS,r.totalAuxiliaryEnergyJ, ...
        r.minimumCellVoltageV,r.maximumIceVolumeFraction];
end
