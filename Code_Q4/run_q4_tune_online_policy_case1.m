function study = run_q4_tune_online_policy_case1()
%RUN_Q4_TUNE_ONLINE_POLICY_CASE1 Screen and verify state-feedback settings.
% All candidates use measured T, V and slopes, never a preset switch time.
% Stage 1 ranks candidate taper rules on a cheap grid. Stage 2 evaluates
% the best rules on the same [2 3 4 4 2] grid as the published comparison.
% This is parameter search within one rule family, not a global optimum.

    codeDir = fileparts(mfilename('fullpath'));
    addpath(codeDir);
    outputDir = fullfile(codeDir,'results_control');
    if ~exist(outputDir,'dir'), mkdir(outputDir); end
    wallClock = tic;

    base.sampleS = 0.5;
    base.centerEndOffC = -1.20;
    base.centerTaperBandC = 0.40;
    base.centerWarmReserveC = 14;
    base.adjacentEndOffC = -0.12;
    base.adjacentTaperBandC = 0.30;
    base.adjacentWarmReserveC = 11;
    base.minimumEndRiseCps = 0.35;
    base.voltageWarningV = 0.45;
    base.voltageLookaheadS = 1;
    base.iceRiskVoltageV = 0.55;
    base.iceRiskDropVps = 0.02;

    % Columns: center shutoff threshold, center taper width, adjacent
    % shutoff threshold, adjacent taper width. Include the current rule.
    candidates = [ ...
        -1.20 0.40 -0.12 0.30; ... % current baseline
        -1.20 0.05 -0.22 0.05; ... % near historical switch temperatures
        -1.20 0.15 -0.22 0.05; ...
        -1.20 0.25 -0.22 0.05; ...
        -1.20 0.05 -0.22 0.15; ...
        -1.20 0.05 -0.22 0.30; ...
        -1.20 0.15 -0.22 0.15; ...
        -1.30 0.05 -0.22 0.05; ...
        -1.10 0.05 -0.22 0.05; ...
        -1.20 0.05 -0.30 0.05; ...
        -1.20 0.05 -0.14 0.05; ...
        -1.20 0.15 -0.14 0.15];
    n = size(candidates,1);
    coarse = nan(n,5); % success, time, energy, Vmin, iceMax
    tune = nan(n,5);
    tuneResults = cell(n,1);

    simOpt.nCellLayer = [1 1 1 1 1];
    simOpt.kFreeze = 0.4;
    simOpt.gammaIce = 3.5;
    simOpt.kMelt = 1;
    simOpt.MaxStep = 0.2;
    simOpt.RelTol = 1e-3;
    simOpt.AbsTol = 1e-7;
    simOpt.successMarginC = 0.01;
    simOpt.recordStepS = 0;
    simOpt.progressEveryS = 0;
    fprintf('Screening %d online rules on the coarse grid...\n',n);
    model = [];
    for j = 1:n
        control = set_candidate(base,candidates(j,:));
        [r,model] = run_rule(control,simOpt,model);
        coarse(j,:) = summary_row(r);
        fprintf(['%2d/%d cOff=%+.2f cBand=%.2f aOff=%+.2f ', ...
            'aBand=%.2f: success=%d, t=%.3f s, E=%.3f J\n'], ...
            j,n,candidates(j,:),r.success,r.stopTimeS, ...
            r.totalAuxiliaryEnergyJ);
    end

    valid = find(coarse(:,1) == 1);
    [~,rank] = sort(coarse(valid,3));
    sorted = valid(rank);
    finalists = unique([1;sorted(1:min(4,numel(sorted)))],'stable');
    simOpt.nCellLayer = [2 3 4 4 2];
    simOpt.MaxStep = 0.05;
    simOpt.RelTol = 1e-4;
    simOpt.AbsTol = 1e-8;
    simOpt.recordStepS = 0.5;
    model = [];
    fprintf('Verifying %d candidates on [2 3 4 4 2]...\n', ...
        numel(finalists));
    for k = 1:numel(finalists)
        j = finalists(k);
        control = set_candidate(base,candidates(j,:));
        [r,model] = run_rule(control,simOpt,model);
        tune(j,:) = summary_row(r);
        tuneResults{j} = r;
        fprintf(['candidate %d: success=%d, t=%.3f s, E=%.3f J, ', ...
            'Vmin=%.4f V, iceMax=%.4f\n'],j,r.success,r.stopTimeS, ...
            r.totalAuxiliaryEnergyJ,r.minimumCellVoltageV, ...
            r.maximumIceVolumeFraction);
    end

    valid = find(tune(:,1) == 1);
    [~,rank] = sort(tune(valid,3));
    bestIndex = valid(rank(1));
    bestControl = set_candidate(base,candidates(bestIndex,:));
    best = tuneResults{bestIndex};
    baseline = tuneResults{1};
    benchmarkFile = fullfile(outputDir,'q4_old_schedule_melt_trial.mat');
    benchmarkEnergyJ = NaN;
    if isfile(benchmarkFile)
        saved = load(benchmarkFile,'r');
        benchmarkEnergyJ = saved.r.totalAuxiliaryEnergyJ;
    end
    study.candidates = candidates;
    study.coarse = coarse;
    study.tune = tune;
    study.finalists = finalists;
    study.bestIndex = bestIndex;
    study.bestControl = bestControl;
    study.best = best;
    study.baseline = baseline;
    study.benchmarkEnergyJ = benchmarkEnergyJ;
    study.simOpt = simOpt;
    study.wallTimeS = toc(wallClock);
    save(fullfile(outputDir,'q4_online_policy_search_case1_melt.mat'), ...
        'study');
    q4_plot_control_trajectory(best,fullfile(outputDir, ...
        'q4_online_policy_search_case1_melt.png'), ...
        'Q4 case 1 - best screened online policy',usejava('desktop'));
    fprintf(['Best online candidate %d: %.3f J / %.3f s; ', ...
        'current rule %.3f J; old fixed schedule %.3f J.\n'], ...
        bestIndex,best.totalAuxiliaryEnergyJ,best.stopTimeS, ...
        baseline.totalAuxiliaryEnergyJ,benchmarkEnergyJ);
    fprintf('Search wall time: %.1f s\n',study.wallTimeS);
end

function control = set_candidate(base,row)
    control = base;
    control.centerEndOffC = row(1);
    control.centerTaperBandC = row(2);
    control.adjacentEndOffC = row(3);
    control.adjacentTaperBandC = row(4);
end

function [r,model] = run_rule(control,simOpt,model)
    schedule.breaksS = [0 80];
    schedule.powerWcm2 = ones(1,5);
    schedule.sampleS = control.sampleS;
    schedule.feedbackFcn = @(~,T,V,dT,dV,~) ...
        q4_online_rule_controller(T,V,dT,dV,control);
    [r,model] = q4_simulate_power_schedule_case1(schedule,simOpt,model);
end

function row = summary_row(r)
    row = [r.success,r.stopTimeS,r.totalAuxiliaryEnergyJ, ...
        r.minimumCellVoltageV,r.maximumIceVolumeFraction];
end
