function summary = trial_q3_freeze_parameters()
% Exploratory coarse-grid check for the requested ice parameters.
% These runs check numerical connectivity, not parameter calibration.

    opt = struct('initialTemperatureC',-30, ...
        'ambientTemperatureC',-30,'nCellLayer',[2 3 4 4 2], ...
        'tMax',180,'outputStep',0.5,'MaxStep',0.05, ...
        'RelTol',1e-4,'AbsTol',1e-8,'plot',false, ...
        'verbose',false,'progressIntervalS',0, ...
        'useJacobianPattern',true,'gammaIce',3.5);

    mode = {'preheat';'coheat';'coheat';'coheat'};
    kFreeze = [0.4;0.3;0.4;0.5];
    power = {[0.90 0.70 0.65 0.70 0.90]; ...
        [0.55 0.35 0.30 0.35 0.55]; ...
        [0.55 0.35 0.30 0.35 0.55]; ...
        [0.55 0.35 0.30 0.35 0.55]};
    durationS = [120;80;80;80];
    n = numel(mode);
    completed = false(n,1);
    success = false(n,1);
    stopTimeS = nan(n,1);
    energyJ = nan(n,1);
    minVoltageV = nan(n,1);
    maxIceVolumeFraction = nan(n,1);
    maxPoreOccupancy = nan(n,1);
    elapsedS = nan(n,1);
    finiteTrajectory = false(n,1);
    reason = strings(n,1);

    for i = 1:n
        opt.kFreeze = kFreeze(i);
        fprintf('Trial %d/%d: %s, kFreeze=%.2f, gammaIce=%.1f\n', ...
            i,n,mode{i},kFreeze(i),opt.gammaIce);
        timer = tic;
        try
            r = pemfc_stack5_simulate_q3(mode{i},power{i},durationS(i),opt);
            completed(i) = ~r.solverTerminatedUnexpectedly;
            success(i) = r.success;
            stopTimeS(i) = r.stopTimeS;
            energyJ(i) = r.totalAuxiliaryEnergyJ;
            minVoltageV(i) = r.minimumCellVoltageV;
            maxIceVolumeFraction(i) = r.maximumIceVolumeFraction;
            maxPoreOccupancy(i) = r.maximumPoreOccupancy;
            finiteTrajectory(i) = all(isfinite(r.x(:))) && ...
                all(isfinite(r.TavgC(:))) && all(isfinite(r.Vcell(:)));
            reason(i) = string(r.stopReason);
        catch ME
            reason(i) = string(ME.identifier) + ": " + string(ME.message);
        end
        elapsedS(i) = toc(timer);
        fprintf('  completed=%d, success=%d, elapsed=%.1f s, reason=%s\n', ...
            completed(i),success(i),elapsedS(i),reason(i));
    end

    summary = table(string(mode),kFreeze,repmat(opt.gammaIce,n,1), ...
        completed,finiteTrajectory,success,stopTimeS,energyJ, ...
        minVoltageV,maxIceVolumeFraction,maxPoreOccupancy,elapsedS,reason, ...
        'VariableNames',{'mode','kFreeze','gammaIce','completed', ...
        'finiteTrajectory','success','stopTimeS','energyJ','minVoltageV', ...
        'maxIceVolumeFraction','maxPoreOccupancy','elapsedS','reason'});
    writetable(summary,fullfile(fileparts(mfilename('fullpath')), ...
        'q3_freeze_parameter_trial.csv'));
    disp(summary);
end
