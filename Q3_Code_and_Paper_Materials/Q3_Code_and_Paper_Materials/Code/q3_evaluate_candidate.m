function record = q3_evaluate_candidate(mode,z,cfg)
% Evaluate six independent normalized variables [q1..q5,t_h/t_h,max].
    z = double(z(:).');
    if numel(z) ~= 6 || any(~isfinite(z)) || any(z < 0) || any(z > 1)
        error('q3_evaluate_candidate:Bounds','z must be six numbers in [0,1].');
    end
    mode = char(mode);
    q = z(1:5);
    th = cfg.timeUpperS*z(6);
    opt = cfg.simOpt;
    if strcmp(mode,'preheat')
        opt.tMax = max(th,1e-6);
    elseif strcmp(mode,'coheat')
        opt.tMax = cfg.coheatObservationS;
    else
        error('q3_evaluate_candidate:Mode','Unknown startup mode.');
    end
    record = struct('mode',mode,'z',z,'q',q,'thRequestedS',th, ...
        'status','numerical_failure','energyJ',inf,'violation',inf, ...
        'safetyMargin',-inf,'successTimeS',inf,'result',[], ...
        'failureReason','');
    if strcmp(mode,'preheat') && th == 0
        record.status = 'infeasible';
        record.energyJ = 0;
        record.violation = 10;
        record.failureReason = 'Zero preheat duration at -30 C';
        return;
    end
    try
        r = pemfc_stack5_simulate_q3(mode,q,th,opt);
        record.result = r;
        record.energyJ = r.totalAuxiliaryEnergyJ;
        record.successTimeS = r.successTimeS;
        record.safetyMargin = min(r.minimumCellVoltageV-0.30, ...
            0.98-r.maximumPoreOccupancy);
        record.failureReason = r.stopReason;
        finite = all(isfinite(r.x(:))) && all(isfinite(r.Vcell(:))) && ...
            all(isfinite(r.TavgC(:))) && isfinite(record.energyJ);
        if finite && ~r.solverTerminatedUnexpectedly
            if r.success
                record.status = 'success';
                record.violation = 0;
            else
                record.status = 'infeasible';
                deficit = max(0,0.01-min(r.TavgC(end,:)))/30;
                voltageDeficit = max(0,0.30-r.minimumCellVoltageV)/0.10;
                iceDeficit = max(0,r.maximumIceVolumeFraction-0.99)/0.10;
                record.violation = 1+deficit+voltageDeficit+iceDeficit;
            end
        end
    catch ME
        record.failureReason = [ME.identifier ': ' ME.message];
    end
end
