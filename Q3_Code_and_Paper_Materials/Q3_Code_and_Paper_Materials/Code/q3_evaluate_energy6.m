function candidate = q3_evaluate_energy6(mode,q,th,opt,simulator)
% Evaluate all six decision variables with an independent observation horizon.
% The coheat horizon follows the 20 C/cm2 charge cap, not heater energy.
    if nargin < 4 || isempty(opt)
        cfg = q3_optimization_config('smoke');
        opt = cfg.simOpt;
    end
    if nargin<5, simulator=@pemfc_stack5_simulate_q3; end
    q = double(q(:).');
    assert(numel(q)==5 && all(isfinite(q)) && all(q>=0 & q<=1));
    assert(isscalar(th) && isfinite(th) && th>=0);
    mode = char(mode);
    if strcmp(mode,'coheat')
        % Q(60)=9 C/cm2, then dQ/dt=0.3 C/(cm2 s).
        opt.tMax = 60+(20-9)/0.3+1e-4;
    elseif strcmp(mode,'preheat')
        opt.tMax = max(th,1e-6);
    else
        error('q3_evaluate_energy6:Mode','Unknown startup mode.');
    end
    opt.enforceChargeCap = true;
    opt.stopOnSuccess = true;
    opt.plot = false;
    opt.verbose = false;
    candidate = struct('mode',mode,'q',q,'th',th, ...
        'status','numerical_unknown','energyJ',Inf,'startupTimeS',NaN, ...
        'violation',Inf,'minimumVoltageV',NaN,'maxIce',NaN, ...
        'chargeCcm2',NaN,'reason','','result',[], ...
        'attempts',0,'elapsedS',0,'observationEndS',opt.tMax);
    tick = tic;
    for attempt = 1:2
        candidate.attempts = attempt;
        if attempt == 2
            opt.RelTol = min(opt.RelTol/10,1e-5);
            opt.AbsTol = min(opt.AbsTol/10,1e-9);
            opt.MaxStep = min(opt.MaxStep/2,0.025);
        end
        try
            r = simulator(mode,q,th,opt);
            candidate.result = r;
            candidate.energyJ = r.totalAuxiliaryEnergyJ;
            candidate.startupTimeS = r.successTimeS;
            candidate.minimumVoltageV = r.minimumCellVoltageV;
            candidate.maxIce = r.maximumIceVolumeFraction;
            candidate.chargeCcm2 = r.chargeUsedCcm2;
            candidate.reason = char(r.stopReason);
            if r.solverTerminatedUnexpectedly
                continue;
            elseif r.success
                candidate.status = 'success';
                candidate.violation = 0;
                break;
            else
                candidate.status = 'physical_infeasible';
                candidate.violation = 1+ ...
                    max(0,0.01-min(r.TavgC(end,:)))/30 + ...
                    max(0,0.3-r.minimumCellVoltageV)/0.1 + ...
                    max(0,r.maximumIceVolumeFraction-0.99)/0.1;
                break;
            end
        catch ME
            candidate.reason = [ME.identifier ': ' ME.message];
        end
    end
    candidate.elapsedS = toc(tick);
end
