function out = q2_evaluate_real(x,kind,initialTemperatureC,simOpt)
%Q2_EVALUATE_REAL Evaluate one candidate using the copied five-cell model.
% Failed runs deliberately have startup_time=NaN (never a penalty time).
if nargin<4, simOpt=struct(); end
strategy = q2_decode_strategy(x,kind);
[r,~] = pemfc_stack5_simulate(strategy,initialTemperatureC,simOpt);
out.x = double(x(:).');
out.strategy = strategy;
out.success = logical(r.success);
out.startup_time = r.successTimeS;
out.V_min = r.minimumCellVoltageV;
out.q_use = r.chargeUsedCcm2;
out.ice_peak = r.maximumIceVolumeFraction;
out.ice_start = NaN;
if r.success
    % The last accepted point is the startup event, not the peak in history.
    out.ice_start = max(r.maxIceVolumeFractionCell(end,:));
end
out.T_min_end = 273.15+min(r.TavgC(end,:));
out.critical_cell = r.minimumVoltageCellIndex;
out.failure_reason = string(r.stopReason);
out.termination_time = r.stopTimeS;
out.is_mock = false;
out.solver_failed = logical(r.solverTerminatedUnexpectedly);
out.time_censored = false;
if isfield(simOpt,'optimizationCutoffActive') && ...
        simOpt.optimizationCutoffActive && ~r.success && ...
        ~r.solverTerminatedUnexpectedly && isempty(r.event.index) && ...
        r.stopTimeS >= simOpt.tMax-1e-6
    % The scheme is only known not to beat the current incumbent;
    % it has NOT been established to be physically infeasible.
    out.time_censored = true;
end
out.raw = r;
end
