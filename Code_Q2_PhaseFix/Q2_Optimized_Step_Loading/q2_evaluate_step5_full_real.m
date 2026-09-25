function out=q2_evaluate_step5_full_real(x,initialTemperatureC,simOpt)
%Q2_EVALUATE_STEP5_FULL_REAL Evaluate all nine step-curve coordinates.
strategy=q2_step5_full_strategy(x);
[r,~]=pemfc_stack5_simulate(strategy,initialTemperatureC,simOpt);
out.x=double(x(:).');
out.strategy=strategy;
out.success=logical(r.success);
out.startup_time=r.successTimeS;
out.V_min=r.minimumCellVoltageV;
out.q_use=r.chargeUsedCcm2;
out.ice_peak=r.maximumIceVolumeFraction;
out.ice_start=NaN;
if r.success, out.ice_start=max(r.maxIceVolumeFractionCell(end,:)); end
out.T_min_end=273.15+min(r.TavgC(end,:));
out.critical_cell=r.minimumVoltageCellIndex;
out.failure_reason=string(r.stopReason);
out.termination_time=r.stopTimeS;
out.is_mock=false;
out.solver_failed=logical(r.solverTerminatedUnexpectedly);
out.time_censored=false;
if isfield(simOpt,'optimizationCutoffActive') && ...
        simOpt.optimizationCutoffActive && ~r.success && ...
        ~r.solverTerminatedUnexpectedly && isempty(r.event.index) && ...
        r.stopTimeS>=simOpt.tMax-1e-6
    out.time_censored=true;
end
out.raw=r;
end
