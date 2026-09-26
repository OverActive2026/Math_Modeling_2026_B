function out=q2_evaluate_schedule(st,temps,opt,timeCapS,oracle)
% Joint feasibility, warm-only objective; budgets are not physical failures.
if nargin<5, oracle=@pemfc_stack5_simulate; end
assert(isrow(temps) && temps(1)==-10 && all(diff(temps)<0));
assert(isfinite(timeCapS) && timeCapS>0);
out=struct('strategy',st,'success',false,'startup_time',NaN, ...
    'warmStartupTimeS',NaN,'V_min',0,'q_use',0,'ice_peak',0, ...
    'ice_start',NaN,'T_min_end',263.15,'termination_time',NaN, ...
    'critical_cell',NaN,'is_mock',false,'solver_failed',false, ...
    'time_censored',false,'robustConstraint',1,'raw',struct(), ...
    'temperaturesC',temps,'status',zeros(size(temps)), ...
    'scenarios',{cell(size(temps))},'failure_reason',"");
% 1 feasible, -1 physical failure, -2 numerical unknown,
% -3 budget/objective censored, 0 not evaluated.
allMargins=[];
for k=1:numel(temps)
    thisOpt=opt;
    if k>=2
        assert(isfield(opt,'coldTMaxS') && opt.coldTMaxS>0, ...
            'Q2:ColdHorizon','Explicit cold computational horizon required.');
        thisOpt.tMax=opt.coldTMaxS;
        if isfield(thisOpt,'optimizationCutoffActive')
            thisOpt=rmfield(thisOpt,'optimizationCutoffActive');
        end
    end
    fprintf('T0=%g C; simulation horizon=%g s\n',temps(k),thisOpt.tMax);
    try
        [r,~]=oracle(st,temps(k),thisOpt);
        if r.solverTerminatedUnexpectedly
            error('Q2:IncompleteIntegration','Solver terminated unexpectedly');
        end
    catch ME
        out.status(k)=-2; out.solver_failed=true;
        out.failure_reason="NUMERICAL_UNKNOWN: "+string(ME.identifier); return
    end
    out.scenarios{k}=r;
    if k==1
        out.raw=r; out.V_min=r.minimumCellVoltageV;
        out.q_use=r.chargeUsedCcm2; out.ice_peak=r.maximumIceVolumeFraction;
        out.T_min_end=273.15+min(r.TavgC(end,:));
        out.termination_time=r.stopTimeS;
        out.critical_cell=r.minimumVoltageCellIndex;
        out.warmStartupTimeS=r.successTimeS;
    end
    fprintf('T0=%g success=%d startup=%g s Vmin=%.8f V\n', ...
        temps(k),r.success,r.successTimeS,r.minimumCellVoltageV);
    if ~r.success && isempty(r.event.index)
        out.status(k)=-3; out.time_censored=true;
        out.failure_reason="BUDGET_CENSORED at T="+temps(k); return
    end
    margins=[(opt.minimumVoltageV-r.minimumCellVoltageV)/opt.minimumVoltageV, ...
        (r.chargeUsedCcm2-opt.qMaxCcm2)/opt.qMaxCcm2, ...
        (r.maximumCurrentDensityAcm2-opt.jMaxAcm2)/opt.jMaxAcm2, ...
        (r.maximumIceVolumeFraction-opt.maximumIceVolumeFraction)/opt.maximumIceVolumeFraction, ...
        (r.maximumPoreOccupancy-opt.maximumPoreOccupancy)/opt.maximumPoreOccupancy];
    assert(all(isfinite(margins)));
    if ~r.success || ~isfinite(r.successTimeS) || max(margins)>1e-6
        out.status(k)=-1; out.robustConstraint=max(1e-5,max(margins));
        out.failure_reason="PHYSICAL_FAILURE at T="+temps(k)+": "+string(r.stopReason); return
    end
    if k==1 && r.successTimeS>timeCapS
        out.status(k)=-3; out.time_censored=true;
        out.failure_reason="WARM_TIME_OBJECTIVE_REJECTED"; return
    end
    out.status(k)=1; allMargins=[allMargins margins]; %#ok<AGROW>
    if k==1, out.ice_start=max(r.maxIceVolumeFractionCell(end,:)); end
end
out.success=all(out.status==1);
if out.success
    out.startup_time=out.warmStartupTimeS;
    out.robustConstraint=min(-1e-5,max(allMargins));
    out.failure_reason="All required temperatures physically feasible";
end
end
