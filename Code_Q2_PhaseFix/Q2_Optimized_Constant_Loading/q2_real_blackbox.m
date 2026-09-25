function [r,trajectory] = q2_real_blackbox(x,cfg)
% 接入现有五片模型。失败目标保持NaN，不将失败伪装成慢速成功。
strategy = q2_decode_parameters(x,cfg);
simOpt = cfg.simOpt;
% 恒流入口与run_q2_constant_current保持相同的时间窗：计算到20 C/cm^2
% 电荷事件之后1 s即可，无须让进度分母始终显示为600 s。
if strcmpi(cfg.strategyType,'constant') && strcmpi(strategy.type,'linear')
    j0 = strategy.initialAcm2;
    jp = strategy.plateauAcm2;
    rampRate = strategy.rampRateAcm2s;
    rampTime = cfg.constantRampTimeS;
    qLimit = simOpt.qMaxCcm2;
    qRamp = 0.5*(j0+jp)*rampTime;
    if qLimit <= qRamp && rampRate > 0
        chargeLimitTime = (-j0+sqrt(j0^2+2*rampRate*qLimit))/rampRate;
    elseif jp > 0
        chargeLimitTime = rampTime+(qLimit-qRamp)/jp;
    else
        chargeLimitTime = simOpt.tMax;
    end
    simOpt.tMax = min(simOpt.tMax,chargeLimitTime+1);
elseif strcmpi(cfg.strategyType,'step') && strcmpi(strategy.type,'step')
    levels = strategy.levelsAcm2;
    switches = strategy.switchTimesS;
    qLimit = simOpt.qMaxCcm2;
    chargeBefore = 0;
    segmentStart = 0;
    chargeLimitTime = Inf;
    for k=1:numel(switches)
        duration = switches(k)-segmentStart;
        chargeAfter = chargeBefore+levels(k)*duration;
        if qLimit<=chargeAfter && levels(k)>0
            chargeLimitTime = segmentStart+(qLimit-chargeBefore)/levels(k);
            break;
        end
        chargeBefore = chargeAfter;
        segmentStart = switches(k);
    end
    if isinf(chargeLimitTime) && levels(end)>0
        chargeLimitTime = segmentStart+(qLimit-chargeBefore)/levels(end);
    end
    if isfinite(chargeLimitTime)
        simOpt.tMax = min(simOpt.tMax,chargeLimitTime+1);
    end
end
try
    [trajectory,~] = pemfc_stack5_simulate( ...
        strategy,cfg.initialTemperatureC,simOpt);
catch simulationError
    if strcmp(simulationError.identifier, ...
            'pemfc_stack5_simulate:WallTimeLimit')
        trajectory = [];
        r = timeout_result(strategy,simOpt,simulationError.message);
        return;
    end
    rethrow(simulationError);
end
a = trajectory;
r.is_mock = false;
r.strategy = strategy;
r.success = logical(a.success) && isfinite(a.successTimeS);
r.startup_time = NaN;
r.ice_start = NaN;
if r.success
    r.startup_time = a.successTimeS;
    r.ice_start = max(a.maxIceVolumeFractionCell(end,:));
end

r.ice_end = max(a.maxIceVolumeFractionCell(end,:));
r.ice_peak = a.maximumIceVolumeFraction;
r.V_min = a.minimumCellVoltageV;
r.q_use = a.chargeUsedCcm2;
r.T_min_end = min(a.TavgC(end,:))+273.15;
r.critical_cell = a.minimumVoltageCellIndex; % 全过程最低电压对应电池
r.coldest_cell = a.coldestCellIndexAtStop;
r.failure_reason = '';
if ~r.success, r.failure_reason = a.stopReason; end
r.stop_time = a.stopTimeS;
r.solver_failed = a.solverTerminatedUnexpectedly;
r.ice_constraint_mode = 'legacy_history_peak_and_terminal_event';
c = a.constraints;
% 包含原模型的电流、孔隙保护约束；最终可行性仍以模型success为准。
margins = [(c.minimumVoltageV-r.V_min)/c.minimumVoltageV, ...
    (r.q_use-c.qMaxCcm2)/c.qMaxCcm2, ...
    (r.ice_peak-c.maximumIceVolumeFraction)/c.maximumIceVolumeFraction, ...
    (a.maximumCurrentDensityAcm2-c.jMaxAcm2)/c.jMaxAcm2, ...
    (a.maximumPoreOccupancy-c.maximumPoreOccupancy)/c.maximumPoreOccupancy, ...
    -min(a.TavgC(end,:))/max(1,abs(cfg.initialTemperatureC))];
assert(all(isfinite(margins)),'Nonfinite physical diagnostics; inspect solver.');
r.constraint_violation = max(margins);
if r.success
    r.constraint_violation = min(0,r.constraint_violation);
else
    r.constraint_violation = max(1e-6,r.constraint_violation);
end
if r.solver_failed
    % 数值失败不是物理不可行的证据，但随机优化不能因此整批中断。
    % 将其标成高违约样本并保留原因，最终候选绝不从这类点中选择。
    r.constraint_violation = max(1,r.constraint_violation);
    r.failure_reason = ['数值求解失败：' char(a.stopReason)];
end
end

function r = timeout_result(strategy,simOpt,message)
% A timed-out numerical solve is never eligible as a feasible optimum.
r = struct();
r.is_mock = false;
r.strategy = strategy;
r.success = false;
r.startup_time = NaN;
r.ice_start = NaN;
r.ice_end = NaN;
r.ice_peak = NaN;
r.V_min = NaN;
r.q_use = NaN;
r.T_min_end = NaN;
r.critical_cell = NaN;
r.coldest_cell = NaN;
r.failure_reason = ['Numerical wall-time limit: ' char(message)];
r.stop_time = NaN;
r.solver_failed = true;
r.ice_constraint_mode = 'legacy_history_peak_and_terminal_event';
r.constraint_violation = 1;
r.wall_time_limit_s = simOpt.maxWallTimeS;
end
