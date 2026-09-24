function [success,failure,reason,criticalCell] = check_startup_state(state,q_use,p)
%CHECK_STARTUP_STATE 逐时检查电压、冰、电荷及五片同时启动。
% q_use [C/cm^2]；T [K]；V [V]；ice_max [-]。
success = false;
failure = false;
reason = "none";
criticalCell = NaN;
if ~isscalar(q_use) || ~isfinite(q_use) || q_use < 0
    error('Q2:InvalidCharge','累计电荷量必须为非负有限数 [C/cm^2]。');
end
if any(~isfinite(state.T)) || any(~isfinite(state.V)) || ...
        any(~isfinite(state.ice_max))
    failure = true;
    reason = "model_incomplete";
    return
end
[vMin,kV] = min(state.V);
[iceMax,kIce] = max(state.ice_max);
if vMin < p.V_min_limit
    failure = true;
    reason = "voltage_limit";
    criticalCell = kV;
elseif iceMax >= p.ice_limit
    failure = true;
    reason = "ice_blockage";
    criticalCell = kIce;
elseif q_use > p.q_max
    failure = true;
    reason = "charge_limit";
elseif all(state.T > p.T_success)
    success = true;
end
end
