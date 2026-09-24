function result = simulate_strategy(strategyType,strategyPar,p)
%SIMULATE_STRATEGY 同一电流下推进五片电池并逐时检查约束。
% 时间 [s]，电流密度 [A/cm^2]，累计电荷 [C/cm^2]。
if ~isstruct(strategyPar) || ~isscalar(strategyPar) || ...
        ~isfinite(p.dt) || p.dt <= 0 || ...
        ~isfinite(p.t_max) || p.t_max <= 0
    error('Q2:InvalidSimulationInput','策略参数、dt 或 t_max 无效。');
end
strategyType = lower(string(strategyType));
if ~isscalar(strategyType) || ...
        ~ismember(strategyType,["constant","linear","step"])
    error('Q2:UnknownStrategy','策略只能为 constant、linear 或 step。');
end
strategyPar.j_max = p.j_max;
state = initial_state(p);
q_use = 0;
time = (0:p.dt:p.t_max).';
Nstep = numel(time);
T_history = NaN(Nstep,p.Ncell);
V_history = NaN(Nstep,p.Ncell);
ice_history = NaN(Nstep,p.Ncell);
lambda_history = NaN(Nstep,p.Ncell);
current_history = NaN(Nstep,1);
T_history(1,:) = state.T.';
V_history(1,:) = state.V.';
ice_history(1,:) = state.ice_max.';
lambda_history(1,:) = state.lambda_mem.';

success = false;
failure = false;
reason = "none";
criticalCell = NaN;
startupTime = NaN;
lastIndex = Nstep;
for n = 1:Nstep-1
    t = time(n);
    switch strategyType
        case "constant"
            j = current_constant(t,strategyPar);
        case "linear"
            j = current_linear(t,strategyPar);
        case "step"
            j = current_step(t,strategyPar);
    end
    if ~isfinite(j) || j < 0 || j > p.j_max
        error('Q2:InvalidCurrent','策略生成了超限电流密度。');
    end
    current_history(n) = j;
    q_use = q_use+j*p.dt; % [A/cm^2]*[s]=[C/cm^2]
    state = stack_step(state,j,p);
    T_history(n+1,:) = state.T.';
    V_history(n+1,:) = state.V.';
    ice_history(n+1,:) = state.ice_max.';
    lambda_history(n+1,:) = state.lambda_mem.';
    current_history(n+1) = j;
    [success,failure,reason,criticalCell] = ...
        check_startup_state(state,q_use,p);
    if success || failure
        lastIndex = n+1;
        if success, startupTime = time(lastIndex); end
        break
    end
end
if ~success && ~failure
    failure = true;
    reason = "time_limit";
end

result.success = success;
result.failure = failure;
result.startup_time = startupTime;
result.q_use = q_use;
result.j_max = max(current_history(1:lastIndex),[],'omitnan');
validV = V_history(1:lastIndex,:);
validV = validV(isfinite(validV));
if isempty(validV), result.V_min = NaN; else, result.V_min = min(validV); end
result.ice_max = max(ice_history(1:lastIndex,:),[],'all');
result.critical_cell = criticalCell;
result.failure_reason = reason;
result.time = time(1:lastIndex);
result.current = current_history(1:lastIndex);
result.T_history = T_history(1:lastIndex,:);
result.V_history = V_history(1:lastIndex,:);
result.ice_history = ice_history(1:lastIndex,:);
result.lambda_history = lambda_history(1:lastIndex,:);
result.final_state = state;
result.strategy_type = strategyType;
result.is_test_placeholder = isfield(p,'test') && p.test.enabled;
end
