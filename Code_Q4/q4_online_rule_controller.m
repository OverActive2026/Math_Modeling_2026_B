function [power,details] = q4_online_rule_controller( ...
    temperatureC,voltageV,dTdt,dVdt,p)
%Q4_ONLINE_RULE_CONTROLLER State-driven five-cell Q4 heating controller.
% Uses only the latest cell-average T, V and their measured slopes. There
% are no preset switch-off times. The thermal thresholds were calibrated
% from a feasible case-1 trajectory and are NOT globally optimal values.

    if nargin < 5 || isempty(p), p = struct(); end
    p = default_value(p,'centerEndOffC',-1.20);
    p = default_value(p,'centerTaperBandC',0.40);
    p = default_value(p,'centerWarmReserveC',14.0);
    p = default_value(p,'adjacentEndOffC',-0.18);
    p = default_value(p,'adjacentTaperBandC',0.30);
    p = default_value(p,'adjacentWarmReserveC',11.0);
    p = default_value(p,'minimumEndRiseCps',0.35);
    p = default_value(p,'voltageWarningV',0.45);
    p = default_value(p,'voltageLookaheadS',1.0);
    p = default_value(p,'iceRiskVoltageV',0.55);
    p = default_value(p,'iceRiskDropVps',0.02);

    vectors = {temperatureC,voltageV,dTdt,dVdt};
    if any(cellfun(@(v) ~isnumeric(v) || ~isequal(size(v),[1,5]) || ...
            any(~isfinite(v)),vectors)) || ...
            p.centerTaperBandC <= 0 || p.adjacentTaperBandC <= 0 || ...
            p.voltageLookaheadS <= 0
        error('q4_online:InvalidInput', ...
            'Measurements must be finite 1x5 rows and bands positive.');
    end

    power = ones(1,5);
    endTemperatureC = min(temperatureC([1,5]));
    endRiseCps = min(dTdt([1,5]));
    details.temperatureLow = false(1,5);
    details.temperatureLow(3) = ...
        temperatureC(3) < p.centerWarmReserveC;
    details.temperatureLow([2,4]) = ...
        temperatureC([2,4]) < p.adjacentWarmReserveC;
    details.insufficientRise = endTemperatureC < 0 && ...
        endRiseCps < p.minimumEndRiseCps;
    details.nearSuccess = endTemperatureC >= p.adjacentEndOffC;

    % The warm center can coast first, but only when the cold ends have
    % nearly recovered and the center has enough stored thermal energy.
    if ~details.temperatureLow(3)
        power(3) = clamp01( ...
            (p.centerEndOffC-endTemperatureC)/p.centerTaperBandC);
    end
    % Cells 2/4 are still valuable heat sources for the cold end cells.
    % Their taper begins much later, close to stack startup.
    for k = [2,4]
        if ~details.temperatureLow(k)
            power(k) = clamp01( ...
                (p.adjacentEndOffC-endTemperatureC)/ ...
                p.adjacentTaperBandC);
        end
    end
    if details.insufficientRise
        power(2:4) = 1;
    end

    predictedVoltageV = voltageV + ...
        p.voltageLookaheadS*min(dVdt,0);
    details.voltageRisk = predictedVoltageV < p.voltageWarningV;
    details.iceRisk = temperatureC < 0 & ...
        predictedVoltageV < p.iceRiskVoltageV & ...
        dVdt < -p.iceRiskDropVps;
    risk = details.voltageRisk | details.iceRisk;
    for k = find(risk)
        % Heat both the threatened cell and its nearest thermal neighbors.
        power(max(1,k-1):min(5,k+1)) = 1;
    end
    details.reducedPower = power < 1-1e-10;
    details.stoppedPower = power <= 1e-10;
    details.endTemperatureC = endTemperatureC;
    details.endRiseCps = endRiseCps;
    details.predictedVoltageV = predictedVoltageV;
end

function x = clamp01(x)
    x = min(1,max(0,x));
end

function s = default_value(s,name,value)
    if ~isfield(s,name) || isempty(s.(name)), s.(name) = value; end
end
