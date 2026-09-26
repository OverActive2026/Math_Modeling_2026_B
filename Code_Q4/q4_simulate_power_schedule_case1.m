function [result,model] = q4_simulate_power_schedule_case1(schedule,simOpt,model)
%Q4_SIMULATE_POWER_SCHEDULE_CASE1 Simulate five-heater open/closed-loop power.
% schedule.breaksS: strictly increasing times starting at zero [s].
% schedule.powerWcm2: one 1x5 power row per interval [W/cm^2].
% Heating ends immediately at the first successful startup event. Passing
% a previously constructed model avoids rebuilding it during optimization.
% With no initialState/File, the default remains uniform -30 C (case 1).
% A supplied precooling state enables Q4 cases 2 and 3 without regridding.

    if nargin < 2 || isempty(simOpt), simOpt = struct(); end
    if nargin < 3, model = []; end
    if ~isstruct(schedule) || ~isscalar(schedule) || ...
            ~isfield(schedule,'breaksS') || ~isfield(schedule,'powerWcm2')
        error('q4_schedule:InvalidSchedule','Provide breaksS and powerWcm2.');
    end
    breaks = schedule.breaksS(:);
    powerRows = schedule.powerWcm2;
    if numel(breaks) < 2 || abs(breaks(1)) > 1e-10 || ...
            any(~isfinite(breaks)) || any(diff(breaks) <= 0) || ...
            ~isequal(size(powerRows),[numel(breaks)-1,5]) || ...
            any(~isfinite(powerRows(:))) || ...
            any(powerRows(:) < 0) || any(powerRows(:) > 1)
        error('q4_schedule:InvalidSchedule', ...
            'Breaks must increase from zero and power must be Nx5 in [0,1].');
    end
    feedbackEnabled = isfield(schedule,'feedbackFcn') && ...
        ~isempty(schedule.feedbackFcn);
    if feedbackEnabled
        if ~isa(schedule.feedbackFcn,'function_handle') || ...
                ~isfield(schedule,'sampleS') || ...
                ~isscalar(schedule.sampleS) || ...
                ~isfinite(schedule.sampleS) || schedule.sampleS <= 0
            error('q4_schedule:InvalidFeedback', ...
                'Feedback requires a function handle and positive sampleS.');
        end
        originalBreaks = breaks;
        originalPower = powerRows;
        breaks = unique([breaks; ...
            (0:schedule.sampleS:breaks(end)).';breaks(end)]);
        powerRows = zeros(numel(breaks)-1,5);
        for j = 1:size(powerRows,1)
            originalIndex = find(originalBreaks <= breaks(j)+1e-9,1,'last');
            powerRows(j,:) = originalPower(originalIndex,:);
        end
    end
    if ~isfield(simOpt,'successMarginC'), simOpt.successMarginC = 0.01; end
    if ~isfield(simOpt,'recordStepS'), simOpt.recordStepS = 0; end
    if ~isfield(simOpt,'progressEveryS'), simOpt.progressEveryS = 0; end
    if ~isscalar(simOpt.successMarginC) || simOpt.successMarginC <= 0 || ...
            ~isfinite(simOpt.successMarginC) || ...
            ~isscalar(simOpt.recordStepS) || simOpt.recordStepS < 0 || ...
            ~isfinite(simOpt.recordStepS) || ...
            ~isscalar(simOpt.progressEveryS) || ...
            simOpt.progressEveryS < 0 || ~isfinite(simOpt.progressEveryS)
        error('q4_schedule:InvalidOptions','Invalid numerical options.');
    end

    if isempty(model)
        % Preserve the case-1 defaults, but honor an explicitly supplied
        % spatial temperature field from 20/40-minute precooling.
        if ~isfield(simOpt,'initialTemperatureC')
            simOpt.initialTemperatureC = -30;
        end
        if ~isfield(simOpt,'ambientTemperatureC')
            simOpt.ambientTemperatureC = -30;
        end
        if ~isfield(simOpt,'initialState'), simOpt.initialState = []; end
        if ~isfield(simOpt,'initialStateFile')
            simOpt.initialStateFile = '';
        end
        if ~isfield(simOpt,'nCellLayer'), simOpt.nCellLayer = [2 3 4 4 2]; end
        if ~isfield(simOpt,'kMelt'), simOpt.kMelt = 1; end
        if ~isfield(simOpt,'MaxStep'), simOpt.MaxStep = 0.05; end
        if ~isfield(simOpt,'RelTol'), simOpt.RelTol = 1e-4; end
        if ~isfield(simOpt,'AbsTol'), simOpt.AbsTol = 1e-8; end
        simOpt.tMax = breaks(end);
        simOpt.verbose = false;
        simOpt.plot = false;
        setupOpt = simOpt;
        setupOpt.tMax = 1e-3;
        setupOpt.outputStep = 1e-3;
        setupOpt.MaxStep = min(simOpt.MaxStep,1e-3);
        [~,model] = pemfc_stack5_simulate_q4( ...
            'coheat',zeros(1,5),0,setupOpt);
        model.options.MaxStep = simOpt.MaxStep;
        model.options.tMax = simOpt.tMax;
    end
    opt = model.options;
    odeBase = odeset('RelTol',opt.RelTol,'AbsTol',opt.AbsTol, ...
        'MaxStep',opt.MaxStep,'NonNegative',model.nonNegativeStates);
    if opt.useJacobianPattern
        odeBase = odeset(odeBase,'JPattern',model.jacobianPattern);
    end
    nextProgressS = simOpt.progressEveryS;
    if simOpt.progressEveryS > 0
        odeBase = odeset(odeBase,'OutputFcn',@progress_output);
    end

    x = model.x0;
    t = 0;
    [~,out] = model.rhsForPower(t,x,zeros(1,5));
    timeS = t;
    temperatureC = out.Tavg.'-273.15;
    voltageV = out.Vcell.';
    iceFraction = out.maxIceVolumeFractionCell.';
    poreOccupancy = out.maxPoreOccupancyCell.';
    powerUsed = zeros(0,5);
    nominalPowerUsed = zeros(0,5);
    intervalStartS = zeros(0,1);
    intervalEndS = zeros(0,1);
    energyJCell = zeros(1,5);
    success = false;
    eventIndex = 0;
    reason = 'Schedule ended before successful startup';
    if min(voltageV) < opt.minimumVoltageV || ...
            max(iceFraction) >= opt.maximumIceVolumeFraction || ...
            max(poreOccupancy) >= opt.maximumPoreOccupancy
        error('q4_schedule:InvalidInitialState','Initial safety limit violated.');
    end

    clock = tic;
    for j = 1:numel(breaks)-1
        if abs(t-breaks(j)) > 1e-7
            error('q4_schedule:TimeMismatch','Schedule intervals are not contiguous.');
        end
        nominalPower = powerRows(j,:);
        power = nominalPower;
        if feedbackEnabled
            if numel(timeS) == 1
                dTdt = zeros(1,5);
                dVdt = zeros(1,5);
            else
                deltaT = timeS(end)-timeS(end-1);
                dTdt = (temperatureC(end,:)-temperatureC(end-1,:))/deltaT;
                dVdt = (voltageV(end,:)-voltageV(end-1,:))/deltaT;
            end
            power = schedule.feedbackFcn(t,temperatureC(end,:), ...
                voltageV(end,:),dTdt,dVdt,nominalPower);
            if ~isnumeric(power) || ~isequal(size(power),[1,5]) || ...
                    any(~isfinite(power)) || any(power < 0) || any(power > 1)
                error('q4_schedule:InvalidFeedbackPower', ...
                    'Feedback must return a finite 1x5 power row in [0,1].');
            end
        end
        rhs = @(tt,xx) model.rhsForPower(tt,xx,power);
        events = @(tt,xx) local_events(tt,xx,rhs,opt, ...
            simOpt.successMarginC);
        segmentOpt = odeset(odeBase,'Events',events);
        if simOpt.recordStepS > 0
            tspan = unique([t:simOpt.recordStepS:breaks(j+1),breaks(j+1)]);
        else
            tspan = [t,breaks(j+1)];
        end
        [tSeg,xSeg,~,~,iEvent] = ode15s(rhs,tspan,x,segmentOpt);
        tNew = tSeg(end);
        x = xSeg(end,:).';
        intervalStartS(end+1,1) = t; %#ok<AGROW>
        intervalEndS(end+1,1) = tNew; %#ok<AGROW>
        powerUsed(end+1,:) = power; %#ok<AGROW>
        nominalPowerUsed(end+1,:) = nominalPower; %#ok<AGROW>
        energyJCell = energyJCell + opt.cellAreaM2*1e4*power*(tNew-t);
        for m = 2:numel(tSeg)
            [~,out] = rhs(tSeg(m),xSeg(m,:).');
            timeS(end+1,1) = tSeg(m); %#ok<AGROW>
            temperatureC(end+1,:) = out.Tavg.'-273.15; %#ok<AGROW>
            voltageV(end+1,:) = out.Vcell.'; %#ok<AGROW>
            iceFraction(end+1,:) = out.maxIceVolumeFractionCell.'; %#ok<AGROW>
            poreOccupancy(end+1,:) = out.maxPoreOccupancyCell.'; %#ok<AGROW>
        end
        t = tNew;
        if ~isempty(iEvent)
            eventIndex = iEvent(end);
            switch eventIndex
                case 1
                    success = true;
                    reason = 'All five cell-average temperatures exceeded 0 C';
                case 2
                    reason = 'Cell voltage reached lower limit';
                case 3
                    reason = 'Local ice fraction reached upper limit';
                case 4
                    reason = 'Current density reached upper limit';
                case 5
                    reason = 'Pore occupancy reached numerical safety limit';
            end
            break;
        end
        if tNew < breaks(j+1)-1e-8
            reason = 'ODE solver stopped early';
            break;
        end
    end
    result.success = success;
    result.stopReason = reason;
    result.eventIndex = eventIndex;
    result.stopTimeS = t;
    result.startupTimeS = NaN;
    if success, result.startupTimeS = t; end
    result.energyJCell = energyJCell;
    result.totalAuxiliaryEnergyJ = sum(energyJCell);
    result.timeS = timeS;
    result.temperatureC = temperatureC;
    result.voltageV = voltageV;
    result.iceVolumeFraction = iceFraction;
    result.poreOccupancy = poreOccupancy;
    result.intervalStartS = intervalStartS;
    result.intervalEndS = intervalEndS;
    result.powerDensityWcm2 = powerUsed;
    result.nominalPowerWcm2 = nominalPowerUsed;
    result.minimumCellVoltageV = min(voltageV(:));
    result.maximumIceVolumeFraction = max(iceFraction(:));
    result.maximumPoreOccupancy = max(poreOccupancy(:));
    result.maximumCellTemperatureSpreadC = max( ...
        max(temperatureC,[],2)-min(temperatureC,[],2));
    result.finalTemperatureC = temperatureC(end,:);
    result.finalState = x;
    result.options = opt;
    result.cpuTimeS = toc(clock);

    function stop = progress_output(currentTime,~,flag)
        stop = false;
        if isempty(flag) && ~isempty(currentTime) && ...
                currentTime(end) >= nextProgressS-1e-9
            fprintf('仿真进度：已计算到 t=%.2f s\n', ...
                currentTime(end));
            nextProgressS = nextProgressS + simOpt.progressEveryS* ...
                (1+floor((currentTime(end)-nextProgressS)/ ...
                simOpt.progressEveryS));
        end
    end
end

function [value,isTerminal,direction] = local_events(t,x,rhs,opt,marginC)
    [~,out] = rhs(t,x);
    value = [min(out.Tavg)-273.15-marginC; ...
        min(out.Vcell)-opt.minimumVoltageV; ...
        opt.maximumIceVolumeFraction-max(out.maxIceVolumeFractionCell); ...
        opt.jMaxAcm2+1e-8-out.jAcm2; ...
        opt.maximumPoreOccupancy-max(out.maxPoreOccupancyCell)];
    isTerminal = ones(5,1);
    direction = [1;-1;-1;-1;-1];
end
