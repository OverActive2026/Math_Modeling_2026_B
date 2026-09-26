function cfg = q3_optimization_config(profile)
% Fixed physical assumptions and separate numerical search settings.
    if nargin < 1, profile = 'smoke'; end
    cfg.profile = char(profile);
    cfg.timeUpperS = 600;
    cfg.coheatObservationS = 600;
    cfg.acceptanceMode = 'question';
    cfg.seed = 20260925;
    cfg.guidance = false;
    cfg.simOpt = struct('initialTemperatureC',-30, ...
        'ambientTemperatureC',-30,'gammaIce',3.5,'kFreeze',0.4, ...
        'plot',false,'verbose',false,'progressIntervalS',0, ...
        'useJacobianPattern',true,'enforceChargeCap',true, ...
        'iceScope','terminal', ...
        'temperatureEventMarginK',0.01);
    switch lower(cfg.profile)
        case 'smoke'
            cfg.timeUpperS = 180;
            cfg.coheatObservationS = 180;
            cfg.maxEvaluations = 12;
            cfg.simOpt.nCellLayer = [2 3 4 4 2];
            cfg.simOpt.outputStep = 0.5;
            cfg.simOpt.RelTol = 1e-4;
            cfg.simOpt.AbsTol = 1e-8;
            cfg.simOpt.MaxStep = 0.05;
        case 'full'
            cfg.maxEvaluations = 2000;
            % Same grid as the base model, per the current comparison scope.
            cfg.simOpt.nCellLayer = [2 3 4 4 2];
            cfg.simOpt.outputStep = 0.5;
            cfg.simOpt.RelTol = 1e-5;
            cfg.simOpt.AbsTol = 1e-9;
            cfg.simOpt.MaxStep = 0.05;
        otherwise
            error('q3_optimization_config:Profile', ...
                'profile must be smoke or full.');
    end
end
