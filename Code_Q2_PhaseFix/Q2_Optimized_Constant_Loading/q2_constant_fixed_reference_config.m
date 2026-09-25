function cfg = q2_constant_fixed_reference_config()
% Freeze the successful constant loading and build a fresh search experiment.
codeDir = fileparts(mfilename('fullpath'));
saved = load(fullfile(codeDir,'optimization_results', ...
    'q2_constant_fixed_reference_gamma3p5_kfreeze03.mat'),'result');
reference = saved.result;
assert(reference.success && strcmp(reference.integrationScheme, ...
    'original_ode15s_single_call_v1'));
assert(reference.options.gammaIce==3.5 && reference.options.kFreeze==0.3 && ...
    reference.options.jMaxAcm2==0.5 && ...
    isequal(reference.options.nCellLayer,[2 3 4 4 2]));
assert(abs(reference.strategy.initialAcm2-0.02)<1e-12 && ...
    abs(reference.strategy.plateauAcm2-0.259937538852887)<1e-12);
cfg = q2_optimization_config();
cfg.strategyType = 'constant';
cfg.experiment = 'constant_fixed_0259937538852887_gamma3p5_kfreeze03_v1';
cfg.seed = 20260927;
cfg.budget = 30;
cfg.nInit = 12;
cfg.currentMax = 0.5;
cfg.constantInitialCurrentAcm2 = 0.02;
cfg.constantRampTimeS = 2;
cfg.constantPlateauBounds = [0.02 0.5];
cfg.initialNormalizedPoint = ...
    (reference.strategy.plateauAcm2-cfg.constantPlateauBounds(1))/ ...
    diff(cfg.constantPlateauBounds);
% LHS initial points focus near the verified feasible reference; subsequent
% FuRBO/LogCEI search still allows the full plateau interval [0.02,0.5].
cfg.initialDesignBounds = ([0.245;0.270]-0.02)/0.48;
assert(cfg.initialNormalizedPoint>=cfg.initialDesignBounds(1) && ...
    cfg.initialNormalizedPoint<=cfg.initialDesignBounds(2));
cfg.simOpt = reference.options;
cfg.simOpt.tMax = 600;
cfg.simOpt.maxWallTimeS = 300;
cfg.initialTemperatureC = reference.initialTemperatureC;
cfg.reference = reference;
cfg.referenceTimeS = reference.successTimeS;
cfg.outputDir = fullfile(codeDir,'optimization_results', ...
    'constant_fixed_reference_kfreeze03');
decoded = q2_decode_parameters(cfg.initialNormalizedPoint,cfg);
assert(abs(decoded.plateauAcm2-reference.strategy.plateauAcm2)<1e-12);
end
