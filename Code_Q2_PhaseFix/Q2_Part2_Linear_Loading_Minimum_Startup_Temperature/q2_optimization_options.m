function simOpt = q2_optimization_options()
%Q2_OPTIMIZATION_OPTIONS Shared model settings for all three strategies.
% Uses the quick-test settings in the three original run_q2_* scripts.
simOpt.tMax = 600;
simOpt.nCellLayer = [2 3 4 4 2];
simOpt.RelTol = 1e-4;
simOpt.AbsTol = 1e-8;
simOpt.MaxStep = 0.05;
simOpt.outputStep = 0.5;
simOpt.plot = false;
simOpt.verbose = false;
simOpt.progressIntervalS = 5; % Print simulated time every 5 s.
simOpt.useJacobianPattern = true;
simOpt.endPlateModel = 'separate';
simOpt.bipolarPlateCounting = 'shared_stack';
simOpt.init.lambdaCCL0 = 2.0;
simOpt.j0Ref = 0.527582014454;
simOpt.Ea = 8000;
simOpt.hVaporCathode = 0.01;
simOpt.hVaporAnode = 0.01;
simOpt.interfaceRateFactor = 0.01;
simOpt.kFreeze = 0.4; % 1/s; midpoint of the requested 0.3-0.5 range.
simOpt.gammaIce = 3.5; % Dimensionless CCL ice-coverage exponent.
simOpt.kMelt = 1;
simOpt.phaseTransitionWidthK = 0.2;
simOpt.directFreezeFraction = 0.01;
simOpt.freezeNucleationLambdaFraction = 0.5;
simOpt.qMaxCcm2 = 20;
simOpt.jMaxAcm2 = 0.5;
simOpt.minimumVoltageV = 0.30;
simOpt.maximumIceVolumeFraction = 0.99;
% Unify the inconsistent 0.98/0.99 settings in the source scripts.
simOpt.maximumPoreOccupancy = 0.98;
end
