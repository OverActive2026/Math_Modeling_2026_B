function test_q3_heater_switch_integration()
% A heater switch must not perturb the preceding constant-heating segment.
    opt = struct('initialTemperatureC',-30,'ambientTemperatureC',-30, ...
        'nCellLayer',[1 1 1 1 1],'tMax',0.2,'outputStep',0.05, ...
        'MaxStep',0.01,'RelTol',1e-5,'AbsTol',1e-9, ...
        'plot',false,'verbose',false,'useJacobianPattern',true, ...
        'gammaIce',3.5,'kFreeze',0.4);
    q = [1 0.7 0.4 0.8 0.9];
    th = 0.1;
    switched = pemfc_stack5_simulate_q3('coheat',q,th,opt);
    prefixOpt = opt;
    prefixOpt.tMax = th;
    prefix = pemfc_stack5_simulate_q3('coheat',q,1,prefixOpt);
    k = find(abs(switched.t-th) < 1e-12);
    assert(isscalar(k),'The switch time must occur exactly once.');
    scaledError = max(abs(switched.x(k,:)-prefix.x(end,:)) ./ ...
        max(1,abs(prefix.x(end,:))));
    assert(scaledError < 1e-12, ...
        'The future heater switch perturbed the continuous-heating prefix.');
    assert(all(diff(switched.t) > 0));
    assert(abs(switched.stopTimeS-0.2) < 1e-12);
    assert(abs(switched.totalAuxiliaryEnergyJ-9.5) < 1e-10);
    assert(max(abs(switched.auxiliaryEnergyJCell-[2.5 1.75 1 2 2.25])) < 1e-10);
    assert(abs(switched.chargeCcm2(end)-0.0001) < 1e-8, ...
        'Restarting the solver must preserve accumulated charge.');

    failureOpt = opt;
    failureOpt.qMaxCcm2 = 1e-6;
    failure = pemfc_stack5_simulate_q3('coheat',q,th,failureOpt);
    assert(~failure.success && any(failure.event.index == 2));
    assert(failure.stopTimeS < th && ...
        abs(failure.chargeCcm2(end)-1e-6) < 1e-9, ...
        'A terminal failure must prevent the heater-off segment from starting.');
    assert(abs(failure.totalAuxiliaryEnergyJ-95*failure.stopTimeS) < 1e-10);

    successOpt = opt;
    successOpt.initialTemperatureC = 0.009;
    successOpt.ambientTemperatureC = 0.009;
    successful = pemfc_stack5_simulate_q3('coheat',ones(1,5),th,successOpt);
    assert(successful.success && successful.stopTimeS < th, ...
        'A terminal success must prevent the heater-off segment from starting.');
    assert(sum(successful.event.index == 1) == 1);
    assert(abs(successful.totalAuxiliaryEnergyJ-125*successful.stopTimeS) < 1e-10);
    fprintf('test_q3_heater_switch_integration: passed\n');
end
