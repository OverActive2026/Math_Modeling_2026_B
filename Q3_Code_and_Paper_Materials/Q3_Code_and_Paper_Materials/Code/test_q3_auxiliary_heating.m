function test_q3_auxiliary_heating
% 问题3热源的快速回归检查；正常完成时打印“全部通过”。

    %% 1. 恒功率开关与功率上限
    heating = struct('powerDensityWcm2',[0.1 0.2 0.3 0.4 0.5], ...
        'durationS',2,'maximumPowerDensityWcm2',1);
    power = q3_auxiliary_heating([0;1;2],heating);
    assert(isequal(size(power),[3 5]));
    assert(max(abs(power(1,:)-heating.powerDensityWcm2)) < 1e-14);
    assert(max(abs(power(2,:)-heating.powerDensityWcm2)) < 1e-14);
    assert(all(power(3,:) == 0));

    %% 2. 面功率-体积热源守恒与解析能耗
    opt = fast_test_options();
    heaterPower = 0.2*ones(1,5);
    heaterDuration = 0.05;
    [heatedResult,heatedModel] = pemfc_stack5_simulate_q3( ...
        'preheat',heaterPower,heaterDuration,opt);
    [~,initialOutput] = heatedModel.rhs(0,heatedModel.x0);

    expectedArealHeatWm2 = 1e4*heaterPower;
    for k = 1:5
        actualArealHeatWm2 = initialOutput.cell{k}.thermal. ...
            auxiliaryHeatFluxEquivalentWm2;
        assert(abs(actualArealHeatWm2-expectedArealHeatWm2(k)) < 1e-9);
    end

    expectedEnergyJCell = 25*heaterPower*heaterDuration;
    assert(max(abs(heatedResult.auxiliaryEnergyJCell- ...
        expectedEnergyJCell)) < 1e-12);
    assert(abs(heatedResult.totalAuxiliaryEnergyJ- ...
        sum(expectedEnergyJCell)) < 1e-12);
    assert(all(heatedResult.TavgC(end,:) > -30));

    %% 3. q_aux=0时退化回问题2方程
    zeroPower = zeros(1,5);
    q3Result = pemfc_stack5_simulate_q3('coheat',zeroPower,0,opt);
    q2Strategy = struct('type','linear','initialAcm2',0, ...
        'rampRateAcm2s',0.005,'plateauAcm2',0.3);
    q2Result = pemfc_stack5_simulate(q2Strategy,-30,opt);

    assert(isequal(q3Result.t,q2Result.t));
    assert(max(abs(q3Result.TavgC(:)-q2Result.TavgC(:))) < 1e-10);
    assert(max(abs(q3Result.Vcell(:)-q2Result.Vcell(:))) < 1e-10);
    assert(max(abs(q3Result.chargeCcm2(:)-q2Result.chargeCcm2(:))) < 1e-12);
    assert(q3Result.totalAuxiliaryEnergyJ == 0);

    fprintf('test_q3_auxiliary_heating：全部通过。\n');

end


function opt = fast_test_options()

    opt.initialTemperatureC = -30;
    opt.ambientTemperatureC = -30;
    opt.nCellLayer = [1 1 1 1 1];
    opt.tMax = 0.05;
    opt.outputStep = 0.05;
    opt.MaxStep = 0.01;
    opt.RelTol = 1e-7;
    opt.AbsTol = 1e-10;
    opt.plot = false;
    opt.verbose = false;
    opt.useJacobianPattern = true;
    opt.progressIntervalS = 0;
    opt.phaseTransitionWidthK = 0.2;

end

