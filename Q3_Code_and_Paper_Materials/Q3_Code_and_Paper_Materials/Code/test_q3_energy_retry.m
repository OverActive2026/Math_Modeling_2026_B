function test_q3_energy_retry()
% A numerical failure must trigger a genuinely different integration attempt.
    seen={};
    cfg=q3_optimization_config('smoke'); opt=cfg.simOpt;
    opt.RelTol=1e-6; opt.AbsTol=1e-10; opt.MaxStep=.0125;
    c=q3_evaluate_energy6('coheat',ones(1,5),30,opt,@fake_simulator);
    assert(strcmp(c.status,'success') && numel(seen)==2);
    assert(seen{2}.RelTol<seen{1}.RelTol && ...
        seen{2}.AbsTol<seen{1}.AbsTol && seen{2}.MaxStep<seen{1}.MaxStep, ...
        'Numerical retry repeated the failed integration settings.');
    assert(seen{1}.tMax>96.6666 && seen{2}.tMax==seen{1}.tMax);
    assert(isequal(seen{1}.nCellLayer,seen{2}.nCellLayer));
    disp('test_q3_energy_retry passed');
    function r=fake_simulator(~,~,~,o)
        seen{end+1}=o;
        r=struct('totalAuxiliaryEnergyJ',100,'successTimeS',31, ...
            'minimumCellVoltageV',.6,'maximumIceVolumeFraction',.1, ...
            'chargeUsedCcm2',2.4,'stopReason','mock', ...
            'solverTerminatedUnexpectedly',numel(seen)==1, ...
            'success',numel(seen)==2);
    end
end
