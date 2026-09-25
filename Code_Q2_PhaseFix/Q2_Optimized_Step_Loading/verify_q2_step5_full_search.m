% Independent, no-cutoff verification of the best nine-variable strategy.
codeDir=fileparts(mfilename('fullpath'));
addpath(codeDir);
saved=load(fullfile(codeDir,'q2_step5_full_report.mat'),'report');
report=saved.report;
assert(~isempty(report.best) && report.best.success);
quickOpt=q2_optimization_options();
quickOpt.progressIntervalS=0;
quick=q2_evaluate_step5_full_real(report.bestX,-10,quickOpt);
strictOpt=quickOpt;
strictOpt.RelTol=5e-5;
strictOpt.AbsTol=5e-9;
strict=q2_evaluate_step5_full_real(report.bestX,-10,strictOpt);
assert(quick.success && strict.success && ...
    ~quick.solver_failed && ~strict.solver_failed);
assert(abs(quick.startup_time-report.bestTimeS)<0.05);
for result=[quick,strict]
    assert(result.V_min>=quickOpt.minimumVoltageV && ...
        result.q_use<=quickOpt.qMaxCcm2 && ...
        result.ice_start<quickOpt.maximumIceVolumeFraction && ...
        result.raw.maximumPoreOccupancy<quickOpt.maximumPoreOccupancy);
end
save(fullfile(codeDir,'q2_step5_full_verification.mat'), ...
    'quick','strict','quickOpt','strictOpt','-v7.3');
fprintf('Nine-variable verified quick %.9f s, strict %.9f s\n', ...
    quick.startup_time,strict.startup_time);
fprintf('Vmin %.6f / %.6f V; charge %.6f / %.6f C/cm^2\n', ...
    quick.V_min,strict.V_min,quick.q_use,strict.q_use);
fprintf('Ice at start %.6f / %.6f; pore occupancy %.6f / %.6f\n', ...
    quick.ice_start,strict.ice_start, ...
    quick.raw.maximumPoreOccupancy,strict.raw.maximumPoreOccupancy);
