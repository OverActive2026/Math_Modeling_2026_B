%RUN_Q2_STEP_LOADING_COMPARISON Re-simulate baseline and optimized steps.
% Runs both curves on the same physical model and quick grid, then saves
% a figure and a MAT result for the team's report. No mock data are used.
codeDir=fileparts(mfilename('fullpath'));
addpath(codeDir);
simOpt=q2_optimization_options();
simOpt.progressIntervalS=0;
baselineStrategy=struct('type','step', ...
    'levelsAcm2',[0.15 0.18 0.30 0.40], ...
    'switchTimesS',[5 15 35]);
published=load(fullfile(codeDir,'q2_step5_full_report.mat'),'report');
optimizedStrategy=published.report.bestStrategy;
optimizedX=q2_step5_full_encode(optimizedStrategy);
[baseline,~]=pemfc_stack5_simulate(baselineStrategy,-10,simOpt);
optimized=q2_evaluate_step5_full_real(optimizedX,-10,simOpt);
assert(baseline.success && optimized.success, ...
    'At least one strategy failed; do not publish this comparison.');
assert(abs(optimized.startup_time-published.report.bestTimeS)<0.05, ...
    'Recomputed optimized result differs from the published quick-grid run.');
assert(optimized.V_min>=simOpt.minimumVoltageV && ...
    optimized.q_use<=simOpt.qMaxCcm2 && ...
    optimized.ice_start<simOpt.maximumIceVolumeFraction && ...
    optimized.raw.maximumPoreOccupancy<simOpt.maximumPoreOccupancy);
r=optimized.raw;
fprintf('Unoptimized four-level startup: %.9f s\n',baseline.successTimeS);
fprintf('Optimized five-level startup: %.9f s\n',optimized.startup_time);
fprintf('Reduction: %.9f s (%.2f%%)\n', ...
    baseline.successTimeS-optimized.startup_time, ...
    100*(baseline.successTimeS-optimized.startup_time)/baseline.successTimeS);
fprintf('Optimized Vmin %.6f V, charge %.6f C/cm^2\n', ...
    optimized.V_min,optimized.q_use);
fig=figure('Name','Q2 step-loading comparison','Color','w', ...
    'Position',[100 100 1100 720]);
layout=tiledlayout(fig,2,2,'TileSpacing','compact');
nexttile(layout);
stairs(baseline.t,baseline.currentDensityAcm2,'k--','LineWidth',1.4);
hold on;
stairs(r.t,r.currentDensityAcm2,'b','LineWidth',1.4);
xlabel('Simulation time / s'); ylabel('Current / A cm^{-2}');
legend('Unoptimized, 4 levels','Optimized, 5 levels', ...
    'Location','southeast');
grid on;
nexttile(layout);
plot(baseline.t,min(baseline.TavgC,[],2),'k--','LineWidth',1.4);
hold on;
plot(r.t,min(r.TavgC,[],2),'b','LineWidth',1.4);
yline(0,'r:'); xlabel('Simulation time / s');
ylabel('Minimum cell temperature / ^\circC'); grid on;
nexttile(layout);
plot(baseline.t,min(baseline.Vcell,[],2),'k--','LineWidth',1.4);
hold on;
plot(r.t,min(r.Vcell,[],2),'b','LineWidth',1.4);
yline(simOpt.minimumVoltageV,'r:');
xlabel('Simulation time / s'); ylabel('Minimum cell voltage / V');
grid on;
nexttile(layout);
plot(baseline.t,max(baseline.maxIceVolumeFractionCell,[],2), ...
    'k--','LineWidth',1.4);
hold on;
plot(r.t,max(r.maxIceVolumeFractionCell,[],2),'b','LineWidth',1.4);
xlabel('Simulation time / s'); ylabel('Maximum cell ice fraction');
grid on;
title(layout,'Step-loading comparison: same five-cell model and quick grid');
exportgraphics(fig,fullfile(codeDir,'Q2_step_loading_comparison.png'), ...
    'Resolution',180);
save(fullfile(codeDir,'Q2_step_loading_comparison.mat'), ...
    'baseline','optimized','baselineStrategy','optimizedStrategy', ...
    'simOpt','-v7.3');
