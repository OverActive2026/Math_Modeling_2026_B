function fig = plot_q2_constant_fixed_comparison(report,outputPath)
% Compare the frozen reference with an independently verified optimum.
assert(report.verified && report.verifiedResult.success);
a = report.config.reference;
b = report.verifiedResult;
fig = figure('Name','Q2 constant loading comparison','Color','w', ...
    'Position',[100 100 1100 760]);
layout = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
nexttile(layout);
ta = [0:0.1:2 a.stopTimeS];
tb = [0:0.1:2 b.stopTimeS];
[~,ja] = q2_current_strategy(ta,a.strategy);
[~,jb] = q2_current_strategy(tb,b.strategy);
plot(ta,ja,'LineWidth',1.6); hold on;
plot(tb,jb,'LineWidth',1.6);
xlabel('Time (s)'); ylabel('Current density (A/cm^2)'); grid on;
legend('Fixed reference','FuRBO + LogCEI','Location','best');
title('2 s rise, then constant current');
nexttile(layout);
plot(a.t,min(a.TavgC,[],2),'LineWidth',1.6); hold on;
plot(b.t,min(b.TavgC,[],2),'LineWidth',1.6);
yline(0,'k--'); xlabel('Time (s)');
ylabel('Coldest cell mean temperature (degC)'); grid on;
title('Startup requires all five cells to reach 0 degC');
nexttile(layout);
plot(a.t,min(a.Vcell,[],2),'LineWidth',1.6); hold on;
plot(b.t,min(b.Vcell,[],2),'LineWidth',1.6);
yline(a.constraints.minimumVoltageV,'k--','Voltage limit');
xlabel('Time (s)'); ylabel('Lowest cell voltage (V)'); grid on;
title('Voltage throughout each simulation');
nexttile(layout);
y = report.search.Y;
y(~report.search.success) = Inf;
best = cummin(y); best(isinf(best)) = NaN;
stairs(1:numel(best),best,'LineWidth',1.6); hold on;
yline(a.successTimeS,'k--','Reference');
xline(report.config.nInit+0.5,':','Start adaptive search');
xlabel('Real model evaluation'); ylabel('Best successful startup time (s)');
grid on; title('Observed convergence');
title(layout,sprintf(['gammaIce=3.5, kFreeze=0.3 | %.4f s -> %.4f s ', ...
    '| reduction %.2f%% | original ode15s, coarse grid'], ...
    a.successTimeS,b.successTimeS,report.improvementPercent));
if nargin>=2 && ~isempty(outputPath)
    exportgraphics(fig,outputPath,'Resolution',160);
end
end
