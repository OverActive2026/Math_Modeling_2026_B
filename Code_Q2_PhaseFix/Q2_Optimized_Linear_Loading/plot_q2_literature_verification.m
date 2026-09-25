function fig = plot_q2_literature_verification(verificationPath,outputPath)
%PLOT_Q2_LITERATURE_VERIFICATION Plot saved real-model verification only.
% This function performs no physical simulation or optimization.
codeDir=fileparts(mfilename('fullpath'));
if nargin<1 || isempty(verificationPath)
    verificationPath=fullfile(codeDir,'q2_literature_verification.mat');
end
if nargin<2 || isempty(outputPath)
    outputPath=fullfile(codeDir,'q2_literature_comparison.png');
end
s=load(verificationPath,'verification');
v=s.verification;
assert(v.passed,'Use a successful same-grid verification first.');
b=v.baseline.raw; o=v.optimized.raw;
fig=figure('Name','Q2 same-grid verified optimization', ...
    'Color','w','Position',[100 100 1000 680]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile;
draw_pair(b.t,b.currentDensityAcm2,o.t,o.currentDensityAcm2);
ylabel('Current density / A cm^{-2}');
title('Linear ramp and plateau');
legend('Unoptimized','Optimized','Location','southeast');
nexttile;
draw_pair(b.t,min(b.TavgC,[],2),o.t,min(o.TavgC,[],2));
yline(0,':','Startup threshold','HandleVisibility','off');
ylabel('Coldest cell mean temperature / degC');
title(sprintf('Startup: %.3f s -> %.3f s', ...
    v.baseline.startup_time,v.optimized.startup_time));
nexttile;
draw_pair(b.t,min(b.Vcell,[],2),o.t,min(o.Vcell,[],2));
yline(v.simOpt.minimumVoltageV,':r','Voltage limit', ...
    'HandleVisibility','off');
ylabel('Minimum cell voltage / V');
title('Five-cell voltage constraint');
nexttile;
draw_pair(b.t,max(b.maxIceVolumeFractionCell,[],2), ...
    o.t,max(o.maxIceVolumeFractionCell,[],2));
ylabel('Maximum local ice volume fraction');
title('Ice history (not experimental data)');
sgtitle(sprintf(['Same quick grid; gammaIce = %.1f; kFreeze = %.1f; ', ...
    'RelTol = %g; improvement %.2f%%'], ...
    v.simOpt.gammaIce,v.simOpt.kFreeze, ...
    v.simOpt.RelTol,v.improvementPercent));
exportgraphics(fig,outputPath,'Resolution',160);
end

function draw_pair(tb,yb,to,yo)
plot(tb,yb,'--','Color',[0.35 0.35 0.35],'LineWidth',1.4);
hold on;
plot(to,yo,'Color',[0 0.35 0.8],'LineWidth',1.6);
grid on;
xlabel('Simulation time / s');
xlim([0,max(tb(end),to(end))]);
end
