function plot_q2_temperature_boundary()
% Rebuild a paper-ready diagnostic plot and CSV from the saved temperature scan.
codeDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(codeDir,'optimization_results','fixed_optimized_temperature');
stored = load(fullfile(dataDir,'scan.mat'),'scan');
scan = stored.scan;
assert(scan.complete && numel(scan.tests)>=2);
successful = load(fullfile(dataDir,'T_-11.07031.mat'),'result');
failed = load(fullfile(dataDir,'T_-11.07812.mat'),'result');
s = successful.result;
f = failed.result;
assert(s.success && ~f.success);
rows = struct2table(scan.tests);
rows = removevars(rows,{'file','stopReason'});
rows = sortrows(rows,'temperatureC','descend');
writetable(rows,fullfile(dataDir,'temperature_scan_results.csv'));

fig = figure('Visible','off','Color','w','Position',[100 100 1100 760]);
tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

nexttile;
plot(s.t,s.TavgC(:,1),'b-','LineWidth',1.5); hold on;
plot(f.t,f.TavgC(:,1),'r--','LineWidth',1.5);
plot(s.t,s.TavgC(:,5),'c-','LineWidth',1.1);
plot(f.t,f.TavgC(:,5),'m--','LineWidth',1.1);
yline(0,'k:');
title('End-cell mean temperature'); xlabel('Time (s)'); ylabel('Temperature (°C)');
legend('Cell 1: success','Cell 1: failure', ...
    'Cell 5: success','Cell 5: failure','Location','northwest'); grid on;

nexttile;
plot(s.t,s.Vcell(:,1),'b-','LineWidth',1.5); hold on;
plot(f.t,f.Vcell(:,1),'r--','LineWidth',1.5);
yline(0.30,'k:','0.30 V limit');
title('Limiting cell voltage'); xlabel('Time (s)'); ylabel('Cell 1 voltage (V)');
legend('−11.0703125 °C: success','−11.078125 °C: failure', ...
    'Location','northeast'); grid on;

nexttile;
plot(s.t,s.maxIceVolumeFractionCell(:,1),'b-','LineWidth',1.5); hold on;
plot(f.t,f.maxIceVolumeFractionCell(:,1),'r--','LineWidth',1.5);
title('Limiting cell local ice fraction'); xlabel('Time (s)');
ylabel('Maximum local ice volume fraction'); grid on;

nexttile;
plot(s.t,s.currentDensityAcm2,'k-','LineWidth',1.5); hold on;
yline(0.5,'r:','0.5 A/cm² limit');
title('Fixed ramp-then-hold current'); xlabel('Time (s)');
ylabel('Current density (A/cm²)'); xlim([0 max(s.t)]); grid on;

print(fig,fullfile(dataDir,'temperature_boundary_comparison.png'), ...
    '-dpng','-r200');
close(fig);
fprintf('Saved temperature_boundary_comparison.png and temperature_scan_results.csv\n');
end
