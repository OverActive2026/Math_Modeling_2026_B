function plot_results(result,p)
%PLOT_RESULTS 四张独立结果图；温度仅绘图时转为摄氏度。
t = result.time;
if size(result.T_history,1) ~= numel(t) || ...
        size(result.T_history,2) ~= p.Ncell
    error('Q2:InvalidHistory','温度历史尺寸与时间或单片数不一致。');
end
labels = arrayfun(@(k)sprintf('Cell %d',k),1:p.Ncell, ...
    'UniformOutput',false);
suffix = '';
if result.is_test_placeholder
    suffix = ' (TEST PLACEHOLDER ONLY)';
end

figure('Name','Q2 current');
plot(t,result.current,'LineWidth',1.5);
xlabel('Time / s'); ylabel('Current density / A cm^{-2}');
title(['Current loading' suffix]); legend('j(t)','Location','best');
grid on;

figure('Name','Q2 temperatures');
plot(t,result.T_history-273.15,'LineWidth',1.2);
xlabel('Time / s'); ylabel('Temperature / ^\circC');
title(['Five-cell temperature' suffix]);
legend(labels,'Location','best'); grid on;

figure('Name','Q2 voltages');
plot(t,result.V_history,'LineWidth',1.2); hold on;
yline(p.V_min_limit,'k--','Voltage limit');
xlabel('Time / s'); ylabel('Voltage / V');
title(['Five-cell voltage' suffix]);
legend([labels,{'Voltage limit'}],'Location','best'); grid on;

figure('Name','Q2 ice');
plot(t,result.ice_history,'LineWidth',1.2); hold on;
yline(p.ice_limit,'k--','Ice limit');
xlabel('Time / s'); ylabel('Ice volume fraction');
title(['Five-cell maximum local ice fraction' suffix]);
legend([labels,{'Ice limit'}],'Location','best'); grid on;
end
