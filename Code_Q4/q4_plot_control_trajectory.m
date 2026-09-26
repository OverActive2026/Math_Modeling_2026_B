function fig = q4_plot_control_trajectory( ...
    result,outputFile,figureTitle,showFigure)
%Q4_PLOT_CONTROL_TRAJECTORY Show/save a six-panel Q4 startup trajectory.
% The figure includes five heater powers, five cell temperatures/voltages,
% maximum local ice fraction per cell, cumulative auxiliary energy, and
% the instantaneous spread of cell-average temperatures.

    if nargin < 3 || isempty(figureTitle)
        figureTitle = 'Q4 dynamic auxiliary heating';
    end
    if nargin < 4 || isempty(showFigure)
        showFigure = usejava('desktop');
    end
    if ~isscalar(showFigure) || ~islogical(showFigure)
        error('q4_plot:InvalidVisibility', ...
            'showFigure must be true or false.');
    end
    if nargin < 2 || isempty(outputFile)
        error('q4_plot:MissingOutputFile', ...
            'Provide a PNG output filename.');
    end

    visibility = 'off';
    if showFigure, visibility = 'on'; end
    fig = figure('Visible',visibility,'Name',figureTitle, ...
        'Color','w','Position',[80 80 1320 800]);
    layout = tiledlayout(fig,3,2,'TileSpacing','compact', ...
        'Padding','compact');

    powerTimeS = [result.intervalStartS;result.stopTimeS];
    powerWcm2 = [result.powerDensityWcm2;zeros(1,5)];
    nexttile(layout);
    stairs(powerTimeS,powerWcm2,'LineWidth',1.3);
    ylim([0 1.05]); grid on;
    xlabel('Time / s'); ylabel('q_k / W cm^{-2}');
    title('Five heater powers');
    legend('Cell 1','Cell 2','Cell 3','Cell 4','Cell 5', ...
        'Location','best');

    nexttile(layout);
    plot(result.timeS,result.temperatureC,'LineWidth',1.3);
    yline(0,'k--','0 C'); grid on;
    xlabel('Time / s'); ylabel('Cell mean temperature / C');
    title('Temperature and startup threshold');

    nexttile(layout);
    plot(result.timeS,result.voltageV,'LineWidth',1.3);
    yline(result.options.minimumVoltageV,'k--','Voltage limit');
    grid on; xlabel('Time / s'); ylabel('Cell voltage / V');
    title('Five cell voltages');

    nexttile(layout);
    plot(result.timeS,result.iceVolumeFraction,'LineWidth',1.3);
    icePeak = max(result.iceVolumeFraction(:));
    ylim([0 max(0.15,min(1.02,1.1*icePeak))]);
    if icePeak > 0.85
        yline(result.options.maximumIceVolumeFraction, ...
            'k--','Ice limit');
    end
    grid on; xlabel('Time / s');
    ylabel('Local max ice volume fraction');
    title(sprintf('Ice fraction (limit %.2f)', ...
        result.options.maximumIceVolumeFraction));

    nexttile(layout);
    intervalDurationS = result.intervalEndS- ...
        result.intervalStartS;
    segmentEnergyJ = result.options.cellAreaM2*1e4* ...
        sum(result.powerDensityWcm2,2).*intervalDurationS;
    cumulativeEnergyJ = [0;cumsum(segmentEnergyJ)];
    plot([0;result.intervalEndS],cumulativeEnergyJ, ...
        'LineWidth',1.8);
    grid on; xlabel('Time / s');
    ylabel('Cumulative auxiliary energy / J');
    title(sprintf('Total %.2f J',result.totalAuxiliaryEnergyJ));

    nexttile(layout);
    spreadC = max(result.temperatureC,[],2)- ...
        min(result.temperatureC,[],2);
    plot(result.timeS,spreadC,'LineWidth',1.8);
    grid on; xlabel('Time / s');
    ylabel('Cell mean temperature spread / C');
    title(sprintf('Maximum spread %.2f C', ...
        result.maximumCellTemperatureSpreadC));

    outcome = 'FAIL';
    if result.success, outcome = 'SUCCESS'; end
    title(layout,sprintf('%s | %s | startup %.3f s | E %.2f J', ...
        figureTitle,outcome,result.stopTimeS, ...
        result.totalAuxiliaryEnergyJ),'Interpreter','none');
    exportgraphics(fig,outputFile,'Resolution',160);
    if ~showFigure, close(fig); end
end
