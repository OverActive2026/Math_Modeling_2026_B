function export_saved_trajectories()
% Read saved original-model solutions only. No ODE solve or optimization.
    root=fileparts(fileparts(mfilename('fullpath')));
    out=fullfile(root,'结果与图表');
    modes={'preheat','coheat','coheat_fullpower'};
    for k=1:3
        s=load(fullfile(out,sprintf('q3_strict_verification_%02d.mat',k)),'checkpoint');
        c=s.checkpoint; r=c.result;
        assert(strcmp(c.status,'success'));
        names={'time_s','current_Acm2','charge_Ccm2'};
        data=[r.t,r.currentDensityAcm2,r.chargeCcm2];
        fields={'TavgC','Vcell','maxIceVolumeFractionCell','cumulativeAuxiliaryEnergyJCell'};
        formats={'T%d_C','V%d_V','ice%d','E%d_J'};
        for j=1:numel(fields)
            data=[data,r.(fields{j})]; %#ok<AGROW>
            for i=1:5, names{end+1}=sprintf(formats{j},i); end %#ok<AGROW>
        end
        data=[data,r.endPlateTemperatureC];
        names=[names,{'leftEndplate_C','rightEndplate_C'}];
        writetable(array2table(data,'VariableNames',names),fullfile(out,[modes{k},'_trajectory.csv']));
        fig=figure('Visible','off','Color','w','Position',[100 100 1200 800]);
        tiledlayout(2,2);
        nexttile; plot(r.t,r.TavgC,'LineWidth',1.4); yline(0,'--');
        xlabel('Time (s)'); ylabel('Mean cell temperature (deg C)'); grid on;
        legend('Cell 1','Cell 2','Cell 3','Cell 4','Cell 5','Location','best');
        nexttile; plot(r.t,r.Vcell,'LineWidth',1.4); yline(.3,'--');
        xlabel('Time (s)'); ylabel('Cell voltage (V)'); grid on;
        nexttile; plot(r.t,r.maxIceVolumeFractionCell,'LineWidth',1.4);
        xlabel('Time (s)'); ylabel('Maximum ice volume fraction'); grid on;
        nexttile; plot(r.t,r.cumulativeAuxiliaryEnergyJCell,'LineWidth',1.4);
        xlabel('Time (s)'); ylabel('Auxiliary energy per cell (J)'); grid on;
        sgtitle(strrep(modes{k},'_',' '));
        exportgraphics(fig,fullfile(out,[modes{k},'_trajectories.png']),'Resolution',180);
        exportgraphics(fig,fullfile(out,[modes{k},'_trajectories.pdf']),'ContentType','vector');
        close(fig);
        fprintf('%s: %d saved samples, E=%.9f J, startup=%.9f s\n',modes{k},height(array2table(data)),r.totalAuxiliaryEnergyJ,r.successTimeS);
    end
end
