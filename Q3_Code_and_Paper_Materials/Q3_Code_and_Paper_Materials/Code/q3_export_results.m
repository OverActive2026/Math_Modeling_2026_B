function summary = q3_export_results(preAudit,coAudit)
% Export task-defined startup results for two strategies on the same grid.
    if ~isequal(preAudit.result.options.nCellLayer, ...
            coAudit.result.options.nCellLayer)
        error('q3_export_results:GridMismatch', ...
            'Both strategies must use the same spatial grid.');
    end
    audits = {preAudit,coAudit};
    mode = strings(2,1);
    q = zeros(2,5);
    thRequestedS = zeros(2,1);
    thUsedS = zeros(2,1);
    ECellsJ = zeros(2,5);
    ETotalJ = zeros(2,1);
    startupTimeS = nan(2,1);
    success = false(2,1);
    assessmentScope = strings(2,1);
    minimumVoltageV = nan(2,1);
    maximumIceVolumeFraction = nan(2,1);
    stopReason = strings(2,1);
    for i = 1:2
        a = audits{i};
        r = a.result;
        mode(i) = string(a.mode);
        q(i,:) = a.q;
        thRequestedS(i) = a.thRequestedS;
        thUsedS(i) = r.actualHeatingTimeS;
        ECellsJ(i,:) = r.auxiliaryEnergyJCell;
        ETotalJ(i) = r.totalAuxiliaryEnergyJ;
        startupTimeS(i) = r.successTimeS;
        success(i) = r.success;
        if strcmp(a.mode,'preheat')
            assessmentScope(i) = "preheat_threshold_only";
        else
            assessmentScope(i) = "loaded_startup";
        end
        minimumVoltageV(i) = r.minimumCellVoltageV;
        maximumIceVolumeFraction(i) = r.maximumIceVolumeFraction;
        stopReason(i) = string(r.stopReason);
    end
    summary = table(mode,q(:,1),q(:,2),q(:,3),q(:,4),q(:,5), ...
        thRequestedS,thUsedS,ECellsJ(:,1),ECellsJ(:,2), ...
        ECellsJ(:,3),ECellsJ(:,4),ECellsJ(:,5),ETotalJ, ...
        startupTimeS,success,assessmentScope,minimumVoltageV,maximumIceVolumeFraction, ...
        stopReason, ...
        'VariableNames',{'mode','q1','q2','q3','q4','q5', ...
        'thRequestedS','thUsedS','E1J','E2J','E3J','E4J','E5J', ...
        'ETotalJ','startupTimeS','success','assessmentScope', ...
        'minimumVoltageV','maximumIceVolumeFraction','stopReason'});
    outdir = fullfile(fileparts(mfilename('fullpath')),'results');
    if ~exist(outdir,'dir'), mkdir(outdir); end
    writetable(summary,fullfile(outdir,'summary_question.csv'), ...
        'Encoding','UTF-8');
    disp(summary);
end
