function history=q3_export_joint_history()
% Export actual physical search records, not mock optimizer test results.
    outDir=fullfile(fileparts(mfilename('fullpath')),'results');
    saved=load(fullfile(outDir,'q3_joint_de_checkpoint.mat'),'run');
    h=saved.run.history;
    if isempty(h), history=table; return; end
    Z=vertcat(h.z);
    history=table((1:numel(h)).',[h.generation].', ...
        'VariableNames',{'evaluation','generation'});
    for k=1:5, history.(sprintf('q%d',k))=Z(:,k); end
    history.th=[h.th].'; history.status=string({h.status}).';
    history.energyJ=[h.energyJ].'; history.startupTimeS=[h.startupTimeS].';
    history.violation=[h.violation].'; history.accepted=[h.accepted].';
    feasibleEnergy=history.energyJ;
    feasibleEnergy(history.status~="success")=Inf;
    history.bestFeasibleEnergyJ=cummin(feasibleEnergy);
    history.bestFeasibleEnergyJ(isinf(history.bestFeasibleEnergyJ))=NaN;
    writetable(history,fullfile(outDir,'q3_joint_de_history.csv'),'Encoding','UTF-8');
end
