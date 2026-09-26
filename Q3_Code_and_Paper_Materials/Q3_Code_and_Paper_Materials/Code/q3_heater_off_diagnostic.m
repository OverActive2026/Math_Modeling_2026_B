function diagnostic=q3_heater_off_diagnostic()
% Local derivative at the verified startup state, with all heaters switched off.
% This diagnoses the immediate heat balance; it does not replace a trajectory.
    out=fullfile(fileparts(mfilename('fullpath')),'results');
    a=load(fullfile(out,'q3_final_energy_table.mat'),'verified');
    c=a.verified{2}; opt=c.result.options; opt.tMax=1e-6;
    opt.plot=false; opt.verbose=false;
    [~,model]=pemfc_stack5_simulate_q3('coheat',zeros(1,5),0,opt);
    [derivative,state]=model.rhs(c.startupTimeS,c.result.x(end,:).');
    rate=zeros(1,5);
    for k=1:5
        index=model.cellStateIndex{k}(model.s.idx.T);
        rate(k)=sum(derivative(index).*model.g.dx)/sum(model.g.dx);
    end
    diagnostic=struct('q',c.q,'startupTimeS',c.startupTimeS, ...
        'temperatureC',state.Tavg.'-273.15, ...
        'heaterOffTemperatureRateKperS',rate, ...
        'endPlateTemperatureC',state.endPlateTemperatureK.'-273.15, ...
        'endCellToEndPlateHeatW', ...
        [state.leftBoundaryOutflux,state.rightBoundaryOutflux]*opt.cellAreaM2);
    diagnostic.intercellHeatFlowW=state.interfaceHeatFlux.'*opt.cellAreaM2;
    diagnostic.reactionHeatW=zeros(1,5);
    diagnostic.phaseHeatW=zeros(1,5);
    diagnostic.cellHeatStorageRateW=zeros(1,5);
    for k=1:5
        thermal=state.cell{k}.thermal;
        index=model.cellStateIndex{k}(model.s.idx.T);
        diagnostic.reactionHeatW(k)=sum(thermal.qgen.*model.g.dx)*opt.cellAreaM2;
        diagnostic.phaseHeatW(k)=sum(thermal.qPhase.*model.g.dx)*opt.cellAreaM2;
        diagnostic.cellHeatStorageRateW(k)=sum( ...
            thermal.rhoCp.*derivative(index).*model.g.dx)*opt.cellAreaM2;
    end
    diagnostic.endPlateHeatStorageRateW=derivative(model.endPlateIndex).'* ...
        model.stackThermal.endPlateArealCapacity*opt.cellAreaM2;
    diagnostic.externalHeatLossW= ...
        [state.leftConvectionOutflux,state.rightConvectionOutflux]*opt.cellAreaM2;
    diagnostic.stackHeatBalanceResidualW= ...
        sum(diagnostic.cellHeatStorageRateW)+sum(diagnostic.endPlateHeatStorageRateW) ...
        -sum(diagnostic.reactionHeatW)-sum(diagnostic.phaseHeatW) ...
        +sum(diagnostic.externalHeatLossW);
    assert(abs(diagnostic.stackHeatBalanceResidualW)<1e-7, ...
        'Instantaneous stack heat balance failed.');
    save(fullfile(out,'q3_heater_off_diagnostic.mat'),'diagnostic');
    fid=fopen(fullfile(out,'q3_heater_off_diagnostic.json'),'w','n','UTF-8');
    cleaner=onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid,'%s',jsonencode(diagnostic,'PrettyPrint',true));
    disp(diagnostic);
end
