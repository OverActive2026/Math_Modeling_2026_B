function c=q3_cached_energy6(mode,q,th,opt)
% Cache physical candidates by all six variables, options and source versions.
    codeDir=fileparts(mfilename('fullpath'));
    outDir=fullfile(codeDir,'results');
    files={'pemfc_stack5_simulate_q3.m','thermal_temperature_state_ice.m', ...
        'water_ice_state_ice.m','gas_transport_state_ice.m','pemfc_setup_ice.m', ...
        'q3_auxiliary_heating.m','q3_evaluate_energy6.m','q2_current_strategy.m'};
    modified=zeros(size(files));
    for k=1:numel(files)
        info=dir(fullfile(codeDir,files{k})); modified(k)=info.datenum;
    end
    signature=struct('files',{files},'modified',modified,'options',opt);
    key=reshape(num2hex([q(:).',th]).',1,[]);
    path=fullfile(outDir,['energy6_' char(mode) '_' key '.mat']);
    if exist(path,'file')
        old=load(path,'c','signature');
        if isequaln(old.signature,signature) && ...
                ~strcmp(old.c.status,'numerical_unknown')
            c=old.c; return;
        end
    end
    % Reuse exact diagnostic points computed with the current split solver.
    probes=dir(fullfile(outDir,'q3_joint_probe_*.mat'));
    compareFields={'initialTemperatureC','ambientTemperatureC', ...
        'gammaIce','kFreeze','nCellLayer','RelTol','AbsTol','MaxStep', ...
        'iceScope','temperatureEventMarginK'};
    for k=1:numel(probes)
        if probes(k).datenum<max(modified), continue; end
        old=load(fullfile(outDir,probes(k).name),'candidate');
        p=old.candidate;
        if ~strcmp(mode,p.mode) || ~isequal(q(:).',p.q) || ...
                abs(th-p.th)>1e-10 || isempty(p.result) || ...
                strcmp(p.status,'numerical_unknown'), continue; end
        same=true;
        for j=1:numel(compareFields)
            f=compareFields{j};
            if ~isfield(opt,f) || ~isequaln(opt.(f),p.result.options.(f))
                same=false; break;
            end
        end
        if same && p.observationEndS>96.6666
            c=p; store(); return;
        end
    end
    c=q3_evaluate_energy6(mode,q,th,opt);
    store();
    function store()
        temporary=[tempname(outDir) '.mat'];
        save(temporary,'c','signature','-v7.3');
        movefile(temporary,path,'f');
    end
end
