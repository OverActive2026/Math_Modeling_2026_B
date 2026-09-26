%-------------------------------------------------------------------------------
% @function: 将均匀温度、预冷结构体或MAT文件统一解析为Q4温度初值
%-------------------------------------------------------------------------------

function [temperatureCellK,endPlateTemperatureK,info] = ...
    q4_resolve_initial_temperature(simOpt,g,nCells)

    hasState = isfield(simOpt,'initialState') && ...
        ~isempty(simOpt.initialState);
    hasStateFile = isfield(simOpt,'initialStateFile') && ...
        ~isempty(simOpt.initialStateFile);
    if hasState && hasStateFile
        error('q4_resolve_initial_temperature:AmbiguousState', ...
            'initialState和initialStateFile不能同时指定。');
    end

    if ~hasState && ~hasStateFile
        if ~isfield(simOpt,'initialTemperatureC') || ...
                ~isscalar(simOpt.initialTemperatureC) || ...
                ~isfinite(simOpt.initialTemperatureC) || ...
                simOpt.initialTemperatureC <= -273.15
            error('q4_resolve_initial_temperature:InvalidUniformState', ...
                '未提供预冷状态时，initialTemperatureC必须是有效标量。');
        end
        temperatureCellK = (simOpt.initialTemperatureC+273.15)* ...
            ones(g.N,nCells);
        endPlateTemperatureK = (simOpt.initialTemperatureC+273.15)* ...
            ones(2,1);
        info.kind = 'uniform_scalar';
        info.sourceFile = '';
        info.coolingDurationS = NaN;
        info.thermalModelId = '';
        return;
    end

    if hasStateFile
        state = q4_load_initial_temperature_state(simOpt.initialStateFile);
        sourceFile = char(simOpt.initialStateFile);
    else
        state = simOpt.initialState;
        sourceFile = '';
        if isstruct(state) && isscalar(state) && ...
                isfield(state,'initialState')
            state = state.initialState;
        end
    end
    if ~isstruct(state) || ~isscalar(state)
        error('q4_resolve_initial_temperature:InvalidState', ...
            '预冷初值必须是标量结构体。');
    end

    if isfield(state,'temperatureCellK')
        temperatureCellK = state.temperatureCellK;
    elseif isfield(state,'temperatureCellC')
        temperatureCellK = state.temperatureCellC+273.15;
    else
        error('q4_resolve_initial_temperature:MissingCellField', ...
            '预冷初值缺少temperatureCellK或temperatureCellC。');
    end
    if ~isequal(size(temperatureCellK),[g.N,nCells]) || ...
            any(~isfinite(temperatureCellK(:))) || ...
            any(temperatureCellK(:) <= 0)
        error('q4_resolve_initial_temperature:InvalidCellField', ...
            'temperatureCellK必须为%d x %d的正有限矩阵。',g.N,nCells);
    end

    if isfield(state,'endPlateTemperatureK')
        endPlateTemperatureK = state.endPlateTemperatureK(:);
    elseif isfield(state,'endPlateTemperatureC')
        endPlateTemperatureK = state.endPlateTemperatureC(:)+273.15;
    else
        error('q4_resolve_initial_temperature:MissingEndPlateField', ...
            '独立端板模型要求预冷初值包含两端板温度。');
    end
    if numel(endPlateTemperatureK) ~= 2 || ...
            any(~isfinite(endPlateTemperatureK)) || ...
            any(endPlateTemperatureK <= 0)
        error('q4_resolve_initial_temperature:InvalidEndPlateField', ...
            '端板温度必须包含两个正有限值。');
    end

    if isfield(state,'nCells') && state.nCells ~= nCells
        error('q4_resolve_initial_temperature:CellCountMismatch', ...
            '预冷初值的单电池数量与当前模型不一致。');
    end
    if isfield(state,'nCellLayer') && ...
            ~isequal(state.nCellLayer(:).',g.nCellLayer(:).')
        error('q4_resolve_initial_temperature:GridMismatch', ...
            '预冷初值的nCellLayer与当前模型不一致，不能直接传递。');
    end
    if isfield(state,'endPlateModel') && ...
            ~strcmpi(char(state.endPlateModel),char(simOpt.endPlateModel))
        error('q4_resolve_initial_temperature:EndPlateModelMismatch', ...
            '预冷初值与冷启动模型的端板模型不一致。');
    end
    if isfield(state,'bipolarPlateCounting') && ...
            ~strcmpi(char(state.bipolarPlateCounting), ...
            char(simOpt.bipolarPlateCounting))
        error('q4_resolve_initial_temperature:BipolarCountingMismatch', ...
            '预冷初值与冷启动模型的双极板热容计数口径不一致。');
    end

    info.kind = 'precooling_temperature_field';
    info.sourceFile = sourceFile;
    if isfield(state,'coolingDurationS')
        info.coolingDurationS = state.coolingDurationS;
    else
        info.coolingDurationS = NaN;
    end
    if isfield(state,'thermalModelId')
        info.thermalModelId = state.thermalModelId;
    else
        info.thermalModelId = '';
    end

end
