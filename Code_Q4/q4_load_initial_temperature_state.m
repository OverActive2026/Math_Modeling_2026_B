%-------------------------------------------------------------------------------
% @function: 从MAT文件读取问题4预冷终态，供冷启动/优化模型作为初值
%-------------------------------------------------------------------------------

function initialState = q4_load_initial_temperature_state(fileName)

    if nargin < 1 || isempty(fileName) || ...
            ~(ischar(fileName) || (isstring(fileName) && isscalar(fileName)))
        error('q4_load_initial_temperature_state:InvalidFileName', ...
            'fileName必须是非空字符向量或字符串标量。');
    end
    fileName = char(fileName);
    if ~isfile(fileName)
        error('q4_load_initial_temperature_state:FileNotFound', ...
            '找不到预冷初值文件：%s',fileName);
    end

    data = load(fileName);
    if isfield(data,'initialState')
        initialState = data.initialState;
    elseif isfield(data,'precoolingResult') && ...
            isstruct(data.precoolingResult) && ...
            isfield(data.precoolingResult,'initialState')
        initialState = data.precoolingResult.initialState;
    else
        error('q4_load_initial_temperature_state:StateNotFound', ...
            ['文件中必须包含initialState，或包含带initialState字段的', ...
             'precoolingResult。']);
    end
    if ~isstruct(initialState) || ~isscalar(initialState)
        error('q4_load_initial_temperature_state:InvalidState', ...
            'initialState必须是标量结构体。');
    end

end
