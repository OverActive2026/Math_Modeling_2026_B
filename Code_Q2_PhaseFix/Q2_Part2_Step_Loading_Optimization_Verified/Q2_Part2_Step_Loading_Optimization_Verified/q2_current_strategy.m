%-------------------------------------------------------------------------------
% @function: 将问题2的三类加载策略统一转换为电流密度
% @author:   PJ, GPT
% @date:     20260924
% @input:    t->时间 [s]
%            strategy->策略结构体
%              .type='constant': .currentAcm2
%              .type='linear':   .initialAcm2(可选), .rampRateAcm2s,
%                                .plateauAcm2
%              .type='step':     .levelsAcm2, .switchTimesS
% @output:   j->电流密度 [A/m^2]
%            jAcm2->电流密度 [A/cm^2]
%
% [Q2-3] 所有策略只在这一个函数中定义，避免优化变量、
%   物理方程和绘图代码互相穿插。函数不自动截断jmax，
%   由上层优化约束负责判定策略是否可行。
%-------------------------------------------------------------------------------

function [j,jAcm2] = q2_current_strategy(t,strategy)

    if ~isstruct(strategy) || ~isscalar(strategy) || ...
            ~isfield(strategy,'type')
        error('q2_current_strategy:InvalidStrategy', ...
            'strategy必须是含type字段的标量结构体。');
    end
    if ~isnumeric(t) || any(~isfinite(t(:))) || any(t(:) < 0)
        error('q2_current_strategy:InvalidTime', ...
            't必须是非负有限数。');
    end

    strategyType = lower(char(strategy.type));
    switch strategyType
        case 'constant'
            require_nonnegative_scalar(strategy,'currentAcm2');
            jAcm2 = strategy.currentAcm2+zeros(size(t));

        case 'linear'
            require_nonnegative_scalar(strategy,'rampRateAcm2s');
            require_nonnegative_scalar(strategy,'plateauAcm2');
            % 未给出初始电流时保持旧脚本行为：从0开始升载。
            initialAcm2 = 0;
            if isfield(strategy,'initialAcm2')
                require_nonnegative_scalar(strategy,'initialAcm2');
                initialAcm2 = strategy.initialAcm2;
            end
            if initialAcm2 > strategy.plateauAcm2
                error('q2_current_strategy:InitialAbovePlateau', ...
                    '线性加载的初始电流不能大于平台电流。');
            end
            jAcm2 = min(initialAcm2+strategy.rampRateAcm2s.*t, ...
                strategy.plateauAcm2);

        case 'step'
            if ~isfield(strategy,'levelsAcm2') || ...
                    ~isfield(strategy,'switchTimesS')
                error('q2_current_strategy:MissingStepParameters', ...
                    '分段阶梯需要levelsAcm2和switchTimesS。');
            end
            levels = strategy.levelsAcm2(:).';
            switchTimes = strategy.switchTimesS(:).';
            if isempty(levels) || any(~isfinite(levels)) || any(levels < 0)
                error('q2_current_strategy:InvalidStepLevels', ...
                    'levelsAcm2必须是非空、非负有限向量。');
            end
            if numel(switchTimes) ~= numel(levels)-1 || ...
                    any(~isfinite(switchTimes)) || any(switchTimes <= 0) || ...
                    any(diff(switchTimes) <= 0)
                error('q2_current_strategy:InvalidSwitchTimes', ...
                    ['switchTimesS长度应比levelsAcm2少1，', ...
                     '且必须严格递增并大于0。']);
            end

            jAcm2 = levels(1)+zeros(size(t));
            for k = 1:numel(switchTimes)
                jAcm2(t >= switchTimes(k)) = levels(k+1);
            end

        otherwise
            error('q2_current_strategy:UnknownStrategy', ...
                '未知加载策略%s。',strategyType);
    end

    j = 1e4*jAcm2;

end


function require_nonnegative_scalar(inputStruct,fieldName)

    if ~isfield(inputStruct,fieldName)
        error('q2_current_strategy:MissingParameter', ...
            '策略缺少参数%s。',fieldName);
    end
    value = inputStruct.(fieldName);
    if ~isscalar(value) || ~isfinite(value) || value < 0
        error('q2_current_strategy:InvalidParameter', ...
            '参数%s必须是非负有限标量。',fieldName);
    end

end
