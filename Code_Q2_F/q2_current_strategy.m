%-------------------------------------------------------------------------------
% @function: 问题2三类电流加载策略的统一实现
% @author:   PJ, DSH
% @date:     20260925
% @input:    t        -> 时间 [s]，可为向量
%            strategy -> 策略结构体
%              .type='constant': .currentAcm2
%                                .rampTimeS (可选) 软启动斜坡时长
%              .type='linear':   .initialAcm2(可选,默认0), .rampRateAcm2s,
%                                .plateauAcm2
%              .type='step':     .levelsAcm2, .switchTimesS
% @output:   j     -> 电流密度 [A/m^2]
%            jAcm2 -> 电流密度 [A/cm^2]
%
% [Q2F-2] 三种策略只在这一个函数里定义，使优化变量、物理方程和绘图代码
%   互不穿插（结构沿用 Code_Q2_PhaseFix 的 q2_current_strategy.m）。
%   函数**不**自动截断 jmax：是否越界由上层约束判定，避免把不可行策略
%   悄悄变成可行策略而掩盖"电流越界"这一真实失败模式。
% [Q2F-13] 恒流策略增加可选软启动斜坡：
%     j(t) = j_c * min(t/rampTimeS, 1)
%   即前 rampTimeS 秒从 0 线性升到 j_c，随后保持恒定。
%   依据：t=0 直接跳到 j_c 时，cCL 离聚物仍处于初始干态（lambda_CCL=1），
%   其质子电阻极大，j_c>=0.16 A/cm^2 会立刻把电压压到 0.30 V 以下而判失败
%   （见 README_Q2F.md 第3节）。给一个短斜坡让离聚物先水合，即可避免这种
%   "一上来就失败"。rampTimeS=0（默认）时严格退化为原来的阶跃恒流。
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
            rampTimeS = 0;
            if isfield(strategy,'rampTimeS') && ~isempty(strategy.rampTimeS)
                require_nonnegative_scalar(strategy,'rampTimeS');
                rampTimeS = strategy.rampTimeS;
            end
            if rampTimeS > 0
                jAcm2 = strategy.currentAcm2.*min(t./rampTimeS,1);
            else
                jAcm2 = strategy.currentAcm2+zeros(size(t));
            end

        case 'linear'
            require_nonnegative_scalar(strategy,'rampRateAcm2s');
            require_nonnegative_scalar(strategy,'plateauAcm2');
            initialAcm2 = 0;
            if isfield(strategy,'initialAcm2') && ~isempty(strategy.initialAcm2)
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
