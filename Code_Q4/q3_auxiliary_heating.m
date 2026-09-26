%-------------------------------------------------------------------------------
% @function: 问题3恒功率电热丝模型
% @author:   PJ, GPT
% @date:     20260924
% @input:    t->时间 [s]，标量或向量
%            heating->辅助加热结构体
%              .powerDensityWcm2: 各片恒定功率密度 [W/cm^2]
%              .durationS: 加热持续时间 [s]
%              .maximumPowerDensityWcm2: 功率密度上限 [W/cm^2]
% @output:   powerDensityWcm2->每个时刻、每片电热丝功率密度
%
% 功率在0<=t<durationS内保持不变，随后同时关闭。输出的每一行对应
% 一个时刻，每一列对应一片单电池。面功率到体积热源的守恒转换由
% pemfc_stack5_simulate_q3完成。
%-------------------------------------------------------------------------------

function powerDensityWcm2 = q3_auxiliary_heating(t,heating)

    if ~isnumeric(t) || any(~isfinite(t(:))) || any(t(:) < 0)
        error('q3_auxiliary_heating:InvalidTime', ...
            't必须是非负有限数。');
    end
    if ~isstruct(heating) || ~isscalar(heating) || ...
            ~isfield(heating,'powerDensityWcm2') || ...
            ~isfield(heating,'durationS')
        error('q3_auxiliary_heating:InvalidHeating', ...
            'heating必须包含powerDensityWcm2和durationS。');
    end

    power = heating.powerDensityWcm2(:).';
    durationS = heating.durationS;
    if isempty(power) || any(~isfinite(power)) || any(power < 0)
        error('q3_auxiliary_heating:InvalidPower', ...
            '电热丝功率密度必须是非负有限向量。');
    end
    if ~isscalar(durationS) || ~isfinite(durationS) || durationS < 0
        error('q3_auxiliary_heating:InvalidDuration', ...
            '加热持续时间必须是非负有限标量。');
    end
    if isfield(heating,'maximumPowerDensityWcm2')
        maximumPower = heating.maximumPowerDensityWcm2;
        if ~isscalar(maximumPower) || ~isfinite(maximumPower) || ...
                maximumPower <= 0 || any(power > maximumPower)
            error('q3_auxiliary_heating:PowerLimitExceeded', ...
                '电热丝功率密度超过允许上限。');
        end
    end

    heaterActive = double(t(:) < durationS);
    powerDensityWcm2 = heaterActive*power;

end
