%-------------------------------------------------------------------------------
% 问题2 工况一：恒流加载（带软启动斜坡）
%
%   j(t) = j_c * min(t / rampTimeS, 1)
%
% 即前 rampTimeS 秒从 0 线性升到 j_c，随后保持恒定。
% rampTimeS = 0 时退化为原来的阶跃恒流 j(t) = j_c。
%
% 为什么必须有软启动（[Q2F-13]）：
%   t=0 直接把电流跳到 j_c 时，cCL 离聚物仍处于初始干态（lambda_CCL = 1），
%   其质子电阻极大（-10 degC 下约 7 ohm cm^2），加上干态下极低的交换电流
%   密度，即使 j_c = 0.16 A/cm^2 也会立刻把电压压到 0.30 V 以下而判失败
%   （实测 V_min = 0.2809 V）。给一个几秒的斜坡让离聚物先水合，即可避免
%   这种"一上来就失败"。
%   扫描结果：阶跃（rampTimeS=0）在整个 j_c = 0.16~0.50 区间**全部失败**；
%   rampTimeS >= 5 s 起出现可行解。
%
% 直接运行即可；optimize=true 时对 (j_c, rampTimeS) 做二维扫描。
%-------------------------------------------------------------------------------

clear; clc;
codeDir = fileparts(mfilename('fullpath'));
addpath(codeDir);

%% 1. 工况设置

currentAcm2  = 0.26;   % [A/cm^2] 恒流电流密度（保持段）
rampTimeS    = 4;      % [s]      软启动斜坡时长；0 = 纯阶跃（会失败）
initialTemperatureC = -10;  % [degC] 题目2(1) 固定初始温度
optimize = false;      % true 时先做二维扫描再复核最优值

%% 2. 求解设置

simOpt = struct();
simOpt.tMax = 300;         % [s] 最长积分时间
simOpt.outputStep = 0.5;   % [s]
simOpt.MaxStep = 0.05;     % [s]
simOpt.RelTol = 1e-6;
simOpt.AbsTol = 1e-8;
simOpt.verbose = true;
simOpt.progressIntervalS = 1;    % 每 1 s 仿真时间打印一次求解进度（0=关闭）
simOpt.plot = true;        % 直接出图

strategy = struct('type','constant','currentAcm2',currentAcm2, ...
    'rampTimeS',rampTimeS);

%% 3. 可选：二维扫描找最优 (j_c, rampTimeS)

if optimize
    fprintf(['\n========== 恒流+软启动 二维扫描（T0=%g degC）==========\n' ...
        '说明：t_ramp=0 即纯阶跃恒流，用于对照"一上来就失败"。\n'], ...
        initialTemperatureC);
    fprintf('%8s %8s | %9s | %11s | %9s | %11s | %s\n', ...
        'j_c','t_ramp','t_s/s','q_use','j_max','V_min','结果');
    scanOpt = simOpt; scanOpt.plot = false; scanOpt.verbose = false;
    scanOpt.tMax = 300; scanOpt.outputStep = 1.0;
    scanOpt.progressIntervalS = 0;   % 扫描时不打印进度
    bestScore = inf; bestCurrent = NaN; bestRamp = NaN; bestTime = NaN;
    for tr = [0 4 5 6 8 10 12 15 20]
        for jc = 0.16:0.02:0.50
            s = struct('type','constant','currentAcm2',jc, ...
                'rampTimeS',tr);
            r = pemfc_stack5_simulate(s,initialTemperatureC,scanOpt);
            if r.success
                fprintf('%8.3f %8.1f | %9.2f | %11.4f | %9.4f | %11.5f | 成功\n', ...
                    jc,tr,r.successTimeS,r.chargeUsedCcm2, ...
                    r.maximumCurrentDensityAcm2,r.minimumCellVoltageV);
                score = r.successTimeS;
            else
                fprintf('%8.3f %8.1f | %9s | %11s | %9s | %11s | 失败\n', ...
                    jc,tr,'-','-','-','-');
                score = 1e6;
            end
            if score < bestScore
                bestScore = score; bestCurrent = jc; bestRamp = tr;
                bestTime = r.successTimeS;
            end
        end
    end
    fprintf(['\n扫描最优：j_c=%.3f A/cm^2，rampTimeS=%g s，t_s=%.3f s\n' ...
        '注意：j(t)=j_c*min(t/t_ramp,1) 与线性升载 min(r*t,j_p) 是同一函数族\n' ...
        '      （r = j_c/t_ramp）。若把 rampTimeS 完全放开，恒流的最优解会\n' ...
        '      退化到线性升载的最优解；要体现"恒流"，应限制 rampTimeS 足够小\n' ...
        '      （本题取 <= 5 s），即"短软启动 + 长恒定保持"。\n'], ...
        bestCurrent,bestRamp,bestTime);
    if isfinite(bestCurrent)
        currentAcm2 = bestCurrent; rampTimeS = bestRamp;
        strategy = struct('type','constant','currentAcm2',currentAcm2, ...
            'rampTimeS',rampTimeS);
        simOpt.plot = true;
    end
end

%% 4. 正式计算

fprintf(['\n========== 工况一：恒流 j_c=%.4f A/cm^2，' ...
    '软启动斜坡 %.3g s ==========\n'],currentAcm2,rampTimeS);
result = pemfc_stack5_simulate(strategy,initialTemperatureC,simOpt);

%% 5. 结果摘要（表3 行）

if rampTimeS > 0
    parameterText = sprintf('j_c=%.3f A/cm^2, 软启动 %.3g s', ...
        currentAcm2,rampTimeS);
else
    parameterText = sprintf('j_c=%.3f A/cm^2 (阶跃)',currentAcm2);
end
fprintf('\n---------- 表3（恒流策略）----------\n');
fprintf('加载参数      : %s\n',parameterText);
if result.success
    fprintf('启动结果      : 成功\n');
    fprintf('启动时间 t_s  : %.3f s\n',result.successTimeS);
else
    fprintf('启动结果      : 失败\n');
    for k = 1:numel(result.failureReasons)
        fprintf('  原因 %d      : %s\n',k,result.failureReasons{k});
    end
end
fprintf('累计电荷量    : %.4f C/cm^2   (上限 20)\n',result.chargeUsedCcm2);
fprintf('最大电流密度  : %.4f A/cm^2   (上限 0.5)\n',result.maximumCurrentDensityAcm2);
fprintf('最低单片电压  : %.5f V        (下限 0.30)\n',result.minimumCellVoltageV);
fprintf('最大冰体积分数: %.4g          (阈值 0.99)\n',result.maximumIceVolumeFraction);
fprintf(['启动判据温度  : %.3f degC  (逐时刻五片最低温的时间最大值；' ...
    '全程最低 %.3f degC = 初始温度)\n'], ...
    result.startupMinimumTemperatureC,result.minimumCellTemperatureC);
fprintf('约束余量      : q=%.3f C/cm^2, j=%.4f A/cm^2, V=%.4f V, ice=%.4f\n', ...
    result.margin.qMaxCcm2,result.margin.jMaxAcm2, ...
    result.margin.minimumVoltageV,result.margin.maximumIceVolumeFraction);

% 逐片末态，用于判断是否出现端部效应
fprintf('\n逐片末态（第1/5片为端部片）：\n');
fprintf('%4s | %10s | %8s | %10s | %10s\n', ...
    '片号','T_end/degC','V_end/V','lambda_cCL','eps_ice,max');
for k = 1:size(result.TavgC,2)
    fprintf('%4d | %10.3f | %8.5f | %10.3f | %10.4g\n', ...
        k,result.TavgC(end,k),result.Vcell(end,k), ...
        result.lambdaCCL(end,k),result.iceVolumeFractionMax(end,k));
end
