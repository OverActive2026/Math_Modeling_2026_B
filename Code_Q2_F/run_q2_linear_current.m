%-------------------------------------------------------------------------------
% 问题2 工况二：线性升载 j(t) = min(r*t, j_p)
%
% 直接运行即可。题目问题3 给出的 r=0.005 A/(cm^2 s)、j_p=0.3 A/cm^2
% 只作为初值，不应预先当成问题2 的最优值（见 README_Q2F.md）。
%
% 若要复现"最优升载"，先把 optimize 设为 true 做二维粗扫描。
%-------------------------------------------------------------------------------

clear; clc;
codeDir = fileparts(mfilename('fullpath'));
addpath(codeDir);

%% 1. 工况设置

rampRateAcm2s = 0.05;   % [A/(cm^2 s)] 升载速率 r（扫描最优）
plateauAcm2   = 0.50;   % [A/cm^2]    平台电流 j_p（扫描最优，等于 j_max）
initialAcm2   = 0;      % [A/cm^2]    起始电流（默认 0）
initialTemperatureC = -10;   % [degC] 题目2(1) 固定初始温度
optimize = false;        % true 时先做二维粗扫描

%% 2. 求解设置

simOpt = struct();
simOpt.tMax = 300;
simOpt.outputStep = 0.5;
simOpt.MaxStep = 0.05;
simOpt.RelTol = 1e-6;
simOpt.AbsTol = 1e-8;
simOpt.verbose = true;
simOpt.progressIntervalS = 1;    % 每 1 s 仿真时间打印一次求解进度（0=关闭）
simOpt.plot = true;

strategy = struct('type','linear','rampRateAcm2s',rampRateAcm2s, ...
    'plateauAcm2',plateauAcm2,'initialAcm2',initialAcm2);

%% 3. 可选：二维粗扫描

if optimize
    fprintf('\n========== 线性升载二维粗扫描（T0=%g degC）==========\n', ...
        initialTemperatureC);
    fprintf('%9s %9s | %9s | %11s | %9s | %s\n', ...
        'r','j_p','t_s/s','q_use','V_min','结果');
    scanOpt = simOpt; scanOpt.plot = false; scanOpt.verbose = false;
    scanOpt.tMax = 240; scanOpt.outputStep = 1.0;
    scanOpt.progressIntervalS = 0;   % 扫描时不打印进度
    bestScore = inf; bestR = NaN; bestP = NaN; bestTime = NaN;
    for rr = [0.02 0.03 0.04 0.05 0.06]
        for jp = [0.35 0.40 0.45 0.50]
            s = struct('type','linear','rampRateAcm2s',rr,'plateauAcm2',jp);
            r = pemfc_stack5_simulate(s,initialTemperatureC,scanOpt);
            if r.success
                fprintf('%9.4f %9.3f | %9.2f | %11.4f | %9.5f | 成功\n', ...
                    rr,jp,r.successTimeS,r.chargeUsedCcm2, ...
                    r.minimumCellVoltageV);
                score = r.successTimeS;
            else
                fprintf('%9.4f %9.3f | %9s | %11s | %9s | 失败\n', ...
                    rr,jp,'-','-','-');
                score = 1e6;
            end
            if score < bestScore
                bestScore = score; bestR = rr; bestP = jp;
                bestTime = r.successTimeS;
            end
        end
    end
    fprintf('扫描最优：r=%.4f A/(cm^2 s)，j_p=%.3f A/cm^2，t_s=%.3f s\n', ...
        bestR,bestP,bestTime);
    if isfinite(bestR)
        rampRateAcm2s = bestR; plateauAcm2 = bestP;
        strategy = struct('type','linear','rampRateAcm2s',rampRateAcm2s, ...
            'plateauAcm2',plateauAcm2,'initialAcm2',initialAcm2);
        simOpt.plot = true;
    end
end

%% 4. 正式计算

fprintf(['\n========== 工况二：线性升载 r=%.4f A/(cm^2 s)，', ...
    'j_p=%.3f A/cm^2 ==========\n'],rampRateAcm2s,plateauAcm2);
result = pemfc_stack5_simulate(strategy,initialTemperatureC,simOpt);

%% 5. 结果摘要（表3 行）

fprintf('\n---------- 表3（线性升载策略）----------\n');
if result.success
    fprintf('启动结果      : 成功\n');
    fprintf('启动时间 t_s  : %.3f s\n',result.successTimeS);
    fprintf('达到平台的时刻: %.3f s\n',plateauAcm2/rampRateAcm2s);
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

fprintf('\n逐片末态（第1/5片为端部片）：\n');
fprintf('%4s | %10s | %8s | %10s | %10s\n', ...
    '片号','T_end/degC','V_end/V','lambda_cCL','eps_ice,max');
for k = 1:size(result.TavgC,2)
    fprintf('%4d | %10.3f | %8.5f | %10.3f | %10.4g\n', ...
        k,result.TavgC(end,k),result.Vcell(end,k), ...
        result.lambdaCCL(end,k),result.iceVolumeFractionMax(end,k));
end
