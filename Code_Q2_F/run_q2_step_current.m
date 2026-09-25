%-------------------------------------------------------------------------------
% 问题2 工况三：分段阶梯加载
%
%   j(t) = j1,  0    <= t < t1
%          j2,  t1   <= t < t2
%          j3,  t2   <= t
%
% 决策变量 (j1,j2,j3,t1,t2)。第一轮建议保持单调升载 j1<=j2<=j3、
% t1<t2；若发现高电流后减载能抑制冰堵，再放开单调限制。
%
% 直接运行即可；optimize=true 时对单调升载做三层粗扫描。
%-------------------------------------------------------------------------------

clear; clc;
codeDir = fileparts(mfilename('fullpath'));
addpath(codeDir);

%% 1. 工况设置

levelsAcm2   = [0.12 0.35 0.50];  % [A/cm^2] 三段电流密度 j1,j2,j3（扫描最优）
switchTimesS = [10 20];           % [s]      切换时刻 t1,t2（扫描最优）
initialTemperatureC = -10;        % [degC]   题目2(1) 固定初始温度
optimize = false;                 % true 时先做三层粗扫描

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

strategy = struct('type','step','levelsAcm2',levelsAcm2, ...
    'switchTimesS',switchTimesS);

%% 3. 可选：三段阶梯粗扫描（单调升载）

if optimize
    fprintf('\n========== 阶梯加载粗扫描（单调升载，T0=%g degC）==========\n', ...
        initialTemperatureC);
    fprintf('%6s %6s %6s %5s %5s | %9s | %11s | %9s | %s\n', ...
        'j1','j2','j3','t1','t2','t_s/s','q_use','V_min','结果');
    scanOpt = simOpt; scanOpt.plot = false; scanOpt.verbose = false;
    scanOpt.tMax = 240; scanOpt.outputStep = 1.0;
    scanOpt.progressIntervalS = 0;   % 扫描时不打印进度
    bestScore = inf; bestLevels = []; bestTimes = []; bestTime = NaN;
    levelGrid = [0.08 0.12 0.16 0.20 0.25 0.30 0.35 0.40 0.45 0.50];
    timeGrid = [5 10 15 20 30 40];
    for j1 = levelGrid
        for j2 = levelGrid(levelGrid >= j1)
            for j3 = levelGrid(levelGrid >= j2)
                for t1 = timeGrid
                    for t2 = timeGrid(timeGrid > t1)
                        s = struct('type','step', ...
                            'levelsAcm2',[j1 j2 j3], ...
                            'switchTimesS',[t1 t2]);
                        r = pemfc_stack5_simulate( ...
                            s,initialTemperatureC,scanOpt);
                        if r.success
                            fprintf('%6.3f %6.3f %6.3f %5d %5d | %9.2f | %11.4f | %9.5f | 成功\n', ...
                                j1,j2,j3,t1,t2,r.successTimeS, ...
                                r.chargeUsedCcm2,r.minimumCellVoltageV);
                            score = r.successTimeS;
                        else
                            score = 1e6;
                        end
                        if score < bestScore
                            bestScore = score;
                            bestLevels = [j1 j2 j3];
                            bestTimes = [t1 t2];
                            bestTime = r.successTimeS;
                        end
                    end
                end
            end
        end
    end
    fprintf('扫描最优：levels=[%s]，switches=[%s]，t_s=%.3f s\n', ...
        num2str(bestLevels),num2str(bestTimes),bestTime);
    if ~isempty(bestLevels)
        levelsAcm2 = bestLevels; switchTimesS = bestTimes;
        strategy = struct('type','step','levelsAcm2',levelsAcm2, ...
            'switchTimesS',switchTimesS);
        simOpt.plot = true;
    end
end

%% 4. 正式计算

fprintf('\n========== 工况三：分段阶梯 [%s] A/cm^2 @ [%s] s ==========\n', ...
    num2str(levelsAcm2),num2str(switchTimesS));
result = pemfc_stack5_simulate(strategy,initialTemperatureC,simOpt);

%% 5. 结果摘要（表3 行）

fprintf('\n---------- 表3（分段阶梯策略）----------\n');
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

fprintf('\n逐片末态（第1/5片为端部片）：\n');
fprintf('%4s | %10s | %8s | %10s | %10s\n', ...
    '片号','T_end/degC','V_end/V','lambda_cCL','eps_ice,max');
for k = 1:size(result.TavgC,2)
    fprintf('%4d | %10.3f | %8.5f | %10.3f | %10.4g\n', ...
        k,result.TavgC(end,k),result.Vcell(end,k), ...
        result.lambdaCCL(end,k),result.iceVolumeFractionMax(end,k));
end
