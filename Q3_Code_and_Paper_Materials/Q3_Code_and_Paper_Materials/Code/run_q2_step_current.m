%-------------------------------------------------------------------------------
% 问题2：五片电堆分段阶梯加载自冷启动快速测试
%
% 使用方法：
%   1. 修改第1节的电流档位currentLevelsAcm2和切换时刻switchTimesS；
%   2. 直接运行本脚本；
%   3. 程序会运行到五片均达到0 degC，或触发电荷/电压/冰约束。
%
% 示例加载关系：
%   0~5 s:   0.05 A/cm^2
%   5~12 s:  0.10 A/cm^2
%   12~20 s: 0.15 A/cm^2
%   20 s后:   0.20 A/cm^2
%
% currentLevelsAcm2长度必须比switchTimesS多1。
% 结果保存在result和model中，并绘制电流、五片温度、
% 电压和冰填充率。
%-------------------------------------------------------------------------------

clearvars;
clc;

codeDir = fileparts(mfilename('fullpath'));
addpath(codeDir);


%% 1. 手动输入分段阶梯加载参数

currentLevelsAcm2 = [0.05 0.10 0.25 0.30]; % 各阶段电流 [A/cm^2]
switchTimesS = [5 10 15];                  % 各次切换时刻 [s]


%% 2. 初始温度、最长运行时间和输入检查

initialTemperatureC = -10;     % 问题2(1)固定为-10 degC
maximumSimulationTimeS = 600;  % 低电流工况的计算上限 [s]

currentLevelsAcm2 = currentLevelsAcm2(:).';
switchTimesS = switchTimesS(:).';

if isempty(currentLevelsAcm2) || any(~isfinite(currentLevelsAcm2)) || ...
        any(currentLevelsAcm2 < 0) || any(currentLevelsAcm2 > 0.5)
    error('run_q2_step_current:InvalidCurrentLevels', ...
        '电流档位必须是[0,0.5] A/cm^2内的非空有限向量。');
end
if all(currentLevelsAcm2 == 0)
    error('run_q2_step_current:AllZeroCurrent', ...
        '至少需要一个大于0的电流档位。');
end
if numel(switchTimesS) ~= numel(currentLevelsAcm2)-1 || ...
        any(~isfinite(switchTimesS)) || any(switchTimesS <= 0) || ...
        any(diff(switchTimesS) <= 0)
    error('run_q2_step_current:InvalidSwitchTimes', ...
        ['switchTimesS长度必须比currentLevelsAcm2少1，', ...
         '并且切换时刻必须正且严格递增。']);
end


%% 3. 整理阶梯策略和求解设置

strategy = struct('type','step', ...
    'levelsAcm2',currentLevelsAcm2, ...
    'switchTimesS',switchTimesS);

% [STEP-1] 按每一个阶梯的“电流×持续时间”累加电荷，
% 解析定位20 C/cm^2上限所在的阶段，再多留1 s供事件定位。
chargeLimitCcm2 = 20;
chargeLimitTimeS = Inf;
chargeBeforeSegmentCcm2 = 0;
segmentStartTimeS = 0;

for k = 1:numel(switchTimesS)
    segmentDurationS = switchTimesS(k)-segmentStartTimeS;
    chargeAfterSegmentCcm2 = chargeBeforeSegmentCcm2+ ...
        currentLevelsAcm2(k)*segmentDurationS;
    if chargeLimitCcm2 <= chargeAfterSegmentCcm2
        chargeLimitTimeS = segmentStartTimeS+ ...
            (chargeLimitCcm2-chargeBeforeSegmentCcm2)/ ...
            currentLevelsAcm2(k);
        break;
    end
    chargeBeforeSegmentCcm2 = chargeAfterSegmentCcm2;
    segmentStartTimeS = switchTimesS(k);
end

if isinf(chargeLimitTimeS) && currentLevelsAcm2(end) > 0
    chargeLimitTimeS = segmentStartTimeS+ ...
        (chargeLimitCcm2-chargeBeforeSegmentCcm2)/currentLevelsAcm2(end);
end

simOpt.tMax = min(maximumSimulationTimeS,chargeLimitTimeS+1);

% 阶梯切换会产生不连续的代数电压和热源，限制最大步长
% 可以减小求解器跨过切换时刻时的数值偏差。
simOpt.nCellLayer = [2 3 4 4 2];
simOpt.RelTol = 1e-4;
simOpt.AbsTol = 1e-8;
simOpt.MaxStep = 0.05;
simOpt.outputStep = 0.5;
simOpt.plot = false;          % 本脚本在第5节统一绘图
simOpt.verbose = false;
simOpt.progressIntervalS = 1; % 每隔1 s仿真时间打印进度；设0可关闭
simOpt.useJacobianPattern = true;

% [Q2-2] 独立端板热节点，并避免重复计算双极板热容。
simOpt.endPlateModel = 'separate';
simOpt.bipolarPlateCounting = 'shared_stack';

% 与问题1零时刻电压标定保持一致。
simOpt.init.lambdaCCL0 = 2.0;
simOpt.j0Ref = 0.527582014454;
simOpt.Ea = 8000;

% 当前水循环和相变参数。
simOpt.hVaporCathode = 0.01;
simOpt.hVaporAnode = 0.01;
simOpt.interfaceRateFactor = 0.01;
simOpt.kFreeze = 0.1;
simOpt.kMelt = 1;
simOpt.phaseTransitionWidthK = 0.2;
simOpt.directFreezeFraction = 0.01;
simOpt.freezeNucleationLambdaFraction = 0.5;

% 题目约束。
simOpt.qMaxCcm2 = chargeLimitCcm2;
simOpt.jMaxAcm2 = 0.5;
simOpt.minimumVoltageV = 0.30;
simOpt.maximumIceVolumeFraction = 0.99;
simOpt.maximumPoreOccupancy = 0.98; % 水+冰占用率堵塞阈值，留2%数值裕度


%% 4. 运行到启动成功或失败

[result,model] = pemfc_stack5_simulate( ...
    strategy,initialTemperatureC,simOpt);

fprintf('\n========== 分段阶梯加载启动结果 ==========\n');
fprintf('电流档位 [A/cm^2]：');
fprintf(' %.4f',currentLevelsAcm2);
fprintf('\n切换时刻 [s]：');
fprintf(' %.4f',switchTimesS);
fprintf('\n终止时间：%.4f s，终止电流：%.4f A/cm^2\n', ...
    result.stopTimeS,result.currentDensityAcm2(end));
fprintf('累积电荷：%.4f C/cm^2\n',result.chargeUsedCcm2);
fprintf('全程最低单片电压：%.4f V（第%d片）\n', ...
    result.minimumCellVoltageV,result.minimumVoltageCellIndex);
fprintf('全程最大局部冰体积分数：%.6f（第%d片）\n', ...
    result.maximumIceVolumeFraction,result.maximumIceCellIndex);
fprintf('全程最大局部冰填充率：%.6f\n', ...
    result.maximumIceFillingRatio);

if result.success
    resultText = '启动成功';
    fprintf('启动结果：成功，ts=%.4f s\n',result.successTimeS);
else
    resultText = '启动失败';
    fprintf('启动结果：失败，原因：%s\n',result.stopReason);
end


%% 5. 绘制加载曲线及五片电池的温度、电压和冰填充率

cellLegend = arrayfun(@(k)sprintf('第%d片',k),1:5, ...
    'UniformOutput',false);
lineColors = lines(5);
if numel(result.t) == 1
    markerStyle = 'o';
else
    markerStyle = 'none';
end

figure('Name','五片电堆分段阶梯加载自冷启动','Color','w');
layout = tiledlayout(4,1,'TileSpacing','compact','Padding','compact');

nexttile;
yyaxis left;
stairs(result.t,result.currentDensityAcm2,'LineWidth',1.5);
ylabel('j / A cm^{-2}');
ylim([0,max(0.05,1.1*max(currentLevelsAcm2))]);
hold on;
for k = 1:numel(switchTimesS)
    xline(switchTimesS(k),':','LineWidth',0.8);
end
yyaxis right;
plot(result.t,result.chargeCcm2,'LineWidth',1.3);
yline(simOpt.qMaxCcm2,'k--','电荷上限');
ylabel('q / C cm^{-2}');
xlabel('时间 / s');
title('分段阶梯电流和累积电荷');
grid on;

nexttile;
temperatureLines = plot(result.t,result.TavgC, ...
    'LineWidth',1.4,'Marker',markerStyle,'MarkerSize',6);
for k = 1:5, temperatureLines(k).Color = lineColors(k,:); end
hold on;
yline(0,'k--','0 ^\circC启动线','LineWidth',1.0);
xlabel('时间 / s');
ylabel('平均温度 / ^\circC');
title('五片电池平均温度');
legend(temperatureLines,cellLegend,'Location','bestoutside');
grid on;

nexttile;
voltageLines = plot(result.t,result.Vcell, ...
    'LineWidth',1.4,'Marker',markerStyle,'MarkerSize',6);
for k = 1:5, voltageLines(k).Color = lineColors(k,:); end
hold on;
yline(simOpt.minimumVoltageV,'k--','电压下限','LineWidth',1.0);
xlabel('时间 / s');
ylabel('单片电压 / V');
title('五片电池电压');
grid on;

nexttile;
iceLines = plot(result.t,result.maxIceFillingRatioCell, ...
    'LineWidth',1.4,'Marker',markerStyle,'MarkerSize',6);
for k = 1:5, iceLines(k).Color = lineColors(k,:); end
hold on;
yline(1,'k--','孔隙完全被冰占据','LineWidth',1.0);
xlabel('时间 / s');
ylabel('最大局部冰填充率 s_i');
title('五片电池孔隙冰填充率');
ylim([0,1.05]);
legend(iceLines,cellLegend,'Location','bestoutside');
grid on;

title(layout,sprintf('分段阶梯加载，最终%.3f A/cm^2：%s', ...
    currentLevelsAcm2(end),resultText));
