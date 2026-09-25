%-------------------------------------------------------------------------------
% 问题2：五片电堆线性升载-恒流自冷启动快速测试
%
% 使用方法：
%   1. 修改下方初始电流、目标恒流和升载时间；
%   2. 直接运行本脚本；
%   3. 程序会运行到启动成功，或者电荷/电压/冰约束失败。
%
% 加载关系：
%   0<=t<rampTimeS: j=initialCurrentAcm2+
%                       (currentDensityAcm2-initialCurrentAcm2)*t/rampTimeS
%   t>=rampTimeS:   j=currentDensityAcm2
%
% 结果保存在result和model中，并绘制五片电池的
% 平均温度、单片电压和最大局部冰体积分数。
%-------------------------------------------------------------------------------

clearvars;
clc;

codeDir = fileparts(mfilename('fullpath'));
addpath(codeDir);


%% 1. 手动输入目标恒流密度和升载时间

% 固定对照组：该电流最初由前轮拉丁超立方初采样发现，现独立固定并复核。
% 与后续优化只比较这个固定加载在同一模型参数下的启动时间。
initialCurrentAcm2 = 0.02; % [A/cm^2]，线性升载的起始电流
currentDensityAcm2 = 0.259937538852887; % [A/cm^2]，固定目标恒流
rampTimeS = 2;             % [s]，从初始电流升到目标恒流的时间


%% 2. 初始温度和最长运行时间

initialTemperatureC = -10;     % 问题2(1)固定为-10 degC
maximumSimulationTimeS = 600;  % 极小电流工况的计算上限 [s]

if ~isscalar(currentDensityAcm2) || ~isfinite(currentDensityAcm2) || ...
        currentDensityAcm2 <= 0 || currentDensityAcm2 > 0.5
    error('run_q2_constant_current:InvalidCurrent', ...
        '恒定电流密度必须在(0,0.5] A/cm^2内。');
end
if ~isscalar(initialCurrentAcm2) || ...
        ~isfinite(initialCurrentAcm2) || initialCurrentAcm2 < 0 || ...
        initialCurrentAcm2 > currentDensityAcm2
    error('run_q2_constant_current:InvalidInitialCurrent', ...
        '初始电流密度必须位于[0,目标恒流]内。');
end
if ~isscalar(rampTimeS) || ~isfinite(rampTimeS) || rampTimeS <= 0
    error('run_q2_constant_current:InvalidRampTime', ...
        '线性升载时间必须是正有限标量。');
end


%% 3. 整理线性升载-恒流策略和求解设置

% 内核的linear策略在到达plateauAcm2后自动保持恒流。
rampRateAcm2s = ...
    (currentDensityAcm2-initialCurrentAcm2)/rampTimeS;
strategy = struct('type','linear', ...
    'initialAcm2',initialCurrentAcm2, ...
    'rampRateAcm2s',rampRateAcm2s, ...
    'plateauAcm2',currentDensityAcm2);

% 把升载段的三角形电荷也计入20 C/cm^2上限。
% 求解器仍会用事件在电荷刚好达到上限时终止。
chargeLimitCcm2 = 20;
chargeAtPlateauCcm2 = ...
    0.5*(initialCurrentAcm2+currentDensityAcm2)*rampTimeS;
if chargeLimitCcm2 <= chargeAtPlateauCcm2
    if rampRateAcm2s > 0
        chargeLimitTimeS = (-initialCurrentAcm2+sqrt( ...
            initialCurrentAcm2^2+2*rampRateAcm2s*chargeLimitCcm2))/ ...
            rampRateAcm2s;
    else
        chargeLimitTimeS = chargeLimitCcm2/currentDensityAcm2;
    end
else
    chargeLimitTimeS = rampTimeS+ ...
        (chargeLimitCcm2-chargeAtPlateauCcm2)/currentDensityAcm2;
end

simOpt.tMax = min(maximumSimulationTimeS,chargeLimitTimeS+1);
% simOpt.outputStep = 0.2;
% simOpt.MaxStep = 0.02;
% simOpt.RelTol = 1e-6;
simOpt.AbsTol = 1e-8;
% simOpt.plot = false;   % 本脚本在第5节统一绘图
% simOpt.verbose = true;

simOpt.nCellLayer = [2 3 4 4 2];
% 升载初期会经过接近零电流区，不再使用恒流快速
% 试算的1e-2松容差，避免Newton试探点越过物理状态。
simOpt.RelTol = 1e-4;
simOpt.MaxStep = 0.05;
simOpt.outputStep = 0.5;
simOpt.plot = false;
simOpt.verbose = false;
simOpt.progressIntervalS = 1; % 每隔1 s仿真时间打印当前进度；设0可关闭
simOpt.useJacobianPattern = true; % 五片电堆使用保守稀疏雅可比结构

% [Q2-2] 端板作为独立热节点，不再把整块端板热容
% 强制压到第1/5片电池上；5片堆的6块双极板只计数一次。
simOpt.endPlateModel = 'separate';
simOpt.bipolarPlateCounting = 'shared_stack';

% [INIT-2] 与问题1零时刻电压重标定保持一致。
simOpt.init.lambdaCCL0 = 2.0;
simOpt.j0Ref = 0.527582014454;
simOpt.Ea = 8000;

% [PHASE-1~5] 修正界面水循环、两侧排水、产水冻结和融化。
simOpt.hVaporCathode = 0.01;
simOpt.hVaporAnode = 0.01;
simOpt.interfaceRateFactor = 0.01;
simOpt.kFreeze = 0.3;
simOpt.gammaIce = 3.5;
simOpt.kMelt = 1;
% [NUM-PHASE] 仅平滑0 degC附近的冻结/融化开关；设0恢复严格分段式。
simOpt.phaseTransitionWidthK = 0.2;
simOpt.directFreezeFraction = 0.01;
simOpt.freezeNucleationLambdaFraction = 0.5;

% 题目约束
simOpt.qMaxCcm2 = 20;
simOpt.jMaxAcm2 = 0.5;
simOpt.minimumVoltageV = 0.30;
simOpt.maximumIceVolumeFraction = 0.99;
simOpt.maximumPoreOccupancy = 0.99; % 水+冰对孔隙的总占用率上限


%% 4. 运行到启动成功或失败

[result,model] = pemfc_stack5_simulate( ...
    strategy,initialTemperatureC,simOpt);

fprintf('\n========== 恒流启动结果 ==========\n');
fprintf('初始电流密度：%.4f A/cm^2\n',initialCurrentAcm2);
fprintf('目标恒流密度：%.4f A/cm^2\n',currentDensityAcm2);
fprintf('线性升载时间：%.4f s，升载速率：%.6f A/(cm^2 s)\n', ...
    rampTimeS,rampRateAcm2s);
fprintf('终止时间：%.4f s\n',result.stopTimeS);
fprintf('累积电荷：%.4f C/cm^2\n',result.chargeUsedCcm2);
fprintf('全程最低单片电压：%.4f V（第%d片）\n', ...
    result.minimumCellVoltageV,result.minimumVoltageCellIndex);
fprintf('全程最大局部冰体积分数：%.6f（第%d片）\n', ...
    result.maximumIceVolumeFraction,result.maximumIceCellIndex);

if result.success
    resultText = '启动成功';
    fprintf('启动结果：成功，ts=%.4f s\n',result.successTimeS);
else
    resultText = '启动失败';
    fprintf('启动结果：失败，原因：%s\n',result.stopReason);
    if result.stopTimeS < rampTimeS && ...
            result.minimumCellVoltageV < simOpt.minimumVoltageV
        fprintf(['提示：电压在线性升载段已越界，', ...
            '当时电流密度约为%.4f A/cm^2；', ...
            '可延长升载时间或降低目标恒流。\n'], ...
            result.currentDensityAcm2(end));
    end
end


%% 5. 绘制五片电池的温度、电压和冰体积分数

cellLegend = arrayfun(@(k)sprintf('第%d片',k),1:5, ...
    'UniformOutput',false);
lineColors = lines(5);
if numel(result.t) == 1
    markerStyle = 'o';
else
    markerStyle = 'none';
end

figure('Name','五片电堆线性升载-恒流自冷启动', ...
    'Color','w');
layout = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

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
yline(0.30,'k--','0.30 V下限','LineWidth',1.0);
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

% title(layout,sprintf( ...
%     '恒流 %.3f A/cm^2，T_0=%.1f ^\circC：%s', ...
%     currentDensityAcm2,initialTemperatureC,resultText));
