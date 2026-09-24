%-------------------------------------------------------------------------------
% 第一版含冰模型调用脚本——电压损失分解绘图版
%
% 本文件是run_ice_model.m的独立副本，不修改原调用脚本。
%
% 使用方法：
%   1. 在下方修改活化参数、冰参数和计算工况；
%   2. 直接运行本脚本；
%   3. 新增图窗绘制Erev、etaAct、etaOhm、etaCon和lambda随时间的变化。
%
% 主要输出：
%   results.T20、results.T25
%   models.T20、models.T25
%   voltageLosses.T20、voltageLosses.T25
%-------------------------------------------------------------------------------

%% 1. 活化损失与水合反应面积参数
% [ACT-2][HYD-3] j0Ref为充分水合状态的参考值。

j0Ref = 0.527582014454;   % [REFIT-1] 298.15 K湿态参考值 [A/m^2]
Ea = 8000;                % [REFIT-1] 双温度整段电压联合标定 [J/mol]
fHydDry = 0.00705033754574; % [REFIT-1] 干态cCL可利用面积比例 [-]
lambdaHydOn = 3.5;        % 反应面积开始恢复的lambda [-]
lambdaHydWet = 8.53122413758; % [REFIT-1] 常温完全恢复阈值 [-]
nHyd = 3.0;               % 水合反应面积恢复指数 [-]
ionomerBruggemanExponent = 1.77972166611; % [REFIT-1] CL离聚物电导指数 [-]


%% 2. cCL离聚物水合参数

lambdaCCL0 = 2.0; % [INIT-2] cCL初始平均离聚物含水量 [-]
tauHyd = 10;      % cCL离聚物水合时间常数 [s]


%% 3. 两侧有限排水和CL/PEM界面参数

hVaporCathode = 0.01; % [WATER-2] cGDL/干空气气道传质系数 [m/s]
hVaporAnode = 0.01;   % [PHASE-2] aGDL/干氢气气道传质系数 [m/s]
interfaceRateFactor = 0.01; % [PHASE-1] CL/PEM界面速率倍率 [-]


%% 4. 三个冰模型标定参数
% 与当前run_ice_model.m保持一致；可在本副本中单独调整。

kFreeze = 0.05;   % 冻结非平衡速率系数 [1/s]，体积平均初值
kMelt = 1;        % [PHASE-5] 融化系数 [1/s]
directFreezeFraction = 0.01; % [PHASE-4] 剩余产水直接冻结比例 [-]
freezeNucleationLambdaFraction = 0.5; % 直接冻结开始的lambda/lambdaSat [-]
gammaIce = 0.1;   % 与主调用脚本保持一致 [-]


%% 5. 计算和绘图设置

temperatureCases = [-20,-25]; % 可改为-20、-25或[-20,-25]
tEnd = [];                     % []表示使用附件2的完整时间
showMainPlot = true;           % 原模型的电流、电压、温度和冰量图
showVoltageLossPlot = true;    % 新增的电压组成图
showInformation = true;        % true显示求解信息


%% 6. 整理求解选项

simOpt.j0Ref = j0Ref;
simOpt.Ea = Ea;
simOpt.fHydDry = fHydDry;
simOpt.lambdaHydOn = lambdaHydOn;
simOpt.lambdaHydWet = lambdaHydWet;
simOpt.nHyd = nHyd;
simOpt.ionomerBruggemanExponent = ionomerBruggemanExponent;
simOpt.tauHyd = tauHyd;
simOpt.hVaporCathode = hVaporCathode;
simOpt.hVaporAnode = hVaporAnode;
simOpt.interfaceRateFactor = interfaceRateFactor;
simOpt.init.lambdaCCL0 = lambdaCCL0;
simOpt.kFreeze = kFreeze;
simOpt.kMelt = kMelt;
simOpt.directFreezeFraction = directFreezeFraction;
simOpt.freezeNucleationLambdaFraction = freezeNucleationLambdaFraction;
simOpt.gammaIce = gammaIce;
simOpt.tEnd = tEnd;
simOpt.plot = showMainPlot;
simOpt.verbose = showInformation;


%% 7. 调用含冰模型

codeDir = fileparts(mfilename('fullpath'));
addpath(codeDir);

results = struct();
models = struct();

for temperatureC = temperatureCases
    caseName = sprintf('T%d',abs(temperatureC));
    fprintf('\n========== 计算%s degC工况 ==========\n',num2str(temperatureC));

    [results.(caseName),models.(caseName)] = ...
        pemfc_calculate_ice(temperatureC,simOpt);
end


%% 8. 汇总最终结果

nCase = numel(temperatureCases);
finalVoltage = zeros(nCase,1);
finalTemperatureC = zeros(nCase,1);
finalIceMassPerArea = zeros(nCase,1);
finalIceSaturationCCL = zeros(nCase,1);
finalIceVolumeFractionCCL = zeros(nCase,1);
maximumIceVolumeFraction = zeros(nCase,1);
finalHydrationAreaFactor = zeros(nCase,1);

for k = 1:nCase
    caseName = sprintf('T%d',abs(temperatureCases(k)));
    currentResult = results.(caseName);

    finalVoltage(k) = currentResult.Vcell(end);
    finalTemperatureC(k) = currentResult.TavgC(end);
    finalIceMassPerArea(k) = currentResult.iceMassPerArea(end);
    finalIceSaturationCCL(k) = currentResult.iceSaturationCCL(end);
    finalIceVolumeFractionCCL(k) = currentResult.iceVolumeFractionCCL(end);
    maximumIceVolumeFraction(k) = max(currentResult.maxIceVolumeFraction);
    finalHydrationAreaFactor(k) = currentResult.hydrationAreaFactor(end);
end

summaryTable = table(temperatureCases(:),finalVoltage,finalTemperatureC, ...
    finalIceMassPerArea,finalIceVolumeFractionCCL, ...
    maximumIceVolumeFraction,finalIceSaturationCCL, ...
    finalHydrationAreaFactor, ...
    'VariableNames',{'CaseTemperatureC','FinalVoltageV', ...
    'FinalTemperatureC','FinalIceMassKgM2', ...
    'FinalCCLIceVolumeFraction','MaximumIceVolumeFraction', ...
    'FinalCCLIceSaturation', ...
    'FinalHydrationAreaFactor'});

fprintf('\n========== 含冰模型计算汇总 ==========\n');
disp(summaryTable);


%% 9. 提取各项电压及损失
% 利用已保存的状态重新调用RHS，只读取代数输出，不重复积分。

voltageLosses = struct();

for k = 1:nCase
    caseName = sprintf('T%d',abs(temperatureCases(k)));
    currentResult = results.(caseName);
    currentModel = models.(caseName);
    nTime = numel(currentResult.t);

    Erev = zeros(nTime,1);
    etaAct = zeros(nTime,1);
    etaOhm = zeros(nTime,1);
    etaOhmPEM = zeros(nTime,1);
    etaOhmACL = zeros(nTime,1);
    etaOhmCCL = zeros(nTime,1);
    etaOhmContact = zeros(nTime,1);
    etaCon = zeros(nTime,1);

    for n = 1:nTime
        [~,rhsOut] = currentModel.rhs( ...
            currentResult.t(n),currentResult.x(n,:).');
        Erev(n) = rhsOut.voltage.Erev;
        etaAct(n) = rhsOut.voltage.etaAct;
        etaOhm(n) = rhsOut.voltage.etaOhm;
        etaOhmPEM(n) = rhsOut.voltage.etaOhmPEM;
        etaOhmACL(n) = rhsOut.voltage.etaOhmACL;
        etaOhmCCL(n) = rhsOut.voltage.etaOhmCCL;
        etaOhmContact(n) = rhsOut.voltage.etaOhmContact;
        etaCon(n) = rhsOut.voltage.etaCon;
    end

    voltageLosses.(caseName).time = currentResult.t;
    voltageLosses.(caseName).Erev = Erev;
    voltageLosses.(caseName).etaAct = etaAct;
    voltageLosses.(caseName).etaOhm = etaOhm;
    % [FIX-5] 保存欧姆损失的四个组成，便于定位快速升流误差。
    voltageLosses.(caseName).etaOhmPEM = etaOhmPEM;
    voltageLosses.(caseName).etaOhmACL = etaOhmACL;
    voltageLosses.(caseName).etaOhmCCL = etaOhmCCL;
    voltageLosses.(caseName).etaOhmContact = etaOhmContact;
    voltageLosses.(caseName).etaCon = etaCon;
    voltageLosses.(caseName).lambdaMean = currentResult.lambdaMean;
    voltageLosses.(caseName).lambdaCCL = currentResult.lambdaCCL;
    voltageLosses.(caseName).hydrationDegree = ...
        currentResult.hydrationDegree;
    voltageLosses.(caseName).hydrationAreaFactor = ...
        currentResult.hydrationAreaFactor;
    voltageLosses.(caseName).VcellReconstructed = ...
        Erev-etaAct-etaOhm-etaCon;
    voltageLosses.(caseName).maxReconstructionError = max(abs( ...
        voltageLosses.(caseName).VcellReconstructed-currentResult.Vcell));
end


%% 10. 新增绘图：可逆电压、三种电压损失和含水量
% 一张figure中使用五个子图；每个子图对比所选温度工况。

if showVoltageLossPlot
    figure('Name','PEMFC voltage and loss components');
    tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

    componentFields = {'Erev','etaAct','etaOhm','etaCon'};
    componentTitles = {'E_{rev}','\eta_{act}','\eta_{ohm}','\eta_{con}'};
    caseLegend = arrayfun(@(value) sprintf('%d ^\circC',value), ...
        temperatureCases,'UniformOutput',false);

    for componentIndex = 1:numel(componentFields)
        nexttile;
        hold on;

        for k = 1:nCase
            caseName = sprintf('T%d',abs(temperatureCases(k)));
            currentLoss = voltageLosses.(caseName);
            plot(currentLoss.time,currentLoss.(componentFields{componentIndex}), ...
                'LineWidth',1.5);
        end

        xlabel('Time / s');
        ylabel('Voltage / V');
        title(componentTitles{componentIndex});
        legend(caseLegend,'Location','best');
        grid on;
    end

    % 底部横跨两列：实线为PEM平均lambda，虚线为cCL离聚物lambda。
    nexttile([1,2]);
    hold on;
    plotColors = lines(nCase);
    lambdaLegend = cell(1,2*nCase);
    for k = 1:nCase
        caseName = sprintf('T%d',abs(temperatureCases(k)));
        currentLoss = voltageLosses.(caseName);
        plot(currentLoss.time,currentLoss.lambdaMean,'-', ...
            'Color',plotColors(k,:),'LineWidth',1.5);
        plot(currentLoss.time,currentLoss.lambdaCCL,'--', ...
            'Color',plotColors(k,:),'LineWidth',1.5);
        lambdaLegend{2*k-1} = sprintf('%d ^\circC PEM',temperatureCases(k));
        lambdaLegend{2*k} = sprintf('%d ^\circC cCL',temperatureCases(k));
    end
    xlabel('Time / s');
    ylabel('\lambda / -');
    title('PEM and cCL ionomer water content \lambda');
    legend(lambdaLegend,'Location','best');
    grid on;

    sgtitle('Reversible voltage, voltage losses and ionomer water content');
end
