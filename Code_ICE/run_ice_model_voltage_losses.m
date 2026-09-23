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

%% 1. 活化损失参数
% [ACT-1] 由-20/-25 degC两个零时刻电压联合标定。

j0Ref = 0.00807579349527; % 298.15 K参考交换电流密度 [A/m^2]
Ea = 22204.0159133;       % 活化能 [J/mol]


%% 2. 三个冰模型标定参数
% 与当前run_ice_model.m保持一致；可在本副本中单独调整。

kFreeze = 10000;  % 冻结速率系数 [1/(K s)]
kMelt = 0;        % 融化速率系数 [1/(K s)]
gammaIce = 0;     % cCL冰覆盖对有效反应面积的修正指数 [-]


%% 3. 计算和绘图设置

temperatureCases = [-20,-25]; % 可改为-20、-25或[-20,-25]
tEnd = [];                     % []表示使用附件2的完整时间
showMainPlot = true;           % 原模型的电流、电压、温度和冰量图
showVoltageLossPlot = true;    % 新增的电压组成图
showInformation = true;        % true显示求解信息


%% 4. 整理求解选项

simOpt.j0Ref = j0Ref;
simOpt.Ea = Ea;
simOpt.kFreeze = kFreeze;
simOpt.kMelt = kMelt;
simOpt.gammaIce = gammaIce;
simOpt.tEnd = tEnd;
simOpt.plot = showMainPlot;
simOpt.verbose = showInformation;


%% 5. 调用含冰模型

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


%% 6. 汇总最终结果

nCase = numel(temperatureCases);
finalVoltage = zeros(nCase,1);
finalTemperatureC = zeros(nCase,1);
finalIceMassPerArea = zeros(nCase,1);
finalIceSaturationCCL = zeros(nCase,1);

for k = 1:nCase
    caseName = sprintf('T%d',abs(temperatureCases(k)));
    currentResult = results.(caseName);

    finalVoltage(k) = currentResult.Vcell(end);
    finalTemperatureC(k) = currentResult.TavgC(end);
    finalIceMassPerArea(k) = currentResult.iceMassPerArea(end);
    finalIceSaturationCCL(k) = currentResult.iceSaturationCCL(end);
end

summaryTable = table(temperatureCases(:),finalVoltage,finalTemperatureC, ...
    finalIceMassPerArea,finalIceSaturationCCL, ...
    'VariableNames',{'CaseTemperatureC','FinalVoltageV', ...
    'FinalTemperatureC','FinalIceMassKgM2','FinalCCLIceSaturation'});

fprintf('\n========== 含冰模型计算汇总 ==========\n');
disp(summaryTable);


%% 7. 提取各项电压及损失
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
    etaCon = zeros(nTime,1);

    for n = 1:nTime
        [~,rhsOut] = currentModel.rhs( ...
            currentResult.t(n),currentResult.x(n,:).');
        Erev(n) = rhsOut.voltage.Erev;
        etaAct(n) = rhsOut.voltage.etaAct;
        etaOhm(n) = rhsOut.voltage.etaOhm;
        etaCon(n) = rhsOut.voltage.etaCon;
    end

    voltageLosses.(caseName).time = currentResult.t;
    voltageLosses.(caseName).Erev = Erev;
    voltageLosses.(caseName).etaAct = etaAct;
    voltageLosses.(caseName).etaOhm = etaOhm;
    voltageLosses.(caseName).etaCon = etaCon;
    voltageLosses.(caseName).lambdaMean = currentResult.lambdaMean;
    voltageLosses.(caseName).VcellReconstructed = ...
        Erev-etaAct-etaOhm-etaCon;
    voltageLosses.(caseName).maxReconstructionError = max(abs( ...
        voltageLosses.(caseName).VcellReconstructed-currentResult.Vcell));
end


%% 8. 新增绘图：可逆电压、三种电压损失和膜含水量
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

    % 底部横跨两列绘制PEM平均含水量lambda。
    nexttile([1,2]);
    hold on;
    for k = 1:nCase
        caseName = sprintf('T%d',abs(temperatureCases(k)));
        currentLoss = voltageLosses.(caseName);
        plot(currentLoss.time,currentLoss.lambdaMean,'LineWidth',1.5);
    end
    xlabel('Time / s');
    ylabel('\lambda / -');
    title('Mean membrane water content \lambda');
    legend(caseLegend,'Location','best');
    grid on;

    sgtitle('Reversible voltage, voltage losses and membrane water content');
end
