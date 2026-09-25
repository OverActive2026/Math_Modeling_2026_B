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
% [F-1] 在 gammaIce=3.5、kFreeze=0.4 固定不变的前提下联合重标定。

Ea = 67000;         % [J/mol] 题目式(40) 给定值，不可标定
j0Ref = 0.28;       % [A/m^2] 标定（lambda_CCL=1、298.15 K 参考点）
mHyd  = 1.70;       % [-] [D-1] j0 水合指数（标定）
lamHydRef = 1.0;    % [-] [E-1] = lambdaCCL0 物理参考点，不参与标定
clResistanceScale = 0.66;  % [-] [FIX-5] cCL 质子电阻缩放（标定）


%% 2. cCL离聚物水合参数

lambdaCCL0 = 1.0; % cCL初始平均离聚物含水量 [-]
tauHyd = 40;      % cCL离聚物水合时间常数 [s]（标定）


%% 3. 低温产冰参数
% [PHASE-4] 就地直接冻结通道；[PHASE-4b] 已去掉把该通道压掉
% 三个数量级的"可动水门槛"，因此冰在本版中真正生成。
% [F-1] kFreeze 固定为 0.4（用户指定）；gammaIce=3.5 为附件1 给定值。

kFreeze = 0.4;     % [1/s] 题目式(9) 冻结非平衡速率系数（固定值）
kMelt = 0;         % [1/s] 当前工况未越过冰点，不可辨识
gammaIce = 3.5;    % [-] 附件1 给定：cCL冰覆盖活性面积指数
directFreezeFraction = 0.15;            % [-] 划给直接冻结通道的产水份额
freezeNucleationLambdaFraction = 0.25;  % [-] 成核门槛 lambda/lambda_sat
kFreezeDirect = 1; % [1/s] 直接冻结速率系数


%% 4. 计算和绘图设置

temperatureCases = [-20,-25]; % 可改为-20、-25或[-20,-25]
tEnd = [];                     % []表示使用附件2的完整时间
showMainPlot = true;           % 原模型的电流、电压、温度和冰量图
showVoltageLossPlot = true;    % 新增的电压组成图
showInformation = true;        % true显示求解信息
progressIntervalS = 1;         % 每 1 s 仿真时间打印一次求解进度；0=关闭


%% 5. 整理求解选项

simOpt.j0Ref = j0Ref;
simOpt.Ea = Ea;
simOpt.mHyd = mHyd;
simOpt.lamHydRef = lamHydRef;
simOpt.clResistanceScale = clResistanceScale;
simOpt.tauHyd = tauHyd;
simOpt.init.lambdaCCL0 = lambdaCCL0;
simOpt.kFreeze = kFreeze;
simOpt.kMelt = kMelt;
simOpt.gammaIce = gammaIce;
simOpt.directFreezeFraction = directFreezeFraction;
simOpt.freezeNucleationLambdaFraction = freezeNucleationLambdaFraction;
simOpt.kFreezeDirect = kFreezeDirect;
simOpt.tEnd = tEnd;
simOpt.plot = showMainPlot;
simOpt.verbose = showInformation;
simOpt.progressIntervalS = progressIntervalS;


%% 6. 调用含冰模型

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


%% 7. 汇总最终结果

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


%% 8. 提取各项电压及损失
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
    % 冰量与本版核心机制：冰覆盖面积修正因子 (1-sIce)^gammaIce
    voltageLosses.(caseName).iceSaturationCCL = currentResult.iceSaturationCCL;
    voltageLosses.(caseName).iceAreaFactor = currentResult.iceAreaFactor;
    voltageLosses.(caseName).iceMassPerArea = currentResult.iceMassPerArea;
    voltageLosses.(caseName).VcellReconstructed = ...
        Erev-etaAct-etaOhm-etaCon;
    voltageLosses.(caseName).maxReconstructionError = max(abs( ...
        voltageLosses.(caseName).VcellReconstructed-currentResult.Vcell));
end


%% 9. 新增绘图：可逆电压、三种电压损失和含水量
% 一张figure中使用五个子图；每个子图对比所选温度工况。

if showVoltageLossPlot
    figure('Name','PEMFC voltage and loss components');
    tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

    componentFields = {'Erev','etaAct','etaOhm','etaCon'};
    componentTitles = {'E_{rev}','\eta_{act}','\eta_{ohm}','\eta_{con}'};
    caseLegend = arrayfun(@(value) sprintf('%d degC',value), ...
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
        lambdaLegend{2*k-1} = sprintf('%d degC PEM',temperatureCases(k));
        lambdaLegend{2*k} = sprintf('%d degC cCL',temperatureCases(k));
    end
    xlabel('Time / s');
    ylabel('\lambda / -');
    title('PEM and cCL ionomer water content \lambda');
    legend(lambdaLegend,'Location','best');
    grid on;

    sgtitle('Reversible voltage, voltage losses and ionomer water content');
end


%% 10. 新增绘图：低温产冰过程
% 左：cCL 冰饱和度与冰面积修正因子；右：冰面密度增长。

if showVoltageLossPlot
    figure('Name','PEMFC ice formation during cold start');
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    caseLegend = arrayfun(@(value) sprintf('%d degC',value), ...
        temperatureCases,'UniformOutput',false);

    nexttile; hold on;
    for k = 1:nCase
        caseName = sprintf('T%d',abs(temperatureCases(k)));
        plot(voltageLosses.(caseName).time, ...
            voltageLosses.(caseName).iceSaturationCCL,'LineWidth',1.5);
    end
    xlabel('Time / s'); ylabel('s_{ice,cCL} / -');
    title('cCL ice saturation (题目式(9) 冰体积分数)');
    legend(caseLegend,'Location','best'); grid on;

    nexttile; hold on;
    for k = 1:nCase
        caseName = sprintf('T%d',abs(temperatureCases(k)));
        plot(voltageLosses.(caseName).time, ...
            voltageLosses.(caseName).iceAreaFactor,'LineWidth',1.5);
    end
    xlabel('Time / s'); ylabel('f_{area}=(1-s_{ice})^{\gamma} / -');
    title(sprintf('Ice area factor, \\gamma = %.1f (附件1)',gammaIce));
    legend(caseLegend,'Location','best'); grid on;

    nexttile; hold on;
    for k = 1:nCase
        caseName = sprintf('T%d',abs(temperatureCases(k)));
        plot(voltageLosses.(caseName).time, ...
            1e4*voltageLosses.(caseName).iceMassPerArea,'LineWidth',1.5);
    end
    xlabel('Time / s'); ylabel('Ice mass per area / (10^{-4} kg m^{-2})');
    title('Areal ice mass'); legend(caseLegend,'Location','best'); grid on;

    nexttile; hold on;
    for k = 1:nCase
        caseName = sprintf('T%d',abs(temperatureCases(k)));
        yyaxis left;
        plot(voltageLosses.(caseName).time, ...
            voltageLosses.(caseName).VcellReconstructed,'-','LineWidth',1.5);
        yyaxis right;
        plot(voltageLosses.(caseName).time, ...
            100*voltageLosses.(caseName).iceSaturationCCL,'--','LineWidth',1.5);
    end
    xlabel('Time / s'); ylabel('V_{cell} / V (left), s_{ice} / % (right)');
    title('Voltage vs ice saturation'); grid on;

    sgtitle('Low-temperature ice formation (Code\_ICE\_E)');
end


