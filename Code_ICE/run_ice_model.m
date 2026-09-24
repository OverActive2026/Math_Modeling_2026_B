%-------------------------------------------------------------------------------
% 第一版含冰模型调用脚本
%
% 使用方法：
%   1. 在下方“标定参数”区域修改参数；
%   2. 设置需要计算的温度工况；
%   3. 直接运行本脚本。
%
% 计算结果保存在：
%   results.T20、results.T25
% 对应模型信息保存在：
%   models.T20、models.T25
%-------------------------------------------------------------------------------

%% 1. 活化损失与水合反应面积参数
% [ACT-2][HYD-3] j0Ref为充分水合状态的参考值。

j0Ref = 0.0827422651388;  % 298.15 K湿态参考交换电流密度 [A/m^2]
Ea = 7028.980902;         % [ACT-2] 温度修正活化能 [J/mol]
fHydDry = 0.06;           % 干态cCL可利用反应面积比例 [-]
lambdaHydOn = 3.5;        % 反应面积开始恢复的lambda [-]
lambdaHydWet = 8.0;       % 反应面积完全恢复的lambda [-]
nHyd = 3.0;               % 水合反应面积恢复指数 [-]


%% 2. cCL离聚物水合参数
% [HYD-1] 吹扫后cCL比PEM更干；tauHyd控制产水进入离聚物的时间尺度。

lambdaCCL0 = 1.0; % cCL初始平均离聚物含水量 [-]
tauHyd = 10;      % cCL离聚物水合时间常数 [s]


%% 3. 阴极有限排水参数
% [WATER-2] cGDL到干空气气道的总体水蒸气传质系数。
% 附件1未给出气道尺寸和流量，因此该值需结合冰量进一步标定。

hVaporCathode = 0.01; % [m/s]；越小表示产水越难排出


%% 4. 三个冰模型标定参数
% 当前数值由用户自行设置；零值表示关闭对应冰机理。

kFreeze = 1;  % 冻结非平衡速率系数 [1/s]，当前采用文献量级
kMelt = 0;    % 融化系数 [1/s]，当前工况未越过冰点，暂不可辨识
gammaIce = 0.1; % 冰面积修正的约束拟合落在下界；仍保留接口供后续标定


%% 5. 计算设置

temperatureCases = [-20,-25]; % 可改为-20、-25或[-20,-25]
tEnd = [];                     % []表示使用附件2的完整时间
showPlot = true;               % true绘图，false不绘图
showInformation = true;        % true显示求解信息


%% 6. 整理求解选项

simOpt.j0Ref = j0Ref;
simOpt.Ea = Ea;
simOpt.fHydDry = fHydDry;
simOpt.lambdaHydOn = lambdaHydOn;
simOpt.lambdaHydWet = lambdaHydWet;
simOpt.nHyd = nHyd;
simOpt.tauHyd = tauHyd;
simOpt.hVaporCathode = hVaporCathode;
simOpt.init.lambdaCCL0 = lambdaCCL0;
simOpt.kFreeze = kFreeze;
simOpt.kMelt = kMelt;
simOpt.gammaIce = gammaIce;
simOpt.tEnd = tEnd;
simOpt.plot = showPlot;
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
