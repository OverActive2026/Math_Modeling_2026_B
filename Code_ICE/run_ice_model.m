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

%% 1. 活化损失参数
% [ACT-1] 由-20/-25 degC两个零时刻电压联合标定。

j0Ref = 0.2; % 298.15 K参考交换电流密度 [A/m^2]
Ea = 67000;       % 活化能 [J/mol]


%% 2. 三个冰模型标定参数
% 当前数值由用户自行设置；零值表示关闭对应冰机理。

kFreeze = 10000;     % 冻结速率系数 [1/(K s)]
kMelt = 0;       % 融化速率系数 [1/(K s)]
gammaIce = 0;     % cCL冰覆盖对有效反应面积的修正指数 [-]


%% 3. 计算设置

temperatureCases = [-20,-25]; % 可改为-20、-25或[-20,-25]
tEnd = [];                     % []表示使用附件2的完整时间
showPlot = true;               % true绘图，false不绘图
showInformation = true;        % true显示求解信息


%% 4. 整理求解选项

simOpt.j0Ref = j0Ref;
simOpt.Ea = Ea;
simOpt.kFreeze = kFreeze;
simOpt.kMelt = kMelt;
simOpt.gammaIce = gammaIce;
simOpt.tEnd = tEnd;
simOpt.plot = showPlot;
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
