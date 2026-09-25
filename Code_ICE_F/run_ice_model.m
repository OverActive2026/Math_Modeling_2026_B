%-------------------------------------------------------------------------------
% Code_ICE_F 调用脚本 —— 固定 gammaIce=3.5、kFreeze=0.4 后重新联合标定
%
% 本目录从 Code_ICE_E 复制而来，E 版及其之前的目录均未修改。
% E 版的全部代码修正（[E-1]~[E-6]）在 F 版中原样保留，本版只改动标定参数。
%
% 与 E 版的区别（详见 README_F.md）：
%   [F-1] kFreeze 固定为 0.4（用户指定），不再作为可标定量
%         （附件1 的 gammaIce=3.5 在 E 版中已是固定值）
%   在 gammaIce=3.5、kFreeze=0.4 的约束下重新联合标定：
%         j0Ref 0.30 -> 0.28, mHyd 1.55 -> 1.70, clResistanceScale 0.64 -> 0.66
%
% 使用方法：
%   1. 直接运行本脚本（参数已固化为标定结果）；
%   2. 如需重新标定，修改下方"标定参数"区域。
%
% 计算结果保存在 results.T20、results.T25；
% 对应模型信息保存在 models.T20、models.T25。
%
% 标定结果（-20/-25 degC 两个工况全时段实验曲线联合拟合）：
%   -20 degC: V_RMSE = 0.01329 V, T_RMSE = 0.0497 K,
%             压降幅度达成率 101.39%，回升幅度达成率 104.62%
%   -25 degC: V_RMSE = 0.01437 V, T_RMSE = 0.1597 K,
%             压降幅度达成率 101.35%，回升幅度达成率 107.69%
%   两工况合并 RMSE = 0.01384 V
%     （E 版 0.01456 V，改善 5.0%；D 版 0.01903 V，改善 27%）
%   最大逐点相对误差(t=0,5,10,...):3.84% (-20 degC) / 3.82% (-25 degC)
%     （E 版同口径 4.69% / 4.22%）
%   末端冰量 1.78e-4 / 3.08e-4 kg/m2，cCL 冰饱和度 4.07% / 7.06%
%-------------------------------------------------------------------------------

%% 1. 题目给定值（不可标定）
% 依据：题目式(39) 电荷传递系数 alpha=0.5；式(40) 活化能 Ea=67000 J/mol；
% 附件1 阴极冰覆盖活性面积指数 3.5。

Ea = 67000;          % [J/mol] 题目式(40) 给定值
% alpha = 0.5 硬编码在 thermal_temperature_state_ice.m 中（题目式(39)）
gammaIce = 3.5;      % [-]     附件1 给定值


%% 2. 标定参数
% 说明：以下 3 个量是在 gammaIce=3.5、kFreeze=0.4 固定不变的前提下，
% 由两个工况的全部实测极化点联合重标定得到
% （标定目标为两工况合并 RMSE sqrt((RMSE20^2+RMSE25^2)/2)）。

j0Ref = 0.28;        % [A/m^2] lambda_CCL=1、T=298.15 K 的参考交换电流密度
mHyd  = 1.70;        % [-]     [D-1] j0 的离聚物水合指数
clResistanceScale = 0.66; % [-] [FIX-5] cCL 离聚物质子电阻缩放
tauHyd = 40;         % [s]     cCL 离聚物水合时间常数

% [E-1] lamHydRef 不是标定量：fHyd=(lambda/lamHydRef)^mHyd 中
% lamHydRef^mHyd 与 j0Ref 完全简并，故固定在物理参考点 lambdaCCL0=1.0，
% 使 fHyd(t=0)=1。
lamHydRef = 1.0;     % [-] = lambdaCCL0，物理参考点，不参与标定

lambdaCCL0 = 1.0;    % [-] cCL 初始平均离聚物含水量

% [D-1] 依据：实验反解的有效 j0 从 t=0 到 36.6 s 上升约 18.3 倍，而式(40)
% 只含温度项，Ea=67000 只能给 2.59 倍，缺约 7 倍；该缺口来自离聚物水合
% 导致的 Pt 有效三相界面面积增加。
% [FIX-5] cCL 离聚物质子电阻是"电流驱动"的欧姆项，负责前段压降；
% 水合因子负责中后段回升，两者作用在不同时间尺度上，因此可以解耦。


%% 3. 低温产冰参数
% [F-1] kFreeze 固定为 0.4（用户指定），gammaIce=3.5 为附件1 给定值。
% 这两个量固定后，冰的"量"和冰的"面积惩罚"都不再靠调参去吸收电压误差，
% 模型只能通过水合/欧姆参数拟合，因此 F 版的标定结果比 E 版更难"作弊"。

directFreezeFraction = 0.15;   % [-] 反应产水中划给直接冻结通道的份额
freezeNucleationLambdaFraction = 0.25; % [-] 成核门槛（lambda/lambda_sat）
kFreezeDirect = 1;             % [1/s] 直接冻结速率系数

% 题目式(9) 的冻结/融化关系保持启用；本组工况全程 T<0，
% 融化系数不可辨识，故取 0。
kFreeze = 0.4;       % [1/s] 冻结速率系数（固定值，见 [F-1]）
kMelt = 0;           % [1/s] 融化系数（全程 T<0，不可辨识）


%% 4. 相变与保水开关（保持 D 版口径，未启用）
kmCond = 0;          % [DSH-1] 有限速率冷凝；0 = 瞬时平衡（题目式(25)(26)）
DwFilmRef = 0;       % [B-1] cCL 离聚物膜水蒸气阻力；0 = 不启用
nSorp = 3;           % [C-1] 吸水等温线指数（本版未启用）
Kcov = 0;            % [C-2] 冰覆盖饱和常数；0 = 退化为 (1-s)^gamma
fRet = 1;            % [DSH-1a] 多孔层水蒸气扩散折减；1 = 原模型
hVaporAnode = 0;     % [PHASE-2] 0 = 干气道零浓度理想汇（原口径）
hVaporCathode = 0;   % [PHASE-2] 同上
phaseTransitionWidthK = 0; % [PHASE-5] 0 = 严格题目分段式（本工况 T 不跨冰点）
interfaceFactor = 30;      % [PHASE-1] CL/PEM 界面水交换倍率（原写死值）


%% 5. 计算设置

temperatureCases = [-20,-25]; % 可改为-20、-25或[-20,-25]
tEnd = [];                     % []表示使用附件2的完整时间
showPlot = true;               % true绘图，false不绘图
showInformation = true;        % true显示求解信息


%% 6. 整理求解选项

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
simOpt.kmCond = kmCond;
simOpt.DwFilmRef = DwFilmRef;
simOpt.nSorp = nSorp;
simOpt.Kcov = Kcov;
simOpt.fRet = fRet;
simOpt.hVaporAnode = hVaporAnode;
simOpt.hVaporCathode = hVaporCathode;
simOpt.phaseTransitionWidthK = phaseTransitionWidthK;
simOpt.interfaceFactor = interfaceFactor;
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
voltageRMSE = zeros(nCase,1);
temperatureRMSE = zeros(nCase,1);
dipAchievement = zeros(nCase,1);
recoveryAchievement = zeros(nCase,1);

for k = 1:nCase
    caseName = sprintf('T%d',abs(temperatureCases(k)));
    currentResult = results.(caseName);

    finalVoltage(k) = currentResult.Vcell(end);
    finalTemperatureC(k) = currentResult.TavgC(end);
    finalIceMassPerArea(k) = currentResult.iceMassPerArea(end);
    finalIceSaturationCCL(k) = currentResult.iceSaturationCCL(end);
    voltageRMSE(k) = currentResult.voltageRMSE;
    temperatureRMSE(k) = currentResult.temperatureRMSE;

    experiment = currentResult.experiment;
    dipAchievement(k) = 100*(interp1(currentResult.t,currentResult.Vcell,6.4)- ...
        interp1(currentResult.t,currentResult.Vcell,13.2))/ ...
        (interp1(experiment.time,experiment.voltage,6.4)- ...
        interp1(experiment.time,experiment.voltage,13.2));
    recoveryAchievement(k) = 100*(interp1(currentResult.t,currentResult.Vcell,24)- ...
        interp1(currentResult.t,currentResult.Vcell,15))/ ...
        (interp1(experiment.time,experiment.voltage,24)- ...
        interp1(experiment.time,experiment.voltage,15));
end

summaryTable = table(temperatureCases(:),finalVoltage,finalTemperatureC, ...
    finalIceMassPerArea,finalIceSaturationCCL,voltageRMSE,temperatureRMSE, ...
    dipAchievement,recoveryAchievement, ...
    'VariableNames',{'CaseTemperatureC','FinalVoltageV', ...
    'FinalTemperatureC','FinalIceMassKgM2','FinalCCLIceSaturation', ...
    'VoltageRMSEV','TemperatureRMSEK','DipAchievementPct','RecoveryAchievementPct'});

fprintf('\n========== 含冰模型计算汇总 ==========\n');
disp(summaryTable);
