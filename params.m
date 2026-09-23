%Suwren  检查于9.23 15：34
function p = params()
%PARAMS  一维 PEMFC 冷启动模型的参数集合
%
%   p = params()
%
%   参数分组说明
%   ----------------
%   p.const   : 通用物理常数
%   p.geom    : 几何参数
%   p.trans   : 气体传输参数
%   p.echem   : 电化学参数
%   p.thermal : 热相关参数
%   p.water   : 水 / 膜参数
%   p.ice     : 冷启动结冰模型参数
%   p.mat     : 各层材料物性
%   p.bc      : 边界条件参数
%
%   题面与附件 1 中明确给出的数值均已填入。NaN 表示仍需建模选择或标定
%   的参数，这些参数列在 p.meta.toBeCalibrated 中。
%
%   所有内部计算均采用国际单位制（SI）。

    %% ========================================================================
    %  1. 通用物理常数
    % ========================================================================
    
    p.const.R = 8.314;          % 通用气体常数 [J/(mol K)]
    p.const.F = 96485;          % 法拉第常数 [C/mol]
    
    p.const.Mw = 0.018;         % 水的摩尔质量 [kg/mol]
    p.const.MH2 = 2.016e-3;     % 氢气的摩尔质量 [kg/mol]
    p.const.MO2 = 31.998e-3;    % 氧气的摩尔质量 [kg/mol]
    p.const.MN2 = 28.014e-3;    % 氮气的摩尔质量 [kg/mol]
    
    p.const.p0   = 101325;      % 参考压力 [Pa]
    p.const.Tref = 298.15;      % 参考温度 [K]
    p.const.Tf   = 273.15;      % 水的凝固点温度 [K]
    
    
    %% ========================================================================
    %  2. 电池几何参数
    % ========================================================================
    %
    % 沿厚度方向的排列顺序：
    %
    %       aGDL | aCL | PEM | cCL | cGDL
    %
    
    p.geom.LaGDL = 150e-6;      % [m]
    p.geom.LaCL  = 3.4e-6;      % [m]
    p.geom.Lpem  = 12e-6;       % [m]
    p.geom.LcCL  = 11.3e-6;     % [m]
    p.geom.LcGDL = 150e-6;      % [m]
    
    p.geom.layerThickness = [ ...
        p.geom.LaGDL, ...
        p.geom.LaCL, ...
        p.geom.Lpem, ...
        p.geom.LcCL, ...
        p.geom.LcGDL ];
    
    p.geom.Lcell = sum(p.geom.layerThickness);
    p.geom.Lmea  = p.geom.Lcell;

    % 导热距离未必如此，这里先注释
    % % 附件 1 给出的其余电堆几何参数。
    % p.geom.LendPlate = 10e-3;   % 端板厚度 [m]
    % p.geom.LaBP      = 2e-3;    % 阳极双极板厚度 [m]
    % p.geom.LcBP      = 2e-3;    % 阴极双极板厚度 [m]

    % % 当两侧双极板均显式建模时的热传导总厚度。
    % p.geom.Lthermal = p.geom.LaBP + p.geom.Lmea + p.geom.LcBP;

    % 电池活化面积，来自附件 1。
    p.geom.Acell = 25e-4;       % 25 cm^2 -> [m^2]
    
    
    %% ========================================================================
    %  3. 气体传输
    % ========================================================================
    %
    % 题面给出的有效扩散系数：
    %
    % D_k,eff =
    %   D_k,ref * (T/298.15)^1.75 * (101325/p) * eps_g^1.5
    %
    
    p.trans.DH2_ref = 1.10e-4;  % H2 参考扩散系数 [m^2/s]
    p.trans.DO2_ref = 2.20e-5;  % O2 参考扩散系数 [m^2/s]
    
    p.trans.Texp = 1.75;        % 温度指数
    p.trans.epsExp = 1.5;       % 孔隙率指数
    
    % 操作压力，来自附件 1。
    % 这里实际上附件 1 只给出一个统一的压力值101325，未分别给出阳极、阴极
    p.bc.pAnode   = 101325;      % [Pa]
    p.bc.pCathode = 101325;      % [Pa]

    % 附件 1 给出的是质量分数，此处先按原值保存，再换算为摩尔分数，
    % 以供 c = y*p/(R*T) 使用。不要把 wO2=0.233 当作氧气摩尔分数。
    p.bc.wH2 = 1.0;
    p.bc.wO2 = 0.233;
    p.bc.wN2 = 0.767;
    p.bc.wVaporAnode   = 0.0;
    p.bc.wVaporCathode = 0.0;

    cathodeMoles = p.bc.wO2/p.const.MO2 + p.bc.wN2/p.const.MN2;
    p.bc.yH2 = 1.0;
    p.bc.yO2 = (p.bc.wO2/p.const.MO2) / cathodeMoles;
    p.bc.yN2 = (p.bc.wN2/p.const.MN2) / cathodeMoles;
    
    
    %% ========================================================================
    %  4. 水传输
    % ========================================================================
    
    % 多孔介质中水的参考扩散系数
    p.water.Dw_ref_anode   = 8.69e-5;  % aGDL / aCL [m^2/s]
    p.water.Dw_ref_cathode = 2.48e-5;  % cGDL / cCL [m^2/s]
    
    % 电渗拖拽（electro-osmotic drag）：
    %
    %       nd = 2.5 * lambda / 22
    %
    p.water.ndCoeff = 2.5 / 22;
    
    % 膜内水扩散系数由题目明确给出，将在 calc_Dwater.m 中实现：
    %
    % Dw_mem =
    % 1e-10 * exp[2416(1/303.15 - 1/T)] *
    % (2.563 - 0.33*lambda + 0.0264*lambda^2 - 0.000671*lambda^3)
    
    p.water.Dmem_prefactor = 1e-10;    % [m^2/s]
    p.water.Dmem_Eterm     = 2416;
    p.water.Dmem_c = [ ...
        2.563, ...
       -0.33, ...
        0.0264, ...
       -0.000671 ];
    
    % 膜当量质量与膜密度来自附件 1，用于下式：
    %
    % lambda = EW * mw / (rho_pem * Mw)
    %
    % 1000 g/mol = 1 kg/mol.
    p.water.EW      = 1.0;       % [kg/mol]
    p.water.rho_pem = 2150;      % [kg/m^3]
    p.water.lambda0 = 3.0;       % 初始膜含水量 [-]
    p.water.clIonomerVolumeFraction = 0.3;

    % 由式 (21) 推导出的初始膜内水质量浓度。
    p.water.mwMembrane0 = p.water.lambda0 * p.water.rho_pem * ...
        p.const.Mw / p.water.EW; % [kg/m^3]

    % 多孔介质参数，保留给后续的液态水模型使用。
    p.porous.contactAngleGDL_deg = 110;
    p.porous.contactAngleCL_deg  = 100;
    p.porous.contactAngleGDL = 110*pi/180; % [rad]
    p.porous.contactAngleCL  = 100*pi/180; % [rad]
    p.porous.permeabilityGDL = 6.2e-12;    % [m^2]
    p.porous.permeabilityCL  = 6.2e-13;    % [m^2]

    % 相变 / 相间转化系数，题目以无量纲模型系数的形式给出。
    % 附件 1 未说明成对数值的方向顺序，因此成对数值保持原样存放，
    % 不做方向上的猜测性分配。
    p.water.phase.memVaporCoeff  = [0.001, 1.0];
    p.water.phase.memLiquidCoeff = 0.5;
    p.water.phase.memIceCoeff    = 1.0;
    p.water.phase.vaporLiquidCoeff = [1.0, 1.0];
    p.water.phase.vaporIceCoeff  = 1.0e-4;
    
    
    %% ========================================================================
    %  5. 电化学参数
    % ========================================================================
    
    p.echem.alpha = 0.5;          % 电荷转移系数

    p.echem.Voc0 = 0.95;          % 初始开路电压 [V]
    
    % 题面给出的初始校准值
    p.echem.j0_ref = 0.01;        % [A/m^2]
    
    p.echem.Ea = 67000;           % 活化能 [J/mol]
    
    % 热中性电压
    p.echem.Eth = 1.48;           % [V]
    
    % 可逆电压常数：
    %
    % Erev =
    % 1.229 - 8.5e-4*(T-298.15)
    % + RT/(2F) ln[ (pH2/p0)*(pO2/p0)^0.5 ]
    %
    p.echem.Erev_ref = 1.229;     % [V]
    p.echem.dEdT     = -8.5e-4;   % [V/K]
    
    % 单位面积接触电阻
    %
    % 题面给出：
    %       Rc = 0.01 ohm cm^2
    %
    % 单位换算：
    %       1 cm^2 = 1e-4 m^2
    %
    p.echem.Rc = 0.01e-4;         % [ohm m^2]
    
    % 冰覆盖指数与电堆浓差极化系数，来自附件 1。
    % 在无冰的单电池最小可用版本中暂不生效。
    p.echem.gammaIceArea = 3.5;
    p.echem.concentrationFactorEnd = 10;
    p.echem.concentrationFactorMiddle = 1;

    % 材料不存在这个公式，Ldiff这个变量不知道是什么
    % % 阴极等效扩散距离。针对采用单元中心值的 cCL 浓度，取阴极边界到 cCL
    % % 代表性中心的距离：阴极 GDL 厚度 + 阴极 CL 厚度的一半。
    % p.echem.Ldiff = p.geom.LcGDL + 0.5*p.geom.LcCL; % [m]
    
    
    %% ========================================================================
    %  6. 热参数
    % ========================================================================
    
    % 对流换热系数
    p.thermal.h = 40;             % [W/(m^2 K)]
    
    % 材料物性，来自附件 1。
    %
    % 排列顺序：
    %       [aGDL, aCL, PEM, cCL, cGDL]
    %
    % rho : 密度         [kg/m^3]
    % cp  : 比热容       [J/(kg K)]
    % k   : 导热系数     [W/(m K)]
    
    p.mat.rho = [185, 970, 2150, 970, 185];
    p.mat.cp  = [545, 240, 1050, 240, 545];
    p.mat.k   = [0.3, 0.27, 0.24, 0.27, 0.3];

    % 未直接映射到五层传输网格上的部件物性。
    p.mat.GDL.rho = 185;          % [kg/m^3]
    p.mat.GDL.cp = 545;           % [J/(kg K)]
    p.mat.GDL.k = 0.3;            % [W/(m K)]
    p.mat.GDL.sigma = 375;        % [S/m]

    p.mat.CL.rho = 970;
    p.mat.CL.cp = 240;
    p.mat.CL.k = 0.27;

    p.mat.ionomer.rho = 2150;
    p.mat.ionomer.cp = 1050;
    p.mat.ionomer.k = 0.24;

    p.mat.BP.rho = 1980;
    p.mat.BP.cp = 766;
    p.mat.BP.k = 95;
    p.mat.BP.sigma = 83000;

    p.mat.endPlate.rho = 7900;
    p.mat.endPlate.cp = 500;
    p.mat.endPlate.k = 15;

    p.phase.vapor.rho = 4.8e-3;
    p.phase.vapor.cp = 2000;
    p.phase.vapor.k = 0.1;

    p.phase.ice.rho = 920;
    p.phase.ice.cp = 2050;
    p.phase.ice.k = 2.3;

    p.phase.liquid.rho = 990;
    p.phase.liquid.cp = 4182;
    p.phase.liquid.k = 0.6;

    p.gas.H2.rho = 0.089;
    p.gas.H2.cp = 14283;
    p.gas.H2.k = 0.1672;

    p.gas.O2.rho = 1.43;
    p.gas.O2.cp = 919.31;
    p.gas.O2.k = 0.0246;

    p.gas.N2.rho = 1.35;
    p.gas.N2.cp = 1041.5;
    p.gas.N2.k = 0.0235;

    p.thermal.Lcondensation = 2.5e6; % [J/kg]
    p.thermal.Lfusion = 333600;       % [J/kg]
    
    
    %% ========================================================================
    %  7. 多孔层的干孔隙率
    % ========================================================================
    %
    % eps0 用于下式：
    %
    %       eps_g = eps0 - eps_l - eps_ice
    %
    % 进而影响气体 / 水的有效扩散系数。
    %
    % 在本简化处理中，质子交换膜按无孔隙对待。
    %
    % 附件 1 没有给膜的孔隙率，这里给了0
    p.mat.eps0 = [ ...
        0.8,    ...     % aGDL
        0.3916, ...     % aCL
        0,   ...       % PEM
        0.4207, ...     % cCL
        0.8 ];          % cGDL
    
    
    %% ========================================================================
    %  8. 结冰 / 冷启动模型
    % ========================================================================
    %
    % 这些不属于常温基线模型的内容！！！
    % 而是我们自行增加的冷启动扩展，最终需要用 -20 ℃ / -25 ℃ 实验数据标定。
    %
    
    % 密度、潜热与活性面积指数来自附件 1。
    p.ice.rho_liquid = p.phase.liquid.rho; % [kg/m^3]
    p.ice.rho_ice    = p.phase.ice.rho;    % [kg/m^3]
    p.ice.Lfusion    = p.thermal.Lfusion;  % [J/kg]
    
    % 初步提出的动力学规律：
    %
    % rFreeze = kFreeze * ml * max(Tf - T,0)
    % rMelt   = kMelt   * mi * max(T - Tf,0)
    %
    % 这两个是模型参数，题目并未给出。
    
    p.ice.kFreeze = NaN;          % 待标定参数
    p.ice.kMelt   = NaN;          % 待标定参数
    
    % 冰堵对电化学有效反应面积的影响：
    %
    %       a_eff = (1 - sIce)^gammaArea
    %
    p.ice.gammaArea = p.echem.gammaIceArea;
    
    % 剩余气相孔隙率的数值下限。
    % 这属于数值保护措施，不是物理参数。
    p.num.epsMin = 1e-8;
    p.num.residualTol = 1e-6;

    % 初始条件与环境温度，来自附件 1。
    p.ic.T0 = 253.15;             % [K]
    p.bc.Tamb = 253.15;           % [K]
    p.ic.liquidWater = 0.0;
    p.ic.iceVolumeFraction = 0.0;
    
    
    %% ========================================================================
    %  9. 膜电导率关联式
    % ========================================================================
    %
    % kappa_pem =
    % (0.5139 lambda - 0.326) *
    % exp[1268(1/303.15 - 1/T)]
    %
    % 该式无需拟合参数，但系数统一存放在此处，
    % 避免在电压计算函数中硬编码。
    
    p.water.kappa_a0 = 0.5139;
    p.water.kappa_a1 = -0.326;
    p.water.kappa_E  = 1268;
    
    
    %% ========================================================================
    %  10. Buck 饱和蒸汽压关联式
    % ========================================================================
    %
    % Tc = T - 273.15 [摄氏度]
    %
    % Tc >= 0:
    % psat = 611.21 exp[(18.678 - Tc/234.5)*Tc/(257.14 + Tc)]
    %
    % Tc < 0:
    % psat = 611.15 exp[(23.036 - Tc/333.7)*Tc/(279.82 + Tc)]
    %
    
    p.water.Buck.pos.A  = 611.21;
    p.water.Buck.pos.B1 = 18.678;
    p.water.Buck.pos.B2 = 234.5;
    p.water.Buck.pos.B3 = 257.14;
    
    p.water.Buck.neg.A  = 611.15;
    p.water.Buck.neg.B1 = 23.036;
    p.water.Buck.neg.B2 = 333.7;
    p.water.Buck.neg.B3 = 279.82;
    
end
