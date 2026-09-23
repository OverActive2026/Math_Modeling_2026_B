function p = params()
%PARAMS  Parameters for the 1-D PEMFC cold-start model
%
%   p = params()
%
%   Parameter policy
%   ----------------
%   p.const   : physical constants
%   p.geom    : geometry
%   p.trans   : gas / water transport parameters
%   p.echem   : electrochemical parameters
%   p.thermal : thermal parameters
%   p.water   : water / membrane parameters
%   p.ice     : cold-start ice model parameters
%   p.mat     : layer material properties
%   p.bc      : boundary-condition parameters
%
%   Values explicitly given in the problem statement and Attachment 1 are
%   filled in. NaN is reserved for parameters that still require model
%   selection or calibration; these are listed in p.meta.toBeCalibrated.
%
%   All internal calculations use SI units.

    %% ========================================================================
    %  1. Universal physical constants
    % ========================================================================
    
    p.const.R = 8.314;          % Universal gas constant [J/(mol K)]
    p.const.F = 96485;          % Faraday constant [C/mol]
    
    p.const.Mw = 0.018;         % Molar mass of water [kg/mol]
    p.const.MH2 = 2.016e-3;     % Molar mass of hydrogen [kg/mol]
    p.const.MO2 = 31.998e-3;    % Molar mass of oxygen [kg/mol]
    p.const.MN2 = 28.014e-3;    % Molar mass of nitrogen [kg/mol]
    
    p.const.p0   = 101325;      % Reference pressure [Pa]
    p.const.Tref = 298.15;      % Reference temperature [K]
    p.const.Tf   = 273.15;      % Freezing temperature of water [K]
    
    
    %% ========================================================================
    %  2. Cell geometry
    % ========================================================================
    %
    % Through-plane sequence:
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
    % % Additional stack geometry from Attachment 1.
    % p.geom.LendPlate = 10e-3;   % End-plate thickness [m]
    % p.geom.LaBP      = 2e-3;    % Anode bipolar-plate thickness [m]
    % p.geom.LcBP      = 2e-3;    % Cathode bipolar-plate thickness [m]

    % % Thermal thickness when both bipolar plates are represented explicitly.
    % p.geom.Lthermal = p.geom.LaBP + p.geom.Lmea + p.geom.LcBP;

    % Cell active area from Attachment 1.
    p.geom.Acell = 25e-4;       % 25 cm^2 -> [m^2]
    
    
    %% ========================================================================
    %  3. Gas transport
    % ========================================================================
    %
    % Problem statement:
    %
    % D_k,eff =
    %   D_k,ref * (T/298.15)^1.75 * (101325/p) * eps_g^1.5
    %
    
    p.trans.DH2_ref = 1.10e-4;  % H2 reference diffusivity [m^2/s]
    p.trans.DO2_ref = 2.20e-5;  % O2 reference diffusivity [m^2/s]
    
    p.trans.Texp = 1.75;        % temperature exponent
    p.trans.epsExp = 1.5;       % porosity exponent
    
    % Operating pressure from Attachment 1.
    % 这里实际上附件 1 只给出一个统一的压力值101325，未分别给出阳极、阴极
    p.bc.pAnode   = 101325;      % [Pa]
    p.bc.pCathode = 101325;      % [Pa]

    % Attachment 1 supplies mass fractions. Store them explicitly, then
    % convert to mole fractions for c = y*p/(R*T). Do not use wO2=0.233 as
    % an oxygen mole fraction.
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
    %  4. Water transport
    % ========================================================================
    
    % Reference diffusivity of water in porous media
    p.water.Dw_ref_anode   = 8.69e-5;  % aGDL / aCL [m^2/s]
    p.water.Dw_ref_cathode = 2.48e-5;  % cGDL / cCL [m^2/s]
    
    % Electro-osmotic drag:
    %
    %       nd = 2.5 * lambda / 22
    %
    p.water.ndCoeff = 2.5 / 22;
    
    % Membrane-water diffusion relation is explicitly specified by the
    % problem and will be implemented in calc_Dwater.m:
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
    
    % Membrane equivalent weight and density from Attachment 1:
    % required by
    %
    % lambda = EW * mw / (rho_pem * Mw)
    %
    % 1000 g/mol = 1 kg/mol.
    p.water.EW      = 1.0;       % [kg/mol]
    p.water.rho_pem = 2150;      % [kg/m^3]
    p.water.lambda0 = 3.0;       % Initial membrane water content [-]
    p.water.clIonomerVolumeFraction = 0.3;

    % Initial membrane water mass concentration derived from Eq. (21).
    p.water.mwMembrane0 = p.water.lambda0 * p.water.rho_pem * ...
        p.const.Mw / p.water.EW; % [kg/m^3]

    % Porous-media inputs retained for the later liquid-water model.
    p.porous.contactAngleGDL_deg = 110;
    p.porous.contactAngleCL_deg  = 100;
    p.porous.contactAngleGDL = 110*pi/180; % [rad]
    p.porous.contactAngleCL  = 100*pi/180; % [rad]
    p.porous.permeabilityGDL = 6.2e-12;    % [m^2]
    p.porous.permeabilityCL  = 6.2e-13;    % [m^2]

    % Phase-change / transfer coefficients supplied as dimensionless model
    % coefficients. Directional ordering is not specified in Attachment 1,
    % so paired values are kept together rather than assigned by guesswork.
    p.water.phase.memVaporCoeff  = [0.001, 1.0];
    p.water.phase.memLiquidCoeff = 0.5;
    p.water.phase.memIceCoeff    = 1.0;
    p.water.phase.vaporLiquidCoeff = [1.0, 1.0];
    p.water.phase.vaporIceCoeff  = 1.0e-4;
    
    
    %% ========================================================================
    %  5. Electrochemistry
    % ========================================================================
    
    p.echem.alpha = 0.5;          % charge-transfer coefficient

    p.echem.Voc0 = 0.95;          % Initial open-circuit voltage [V]
    
    % Initial calibration value specified in the statement
    p.echem.j0_ref = 0.01;        % [A/m^2]
    
    p.echem.Ea = 67000;           % activation energy [J/mol]
    
    % Thermoneutral voltage
    p.echem.Eth = 1.48;           % [V]
    
    % Reversible-voltage constants:
    %
    % Erev =
    % 1.229 - 8.5e-4*(T-298.15)
    % + RT/(2F) ln[ (pH2/p0)*(pO2/p0)^0.5 ]
    %
    p.echem.Erev_ref = 1.229;     % [V]
    p.echem.dEdT     = -8.5e-4;   % [V/K]
    
    % Area-specific contact resistance
    %
    % Problem statement:
    %       Rc = 0.01 ohm cm^2
    %
    % Convert:
    %       1 cm^2 = 1e-4 m^2
    %
    p.echem.Rc = 0.01e-4;         % [ohm m^2]
    
    % Ice coverage exponent and stack concentration-polarization factors
    % from Attachment 1. They are inactive in the no-ice single-cell MVP.
    p.echem.gammaIceArea = 3.5;
    p.echem.concentrationFactorEnd = 10;
    p.echem.concentrationFactorMiddle = 1;

    % 材料不存在这个公式，Ldiff这个变量不知道是什么
    % % Equivalent cathode diffusion distance. For the cell-centered cCL
    % % concentration, use the distance from the cathode boundary to the cCL
    % % representative center: cGDL thickness + half cCL thickness.
    % p.echem.Ldiff = p.geom.LcGDL + 0.5*p.geom.LcCL; % [m]
    
    
    %% ========================================================================
    %  6. Thermal parameters
    % ========================================================================
    
    % Convective heat-transfer coefficient
    p.thermal.h = 40;             % [W/(m^2 K)]
    
    % Material properties from Attachment 1.
    %
    % order:
    %       [aGDL, aCL, PEM, cCL, cGDL]
    %
    % rho : density         [kg/m^3]
    % cp  : heat capacity   [J/(kg K)]
    % k   : conductivity    [W/(m K)]
    
    p.mat.rho = [185, 970, 2150, 970, 185];
    p.mat.cp  = [545, 240, 1050, 240, 545];
    p.mat.k   = [0.3, 0.27, 0.24, 0.27, 0.3];

    % Properties not mapped directly onto the five-layer transport grid.
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
    %  7. Dry porosity of porous layers
    % ========================================================================
    %
    % eps0 is needed by:
    %
    %       eps_g = eps0 - eps_l - eps_ice
    %
    % and hence by gas/water effective diffusivities.
    %
    % PEM is non-porous for this simplified treatment.
    %
    % 附件 1 没有给膜的孔隙率，这里给了0
    p.mat.eps0 = [ ...
        0.8,    ...     % aGDL
        0.3916, ...     % aCL
        0,   ...       % PEM
        0.4207, ...     % cCL
        0.8 ];          % cGDL
    
    
    %% ========================================================================
    %  8. Ice / cold-start model
    % ========================================================================
    %
    % These are NOT prescribed by the baseline model.
    % They belong to our cold-start extension and should eventually be
    % calibrated against the -20 C / -25 C experimental data.
    %
    
    % Densities, latent heat and active-area exponent from Attachment 1.
    p.ice.rho_liquid = p.phase.liquid.rho; % [kg/m^3]
    p.ice.rho_ice    = p.phase.ice.rho;    % [kg/m^3]
    p.ice.Lfusion    = p.thermal.Lfusion;  % [J/kg]
    
    % Initial proposed kinetic law:
    %
    % rFreeze = kFreeze * ml * max(Tf - T,0)
    % rMelt   = kMelt   * mi * max(T - Tf,0)
    %
    % These are model parameters, NOT given by the problem.
    
    p.ice.kFreeze = NaN;          % calibration parameter
    p.ice.kMelt   = NaN;          % calibration parameter
    
    % Ice-blockage effect on electrochemically active area:
    %
    %       a_eff = (1 - sIce)^gammaArea
    %
    p.ice.gammaArea = p.echem.gammaIceArea;
    
    % Numerical lower bound for remaining gas porosity.
    % This is a numerical safeguard, not a physical parameter.
    p.num.epsMin = 1e-8;
    p.num.residualTol = 1e-6;

    % Initial and ambient conditions from Attachment 1.
    p.ic.T0 = 253.15;             % [K]
    p.bc.Tamb = 253.15;           % [K]
    p.ic.liquidWater = 0.0;
    p.ic.iceVolumeFraction = 0.0;
    
    
    %% ========================================================================
    %  9. Membrane conductivity correlation
    % ========================================================================
    %
    % kappa_pem =
    % (0.5139 lambda - 0.326) *
    % exp[1268(1/303.15 - 1/T)]
    %
    % No fitted parameter is required, but coefficients are stored here
    % instead of hard-coded in the voltage function.
    
    p.water.kappa_a0 = 0.5139;
    p.water.kappa_a1 = -0.326;
    p.water.kappa_E  = 1268;
    
    
    %% ========================================================================
    %  10. Buck saturation-pressure correlation
    % ========================================================================
    %
    % Tc = T - 273.15 [degC]
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
    
    
    %% ========================================================================
    %  11. Current constraints used later in Q2
    % ========================================================================
    %
    % Keep these here although Q1 does not need the limits.
    
    p.limit.jMax_Acm2 = 0.5;          % [A/cm^2]
    p.limit.jMax      = 0.5 * 1e4;    % [A/m^2]
    
    p.limit.qChargeMax_Acm2 = 20;     % [C/cm^2]
    p.limit.qChargeMax      = 20*1e4; % [C/m^2]
    
    
    %% ========================================================================
    %  12. Auxiliary-heater constraint used later in Q3 / Q4
    % ========================================================================
    
    p.limit.qAuxMax_Wcm2 = 1;         % [W/cm^2]
    p.limit.qAuxMax      = 1e4;       % [W/m^2]
    
    
    %% ========================================================================
    %  13. Parameter status bookkeeping
    % ========================================================================
    %
    % Useful during development: immediately know what has not yet been
    % populated from Attachment 1 / literature.
    
    p.meta.missingFromAttachment1 = {};

    p.meta.toBeCalibrated = { ...
        'echem.j0_ref', ...
        'ice.kFreeze', ...
        'ice.kMelt' ...
        };

    p.meta.derivedAssumptions = { ...
        'bc.yO2 and bc.yN2 are converted from attachment mass fractions', ...
        'echem.Ldiff = LcGDL + 0.5*LcCL for a cell-centered cCL value' ...
        };
    
end
