function [src, info] = calc_reaction_sources(j, g, p)
%CALC_REACTION_SOURCES
% Electrochemical reaction source terms for the 1-D PEMFC model.
%
%   [src,info] = calc_reaction_sources(j,g,p)
%
%
% INPUT
% -----
% j :
%     Applied current density [A/m^2]
%
% g :
%     Grid structure from build_grid.m
%
% p :
%     Parameter structure from pemfc_params.m
%
%
% OUTPUT
% ------
% src.SH2 :
%     H2 molar source term [mol/(m^3 s)]
%     Negative = consumption
%
% src.SO2 :
%     O2 molar source term [mol/(m^3 s)]
%     Negative = consumption
%
% src.Sw :
%     Water mass source term [kg/(m^3 s)]
%     Positive = generation
%
%
% Reactions:
%
%   H2  -> 2H+ + 2e-
%
%   O2 + 4H+ + 4e- -> 2H2O
%
%
% Under the problem's uniform-reaction assumption:
%
%   aCL:
%
%       SH2 = -j/(2*F*LaCL)
%
%   cCL:
%
%       SO2 = -j/(4*F*LcCL)
%
%       Sw  = Mw*j/(2*F*LcCL)


%% ------------------------------------------------------------------------
% 1. Input validation
% -------------------------------------------------------------------------

if ~isscalar(j) || ~isfinite(j)
    error('calc_reaction_sources:InvalidCurrent', ...
        'j must be a finite scalar.');
end

if j < 0
    error('calc_reaction_sources:NegativeCurrent', ...
        ['Current implementation assumes PEM fuel-cell operation ', ...
         'with j >= 0.']);
end


F  = p.const.F;
Mw = p.const.Mw;

LaCL = p.geom.LaCL;
LcCL = p.geom.LcCL;


if LaCL <= 0 || LcCL <= 0
    error('Catalyst-layer thickness must be positive.');
end


%% ------------------------------------------------------------------------
% 2. Allocate global source arrays
% -------------------------------------------------------------------------

src.SH2 = zeros(g.N,1);
src.SO2 = zeros(g.N,1);
src.Sw  = zeros(g.N,1);


%% ------------------------------------------------------------------------
% 3. Anode H2 consumption
% -------------------------------------------------------------------------

SH2_aCL = ...
    -j ./ ...
    (2 * F * LaCL);

src.SH2(g.idx_aCL) = SH2_aCL;


%% ------------------------------------------------------------------------
% 4. Cathode O2 consumption
% -------------------------------------------------------------------------

SO2_cCL = ...
    -j ./ ...
    (4 * F * LcCL);

src.SO2(g.idx_cCL) = SO2_cCL;


%% ------------------------------------------------------------------------
% 5. Cathode water generation
% -------------------------------------------------------------------------

Sw_cCL = ...
    Mw .* j ./ ...
    (2 * F * LcCL);

src.Sw(g.idx_cCL) = Sw_cCL;


%% ------------------------------------------------------------------------
% 6. Diagnostics
% -------------------------------------------------------------------------

info.SH2_aCL = SH2_aCL;
info.SO2_cCL = SO2_cCL;
info.Sw_cCL  = Sw_cCL;


% -------------------------------------------------------------------------
% Integrated reaction rates per unit active area
%
% Because the source is volumetric:
%
%   integral S dx
%
% gives the source per unit electrode area.
% -------------------------------------------------------------------------

info.H2Consumption_Area = ...
    -sum(src.SH2 .* g.dx);         % [mol/(m^2 s)]

info.O2Consumption_Area = ...
    -sum(src.SO2 .* g.dx);         % [mol/(m^2 s)]

info.WaterGeneration_Area = ...
     sum(src.Sw .* g.dx);          % [kg/(m^2 s)]


%% ------------------------------------------------------------------------
% 7. Theoretical Faraday-law values
% -------------------------------------------------------------------------

info.H2Expected = ...
    j / (2*F);

info.O2Expected = ...
    j / (4*F);

info.WaterExpected = ...
    Mw*j / (2*F);


%% ------------------------------------------------------------------------
% 8. Conservation / consistency errors
% -------------------------------------------------------------------------

info.errorH2 = ...
    info.H2Consumption_Area - ...
    info.H2Expected;

info.errorO2 = ...
    info.O2Consumption_Area - ...
    info.O2Expected;

info.errorWater = ...
    info.WaterGeneration_Area - ...
    info.WaterExpected;


scale = max(abs(j/F),1);

tol = max(1e-14,1e-10*scale);

info.consistent = ...
    abs(info.errorH2) < tol && ...
    abs(info.errorO2) < tol && ...
    abs(info.errorWater) < tol;

end