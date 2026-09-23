function [Vcell, out] = calc_cell_voltage( ...
    T, pH2, pO2, j, lambda, DO2eff, cO2_cCL, Ldiff, p)
    
    %   PEMFC cell voltage and loss components
    %
    %   [Vcell,out] = calc_cell_voltage( ...
    %       T,pH2,pO2,j,lambda,DO2eff,cO2_cCL,Ldiff,p)
    %
    %   All quantities use SI units.
    %
    %   INPUT
    %   -----
    %   T          : representative cell/electrode temperature [K]
    %   pH2        : H2 partial pressure at aCL [Pa]
    %   pO2        : O2 partial pressure at cCL [Pa]
    %   j          : current density [A/m^2]
    %   lambda     : representative membrane water content [-]
    %   DO2eff     : effective O2 diffusivity [m^2/s]
    %   cO2_cCL    : O2 concentration at cCL [mol/m^3]
    %   Ldiff      : effective cathode diffusion length [m]
    %
    %   OUTPUT
    %   ------
    %   Vcell
    %
    %   out.Erev
    %   out.j0
    %   out.etaAct
    %   out.kappaPem
    %   out.etaOhm
    %   out.jLim
    %   out.etaCon
    %   out.transportValid
    
    R = p.const.R;
    F = p.const.F;
    
    %% ------------------------------------------------------------------------
    % 1. Reversible voltage
    % -------------------------------------------------------------------------
    
    if pH2 <= 0 || pO2 <= 0
        error('Gas partial pressures must be positive.');
    end
    
    Erev = ...
        p.echem.Erev_ref + ...
        p.echem.dEdT .* (T - p.const.Tref) + ...
        R*T/(2*F) .* ...
        log( ...
            (pH2/p.const.p0) .* ...
            sqrt(pO2/p.const.p0) ...
            );
    
    
    %% ------------------------------------------------------------------------
    % 2. Exchange current density
    % -------------------------------------------------------------------------
    
    j0 = ...
        p.echem.j0_ref .* ...
        exp( ...
            -(p.echem.Ea/R) .* ...
            (1/T - 1/p.const.Tref) ...
            );
    
    
    %% ------------------------------------------------------------------------
    % 3. Activation loss
    % -------------------------------------------------------------------------
    
    etaAct = ...
        R*T/(p.echem.alpha*F) .* ...
        asinh( ...
            j ./ (2*j0) ...
            );
    
    
    %% ------------------------------------------------------------------------
    % 4. Ohmic loss
    % -------------------------------------------------------------------------
    
    kappaPem = membrane_conductivity(lambda,T,p);
    
    if kappaPem <= 0
        etaOhm = Inf;
    else
        etaOhm = ...
            j .* ...
            ( ...
            p.geom.Lpem ./ kappaPem ...
            + p.echem.Rc ...
            );
    end
    
    
    %% ------------------------------------------------------------------------
    % 5. Limiting current density
    % -------------------------------------------------------------------------
    
    if Ldiff <= 0
        error('Ldiff must be positive.');
    end
    
    jLim = ...
        4 * F .* ...
        DO2eff .* ...
        cO2_cCL ./ ...
        Ldiff;
    
    
    %% ------------------------------------------------------------------------
    % 6. Concentration loss
    % -------------------------------------------------------------------------
    
    transportValid = ...
        jLim > 0 && ...
        j >= 0 && ...
        j < jLim;
    
    if j == 0
    
        etaCon = 0;
    
    elseif transportValid
    
        etaCon = ...
            -R*T/(4*F) .* ...
            log( ...
                1 - j/jLim ...
                );
    
    else
    
        % Physical transport limit has been exceeded.
        % Do NOT hide this using arbitrary clipping at this level.
        etaCon = Inf;
    
    end
    
    
    %% ------------------------------------------------------------------------
    % 7. Cell voltage
    % -------------------------------------------------------------------------
    
    Vcell = ...
        Erev ...
        - etaAct ...
        - etaOhm ...
        - etaCon;
    
    
    %% ------------------------------------------------------------------------
    % Diagnostics
    % -------------------------------------------------------------------------
    
    out.Erev = Erev;
    
    out.j0 = j0;
    
    out.etaAct = etaAct;
    
    out.kappaPem = kappaPem;
    
    out.etaOhm = etaOhm;
    
    out.jLim = jLim;
    
    out.etaCon = etaCon;
    
    out.transportValid = transportValid;

end