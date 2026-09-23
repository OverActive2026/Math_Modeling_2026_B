function [qCond, qFace, info] = ...
    heat_conduction_fv(T, k, Tamb, h, g)
 
    % 1-D finite-volume heat conduction with convective boundaries.
    %
    %   [qCond,qFace,info] =
    %       heat_conduction_fv(T,k,Tamb,h,g)
    %
    %
    % Governing equation contribution:
    %
    %       qCond =
    %       d/dx ( k*dT/dx )
    %
    % Define conductive heat flux in +x direction:
    %
    %       q_x = -k*dT/dx
    %
    % Therefore:
    %
    %       qCond = -d(q_x)/dx
    %
    %
    % Boundary condition:
    %
    %       -k*dT/dn = h*(T_surface - Tamb)
    %
    %
    % INPUT
    % -----
    % T :
    %     cell-centered temperature [K]
    %     length = g.N
    %
    % k :
    %     cell-centered thermal conductivity [W/(m K)]
    %     scalar or length g.N
    %
    % Tamb :
    %     ambient temperature [K]
    %
    % h :
    %     convective heat-transfer coefficient [W/(m^2 K)]
    %
    % g :
    %     grid structure
    %
    %
    % OUTPUT
    % ------
    % qCond :
    %     volumetric heat-conduction term [W/m^3]
    %
    % qFace :
    %     conductive heat flux on faces [W/m^2]
    %
    %     positive value = heat flow in +x direction
    %
    % info :
    %     diagnostic quantities
    
    
    %% ------------------------------------------------------------------------
    % 1. Input preparation
    % -------------------------------------------------------------------------
    
    T = T(:);
    
    if numel(T) ~= g.N
        error('heat_conduction_fv:SizeMismatch', ...
            'T must contain exactly g.N values.');
    end
    
    if isscalar(k)
    
        k = repmat(k,g.N,1);
    
    else
    
        k = k(:);
    
    end
    
    if numel(k) ~= g.N
        error('heat_conduction_fv:SizeMismatch', ...
            'k must be scalar or contain exactly g.N values.');
    end
    
    
    %% ------------------------------------------------------------------------
    % 2. Validation
    % -------------------------------------------------------------------------
    
    if any(~isfinite(T)) || any(T <= 0)
        error('Temperature must contain finite values in Kelvin.');
    end
    
    if any(~isfinite(k)) || any(k <= 0)
        error('Thermal conductivity must be finite and positive.');
    end
    
    if ~isscalar(Tamb) || ~isfinite(Tamb) || Tamb <= 0
        error('Tamb must be a finite scalar temperature in Kelvin.');
    end
    
    if ~isscalar(h) || ~isfinite(h) || h <= 0
        error('h must be a finite positive scalar.');
    end
    
    
    %% ------------------------------------------------------------------------
    % 3. Allocate face heat flux
    % -------------------------------------------------------------------------
    
    N = g.N;
    
    qFace = zeros(N+1,1);
    
    
    %% ------------------------------------------------------------------------
    % 4. Internal conductive faces
    % -------------------------------------------------------------------------
    %
    % For internal face:
    %
    %       cell L | face | cell R
    %
    % Heat flux:
    %
    %       q =
    %       -(TR - TL) /
    %       (dL/kL + dR/kR)
    %
    % This is equivalent to two thermal resistances in series.
    %
    
    for f = 2:N
    
        iL = f-1;
        iR = f;
    
        xFace = g.xf(f);
    
        dL = ...
            xFace - g.x(iL);
    
        dR = ...
            g.x(iR) - xFace;
    
        Rth = ...
            dL ./ k(iL) + ...
            dR ./ k(iR);
    
        qFace(f) = ...
            -(T(iR)-T(iL)) ./ Rth;
    
    end
    
    
    %% ------------------------------------------------------------------------
    % 5. Left convective boundary
    % -------------------------------------------------------------------------
    %
    % +x points from anode -> cathode.
    %
    % Heat leaving through the LEFT surface travels in -x direction.
    %
    % Series resistance from cell center to ambient:
    %
    %       R'' = d/k + 1/h
    %
    % Outward heat loss:
    %
    %       q_out =
    %       (Tcell - Tamb)/R''
    %
    % Therefore heat flux in +x direction:
    %
    %       q_x,left = -q_out
    %
    
    dLeft = ...
        g.x(1) - g.xf(1);
    
    Rleft = ...
        dLeft ./ k(1) + ...
        1 ./ h;
    
    qOutLeft = ...
        (T(1)-Tamb) ./ Rleft;
    
    qFace(1) = ...
        -qOutLeft;
    
    
    %% ------------------------------------------------------------------------
    % 6. Right convective boundary
    % -------------------------------------------------------------------------
    %
    % Heat leaving through the RIGHT surface travels in +x direction.
    %
    
    dRight = ...
        g.xf(end) - g.x(end);
    
    Rright = ...
        dRight ./ k(end) + ...
        1 ./ h;
    
    qOutRight = ...
        (T(end)-Tamb) ./ Rright;
    
    qFace(end) = ...
        qOutRight;
    
    
    %% ------------------------------------------------------------------------
    % 7. Finite-volume heat-flux divergence
    % -------------------------------------------------------------------------
    %
    %       qCond =
    %       -dq/dx
    %
    %       =
    %       (q_left - q_right)/dx
    %
    
    [qCond, divInfo] = ...
        flux_divergence_fv( ...
            qFace, ...
            g, ...
            (1:g.N).' ...
            );
    
    
    %% ------------------------------------------------------------------------
    % 8. Diagnostics
    % -------------------------------------------------------------------------
    
    info.qOutLeft  = qOutLeft;
    info.qOutRight = qOutRight;
    
    info.totalHeatLossPerArea = ...
        qOutLeft + qOutRight;
    
    info.leftSurfaceResistance = ...
        Rleft;
    
    info.rightSurfaceResistance = ...
        Rright;
    
    info.divergenceInfo = ...
        divInfo;
    
    
    %% ------------------------------------------------------------------------
    % 9. Global energy check
    % -------------------------------------------------------------------------
    %
    % For unit cross-sectional area:
    %
    %   integral(qCond dx)
    %
    % should equal
    %
    %   -(qOutLeft + qOutRight)
    %
    
    info.integratedConduction = ...
        sum(qCond .* g.dx);
    
    info.expectedIntegratedConduction = ...
        -(qOutLeft + qOutRight);
    
    info.energyBalanceError = ...
        info.integratedConduction - ...
        info.expectedIntegratedConduction;
    
    tol = ...
        max(1e-10, ...
        1e-9 * max(abs(info.totalHeatLossPerArea),1));
    
    info.energyBalanced = ...
        abs(info.energyBalanceError) <= tol;

end