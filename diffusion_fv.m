function [divDiff, Jface, info] = ...
    diffusion_fv(phi, D, g, idx, bc)

%DIFFUSION_FV
% 1-D finite-volume diffusion operator on a contiguous sub-domain
%
%   [divDiff,Jface,info] =
%       diffusion_fv(phi,D,g,idx,bc)
%
%
% Governing flux:
%
%       J = -D * d(phi)/dx
%
% Diffusion term:
%
%       divDiff =
%       d/dx (D*d(phi)/dx)
%       =
%       -dJ/dx
%
%
% INPUT
% -----
% phi :
%     Cell-centered transported quantity.
%
%     Examples:
%       gas concentration    [mol/m^3]
%       water concentration  [kg/m^3]
%
% D :
%     Cell-centered effective diffusivity [m^2/s]
%
% g :
%     Grid structure from build_grid.m
%
% idx :
%     Global cell indices corresponding to this transport domain.
%     Must form one contiguous region.
%
% bc.left.type / bc.right.type:
%
%     'dirichlet'
%         prescribed phi at boundary face
%
%     'neumann'
%         prescribed physical flux J in +x direction
%
%     'noflux'
%         J = 0
%
% bc.left.value / bc.right.value:
%     boundary value for dirichlet/neumann condition
%
%
% OUTPUT
% ------
% divDiff :
%
%     finite-volume approximation of
%
%       d/dx(D*dphi/dx)
%
%     same length as phi
%
% Jface :
%     diffusive flux on all faces
%
%       length = Ncell + 1
%
%     positive J means transport in +x direction
%
% info :
%     diagnostic quantities


%% ------------------------------------------------------------------------
% 1. Prepare input
% -------------------------------------------------------------------------

phi = phi(:);
D   = D(:);
idx = idx(:);

N = numel(idx);

if numel(phi) ~= N
    error('diffusion_fv:SizeMismatch', ...
        'phi and idx must have the same number of elements.');
end

if numel(D) ~= N
    error('diffusion_fv:SizeMismatch', ...
        'D and idx must have the same number of elements.');
end

if any(D < 0) || any(~isfinite(D))
    error('diffusion_fv:InvalidDiffusivity', ...
        'Diffusivity must be finite and non-negative.');
end

if any(~isfinite(phi))
    error('diffusion_fv:InvalidPhi', ...
        'phi contains non-finite values.');
end


%% ------------------------------------------------------------------------
% 2. Domain must be contiguous
% -------------------------------------------------------------------------

if N > 1 && any(diff(idx) ~= 1)
    error('diffusion_fv:NonContiguousDomain', ...
        'idx must describe one contiguous spatial region.');
end

if idx(1) < 1 || idx(end) > g.N
    error('diffusion_fv:InvalidIndex', ...
        'idx lies outside the global grid.');
end


%% ------------------------------------------------------------------------
% 3. Extract local geometry
% -------------------------------------------------------------------------

x  = g.x(idx);
dx = g.dx(idx);

% If local cells are global cells i0 ... i1,
% corresponding faces are i0 ... i1+1.

i0 = idx(1);
i1 = idx(end);

xf = g.xf(i0:i1+1);


%% ------------------------------------------------------------------------
% 4. Allocate face fluxes
% -------------------------------------------------------------------------

Jface = zeros(N+1,1);


%% ------------------------------------------------------------------------
% 5. Internal faces
% -------------------------------------------------------------------------
%
% Face f lies between local cells f-1 and f:
%
%       cell L | face | cell R
%
% Resistance formulation:
%
%       Rdiff = dL/DL + dR/DR
%
%       J =
%       -(phiR - phiL)/Rdiff
%

for f = 2:N

    iL = f-1;
    iR = f;

    xFace = xf(f);

    dL = xFace - x(iL);
    dR = x(iR) - xFace;

    DL = D(iL);
    DR = D(iR);

    if DL == 0 || DR == 0

        % Zero diffusivity on either side blocks transport
        Jface(f) = 0;

    else

        Rdiff = ...
            dL / DL + ...
            dR / DR;

        Jface(f) = ...
            -(phi(iR) - phi(iL)) / Rdiff;

    end

end


%% ------------------------------------------------------------------------
% 6. Left boundary
% -------------------------------------------------------------------------

Jface(1) = boundary_flux( ...
    'left', ...
    phi(1), ...
    D(1), ...
    x(1), ...
    xf(1), ...
    bc.left);


%% ------------------------------------------------------------------------
% 7. Right boundary
% -------------------------------------------------------------------------

Jface(end) = boundary_flux( ...
    'right', ...
    phi(end), ...
    D(end), ...
    x(end), ...
    xf(end), ...
    bc.right);


%% ------------------------------------------------------------------------
% 8. Flux divergence
% -------------------------------------------------------------------------
%
% PDE:
%
%       divDiff =
%       -dJ/dx
%
% Finite volume:
%
%       divDiff_i =
%       (J_left - J_right)/dx_i
%

divDiff = ...
    (Jface(1:end-1) - Jface(2:end)) ./ dx;


%% ------------------------------------------------------------------------
% 9. Conservation diagnostics
% -------------------------------------------------------------------------
%
% For unit cross-sectional area:
%
% integral(divDiff dx)
% =
% J_left - J_right
%

integralDiv = ...
    sum(divDiff .* dx);

boundaryNet = ...
    Jface(1) - Jface(end);

info.balanceError = ...
    integralDiv - boundaryNet;

info.Jleft  = Jface(1);
info.Jright = Jface(end);

info.integralDiv = integralDiv;

info.conserved = ...
    abs(info.balanceError) < ...
    max(1e-12, ...
        1e-10*max(abs(boundaryNet),1));


end


%% =========================================================================
% Local helper: boundary flux
% =========================================================================

function J = boundary_flux(side, phiCell, Dcell, ...
                           xCell, xFace, bc)

if ~isfield(bc,'type')
    error('Boundary condition must contain field "type".');
end

type = lower(string(bc.type));

switch type

    %% --------------------------------------------------------------------
    % No flux
    % ---------------------------------------------------------------------

    case "noflux"

        J = 0;


    %% --------------------------------------------------------------------
    % Prescribed physical flux
    %
    % IMPORTANT:
    %
    % bc.value is always defined as flux in the +x direction.
    % ---------------------------------------------------------------------

    case "neumann"

        if ~isfield(bc,'value')
            error('Neumann boundary requires bc.value.');
        end

        J = bc.value;


    %% --------------------------------------------------------------------
    % Prescribed phi at the boundary face
    % ---------------------------------------------------------------------

    case "dirichlet"

        if ~isfield(bc,'value')
            error('Dirichlet boundary requires bc.value.');
        end

        phiBoundary = bc.value;

        if Dcell == 0
            J = 0;
            return;
        end

        switch lower(side)

            case 'left'

                d = xCell - xFace;

                % derivative from boundary -> cell:
                %
                % dphi/dx =
                % (phiCell - phiBoundary)/d

                J = ...
                    -Dcell .* ...
                    (phiCell - phiBoundary) ./ d;


            case 'right'

                d = xFace - xCell;

                % derivative from cell -> boundary:
                %
                % dphi/dx =
                % (phiBoundary - phiCell)/d

                J = ...
                    -Dcell .* ...
                    (phiBoundary - phiCell) ./ d;


            otherwise

                error('Unknown boundary side.');

        end


    otherwise

        error('diffusion_fv:UnknownBC', ...
            'Unknown boundary type "%s".',bc.type);

end

end