function g = build_grid(p,nCellLayer)
%BUILD_GRID  Build 1-D through-plane mesh for PEMFC cold-start model
%
%   g = build_grid(p)
%   g = build_grid(p,nCellLayer)
%
%   The through-plane domain is:
%
%       aGDL | aCL | PEM | cCL | cGDL
%
%   A cell-centered finite-volume mesh is used.
%   All internal calculations use SI units.
%
%   INPUT
%   -----
%   p : parameter structure returned by params()
%
%   nCellLayer : 1x5 integer vector
%       Number of control volumes in
%       [aGDL, aCL, PEM, cCL, cGDL]
%
%       Default:
%           [4 2 3 3 4]
%
%   OUTPUT
%   ------
%   g : structure containing
%       .N                 total number of cells
%       .x                 cell-center coordinates [m]
%       .xf                cell-face coordinates   [m]
%       .dx                cell width              [m]
%       .dxc               center-to-center distance [m]
%       .layerId           material ID of each cell
%
%       .idx_aGDL
%       .idx_aCL
%       .idx_PEM
%       .idx_cCL
%       .idx_cGDL
%
%       .idx_H2            cells containing H2 transport
%       .idx_O2            cells containing O2 transport
%       .idx_porous        porous-media cells
%
%       .layerThickness
%       .layerBoundary
%       .interfaceFace
%
%   Patrick / 2026 Math Modeling Contest B

%% ------------------------------------------------------------------------
%  1. Layer definitions
% -------------------------------------------------------------------------

    g.layerNames = {'aGDL','aCL','PEM','cCL','cGDL'};
    
    if nargin < 1 || isempty(p) || ~isstruct(p) || ...
            ~isfield(p,'geom') || ~isfield(p.geom,'layerThickness')
        error('build_grid requires the parameter structure returned by params().');
    end

    % Use params.m as the single source of truth for geometry.
    g.layerThickness = p.geom.layerThickness(:).';

    if numel(g.layerThickness) ~= 5 || ...
            any(~isfinite(g.layerThickness)) || ...
            any(g.layerThickness <= 0)
        error('p.geom.layerThickness must contain five positive finite values.');
    end
    
    g.nLayer = numel(g.layerThickness);
    
    
    %% ------------------------------------------------------------------------
    %  2. Number of cells in each layer
    % -------------------------------------------------------------------------
    
    if nargin < 2 || isempty(nCellLayer)
        % Coarse first-version mesh:
        % GDL relatively coarse, CL / PEM relatively fine
        nCellLayer = [4 2 3 3 4];
    end
    
    if numel(nCellLayer) ~= g.nLayer
        error('nCellLayer must contain exactly 5 elements.');
    end
    
    if any(~isfinite(nCellLayer))
        error('All entries of nCellLayer must be finite.');
    end

    if any(nCellLayer <= 0) || any(mod(nCellLayer,1) ~= 0)
        error('All entries of nCellLayer must be positive integers.');
    end
    
    nCellLayer = nCellLayer(:).';
    
    g.nCellLayer = nCellLayer;
    g.N = sum(nCellLayer);
    
    
    %% ------------------------------------------------------------------------
    %  3. Layer boundaries
    % -------------------------------------------------------------------------
    %
    %  x = 0 is the outer boundary of the anode GDL.
    %
    %  layerBoundary:
    %
    %     0 | aGDL | aCL | PEM | cCL | cGDL |
    %
    
    g.layerBoundary = [0, cumsum(g.layerThickness)];
    
    g.Ltotal = g.layerBoundary(end);
    
    
    %% ------------------------------------------------------------------------
    %  4. Build finite-volume cell faces
    % -------------------------------------------------------------------------
    
    xf = zeros(g.N + 1,1);
    
    faceCounter = 1;
    
    for k = 1:g.nLayer
    
        xL = g.layerBoundary(k);
        xR = g.layerBoundary(k+1);
    
        nk = nCellLayer(k);
    
        % Local faces for this material layer
        xf_local = linspace(xL,xR,nk+1).';
    
        if k == 1
            xf(1:nk+1) = xf_local;
            faceCounter = nk + 1;
        else
            % First face already exists as the previous layer interface
            xf(faceCounter+1:faceCounter+nk) = xf_local(2:end);
            faceCounter = faceCounter + nk;
        end
    end
    
    g.xf = xf;
    
    
    %% ------------------------------------------------------------------------
    %  5. Cell centers and cell widths
    % -------------------------------------------------------------------------
    
    g.dx = diff(g.xf);
    
    g.x = 0.5 .* ( ...
        g.xf(1:end-1) + ...
        g.xf(2:end) );
    
    % Distance between neighboring cell centers
    g.dxc = diff(g.x);
    
    % Half-cell distances at external boundaries
    g.dxLeftBoundary  = g.x(1)   - g.xf(1);
    g.dxRightBoundary = g.xf(end)- g.x(end);
    
    
    %% ------------------------------------------------------------------------
    %  6. Assign material ID to each cell
    % -------------------------------------------------------------------------
    %
    %  layerId:
    %
    %      1 = aGDL
    %      2 = aCL
    %      3 = PEM
    %      4 = cCL
    %      5 = cGDL
    %
    
    g.layerId = repelem(1:g.nLayer,nCellLayer).';
    
    
    %% ------------------------------------------------------------------------
    %  7. Cell indices for each material region
    % -------------------------------------------------------------------------
    
    g.idx_aGDL = find(g.layerId == 1);
    g.idx_aCL  = find(g.layerId == 2);
    g.idx_PEM  = find(g.layerId == 3);
    g.idx_cCL  = find(g.layerId == 4);
    g.idx_cGDL = find(g.layerId == 5);
    
    
    %% ------------------------------------------------------------------------
    %  8. Physical-domain indices
    % -------------------------------------------------------------------------
    
    % Hydrogen transport:
    % anode GDL + anode catalyst layer
    g.idx_H2 = [ ...
        g.idx_aGDL;
        g.idx_aCL ];
    
    % Oxygen transport:
    % cathode catalyst layer + cathode GDL
    g.idx_O2 = [ ...
        g.idx_cCL;
        g.idx_cGDL ];
    
    % Porous regions where liquid water / ice blockage may be considered
    g.idx_porous = [ ...
        g.idx_aGDL;
        g.idx_aCL;
        g.idx_cCL;
        g.idx_cGDL ];
    
    % Catalyst layers
    g.idx_CL = [ ...
        g.idx_aCL;
        g.idx_cCL ];
    
    
    %% ------------------------------------------------------------------------
    %  9. Boolean masks
    % -------------------------------------------------------------------------
    
    g.is_aGDL = false(g.N,1);
    g.is_aCL  = false(g.N,1);
    g.is_PEM  = false(g.N,1);
    g.is_cCL  = false(g.N,1);
    g.is_cGDL = false(g.N,1);
    
    g.is_aGDL(g.idx_aGDL) = true;
    g.is_aCL (g.idx_aCL ) = true;
    g.is_PEM (g.idx_PEM ) = true;
    g.is_cCL (g.idx_cCL ) = true;
    g.is_cGDL(g.idx_cGDL) = true;
    
    g.isPorous = false(g.N,1);
    g.isPorous(g.idx_porous) = true;
    
    g.isCL = false(g.N,1);
    g.isCL(g.idx_CL) = true;
    
    
    %% ------------------------------------------------------------------------
    %  10. Material-interface face indices
    % -------------------------------------------------------------------------
    %
    %  Example:
    %
    %       cell 4 | face 5 | cell 5
    %
    %  Therefore the interface after layer k is:
    %
    %       cumulative cell number + 1
    %
    
    g.interfaceFace = cumsum(nCellLayer(1:end-1)) + 1;
    
    g.face_aGDL_aCL = g.interfaceFace(1);
    g.face_aCL_PEM  = g.interfaceFace(2);
    g.face_PEM_cCL  = g.interfaceFace(3);
    g.face_cCL_cGDL = g.interfaceFace(4);
    
    
    %% ------------------------------------------------------------------------
    %  11. Per-unit-area control-volume volume
    % -------------------------------------------------------------------------
    %
    %  Using unit cross-sectional area:
    %
    %       V_i / A = dx_i
    %
    %  This will later be useful for finite-volume balances.
    %
    
    g.volPerArea = g.dx;
    
    
    %% ------------------------------------------------------------------------
    %  12. Convenient micrometer coordinates for plotting
    % -------------------------------------------------------------------------
    
    g.x_um  = g.x  * 1e6;
    g.xf_um = g.xf * 1e6;
    g.layerBoundary_um = g.layerBoundary * 1e6;
    
    
    %% ------------------------------------------------------------------------
    %  13. Sanity checks
    % -------------------------------------------------------------------------
    
    tol = 1e-12;
    
    if abs(sum(g.dx) - g.Ltotal) > tol
        error('Grid length does not match total PEMFC thickness.');
    end
    
    if any(g.dx <= 0)
        error('Invalid grid: cell width must be positive.');
    end

    if numel(g.xf) ~= g.N + 1 || numel(g.x) ~= g.N
        error('Invalid grid: inconsistent numbers of faces and cells.');
    end

    if any(abs(g.xf(g.interfaceFace) - ...
            g.layerBoundary(2:end-1).') > tol)
        error('Material-interface faces do not match layer boundaries.');
    end

end
