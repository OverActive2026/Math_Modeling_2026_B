%-------------------------------------------------------------------------------
% @function: 构建第一版含冰PEMFC模型的网格、状态编号和初始状态
% @author:   PJ, GPT
% @date:     20260923
% @input:    nCellLayer->五层网格数，默认[6 13 19 19 6]
%            opt->可选初值：T0、lambda0、lambdaCCL0、cH20、cO20、mw0、mi0
% @output:   g->空间网格
%            s->状态编号
%            x0->初始状态
%            init->初值汇总
%
% [ICE-1] 初始冰处理：
%   吹扫后的基准工况令所有多孔层mi0=0；PEM仍保留lambda0=3的结合水。
%   冰状态只定义在aGDL、aCL、cCL和cGDL中，不在PEM内设置孔隙冰。
% [HYD-1] cCL离聚物保留一个层平均lambda状态。总水mw中包含
%   离聚物结合水；吹扫只清除孔隙自由水。
%   [INIT-2] 综合问题1初始电压和问题2恒流可行性，
%   cCL初始值从未经题目给定的1修正为2；PEM初始lambda0仍为3。
%-------------------------------------------------------------------------------

function [g,s,x0,init] = pemfc_setup_ice(nCellLayer,opt)

    %% 1. 输入处理

    if nargin < 1 || isempty(nCellLayer)
        nCellLayer = [6 13 19 19 6];
    end
    if nargin < 2 || isempty(opt)
        opt = struct();
    end
    if ~isnumeric(nCellLayer) || numel(nCellLayer) ~= 5 || ...
            any(~isfinite(nCellLayer)) || any(nCellLayer <= 0) || ...
            any(mod(nCellLayer,1) ~= 0)
        error('pemfc_setup_ice:InvalidGrid','nCellLayer必须包含5个正整数。');
    end
    if ~isstruct(opt) || ~isscalar(opt)
        error('pemfc_setup_ice:InvalidOptions','opt必须是标量结构体。');
    end
    if isfield(opt,'lambda0') && isfield(opt,'mw0')
        error('pemfc_setup_ice:AmbiguousWaterInitialization', ...
            'lambda0和mw0不能同时指定。');
    end
    nCellLayer = nCellLayer(:).';


    %% 2. 五层有限体积网格

    g.layerNames = {'aGDL','aCL','PEM','cCL','cGDL'};
    g.layerThickness = [150e-6,3.4e-6,12e-6,11.3e-6,150e-6];
    g.nLayer = 5;
    g.nCellLayer = nCellLayer;
    g.N = sum(nCellLayer);
    g.layerBoundary = [0,cumsum(g.layerThickness)];
    g.Ltotal = g.layerBoundary(end);

    g.xf = zeros(g.N+1,1);
    faceCounter = 1;
    for k = 1:g.nLayer
        localFace = linspace(g.layerBoundary(k), ...
            g.layerBoundary(k+1),nCellLayer(k)+1).';
        if k == 1
            g.xf(1:nCellLayer(k)+1) = localFace;
            faceCounter = nCellLayer(k)+1;
        else
            g.xf(faceCounter+1:faceCounter+nCellLayer(k)) = localFace(2:end);
            faceCounter = faceCounter+nCellLayer(k);
        end
    end

    g.dx = diff(g.xf);
    g.x = 0.5*(g.xf(1:end-1)+g.xf(2:end));
    g.dxc = diff(g.x);
    g.layerId = repelem(1:g.nLayer,nCellLayer).';
    g.idx_aGDL = find(g.layerId == 1);
    g.idx_aCL = find(g.layerId == 2);
    g.idx_PEM = find(g.layerId == 3);
    g.idx_cCL = find(g.layerId == 4);
    g.idx_cGDL = find(g.layerId == 5);
    g.idx_H2 = [g.idx_aGDL;g.idx_aCL];
    g.idx_O2 = [g.idx_cCL;g.idx_cGDL];
    g.idx_porous = [g.idx_aGDL;g.idx_aCL;g.idx_cCL;g.idx_cGDL];
    g.idx_CL = [g.idx_aCL;g.idx_cCL];
    g.isPorous = ismember(g.layerId,[1 2 4 5]);
    g.isCL = ismember(g.layerId,[2 4]);

    g.interfaceFace = cumsum(nCellLayer(1:end-1))+1;
    g.face_aGDL_aCL = g.interfaceFace(1);
    g.face_aCL_PEM = g.interfaceFace(2);
    g.face_PEM_cCL = g.interfaceFace(3);
    g.face_cCL_cGDL = g.interfaceFace(4);
    g.volPerArea = g.dx;
    g.x_um = g.x*1e6;
    g.xf_um = g.xf*1e6;
    g.layerBoundary_um = g.layerBoundary*1e6;


    %% 3. 状态空间编号
    % x=[T;H2;O2;mw;mi;lambdaCCL]，其中mi只覆盖四个多孔层。

    s.n.T = g.N;
    s.n.H2 = numel(g.idx_H2);
    s.n.O2 = numel(g.idx_O2);
    s.n.mw = g.N;
    s.n.mi = numel(g.idx_porous);
    s.n.lambdaCCL = 1;

    cursor = 0;
    s.idx.T = cursor+(1:s.n.T);       cursor = cursor+s.n.T;
    s.idx.H2 = cursor+(1:s.n.H2);     cursor = cursor+s.n.H2;
    s.idx.O2 = cursor+(1:s.n.O2);     cursor = cursor+s.n.O2;
    s.idx.mw = cursor+(1:s.n.mw);     cursor = cursor+s.n.mw;
    s.idx.mi = cursor+(1:s.n.mi);     cursor = cursor+s.n.mi;
    s.idx.lambdaCCL = cursor+1;        cursor = cursor+s.n.lambdaCCL;
    s.Nx = cursor;

    s.globalToLocal.H2 = zeros(g.N,1);
    s.globalToLocal.O2 = zeros(g.N,1);
    s.globalToLocal.mi = zeros(g.N,1);
    s.globalToLocal.H2(g.idx_H2) = (1:s.n.H2).';
    s.globalToLocal.O2(g.idx_O2) = (1:s.n.O2).';
    s.globalToLocal.mi(g.idx_porous) = (1:s.n.mi).';
    s.local.H2_aCL = s.globalToLocal.H2(g.idx_aCL);
    s.local.O2_cCL = s.globalToLocal.O2(g.idx_cCL);
    s.local.mi_aGDL = s.globalToLocal.mi(g.idx_aGDL);
    s.local.mi_aCL = s.globalToLocal.mi(g.idx_aCL);
    s.local.mi_cCL = s.globalToLocal.mi(g.idx_cCL);
    s.local.mi_cGDL = s.globalToLocal.mi(g.idx_cGDL);
    s.names = {'T','H2','O2','mw','mi','lambdaCCL'};


    %% 4. 温度和气体初值

    if isfield(opt,'T0'), T0 = opt.T0; else, T0 = 253.15; end
    if isscalar(T0), T0 = T0*ones(g.N,1); else, T0 = T0(:); end
    if numel(T0) ~= g.N || any(~isfinite(T0)) || any(T0 <= 0)
        error('pemfc_setup_ice:InvalidTemperature','T0的长度或数值不合法。');
    end

    cathodeMoles = 0.233/31.998e-3+0.767/28.014e-3;
    yO2 = (0.233/31.998e-3)/cathodeMoles;

    if isfield(opt,'cH20')
        cH20 = opt.cH20;
        if isscalar(cH20), cH20 = cH20*ones(s.n.H2,1); else, cH20 = cH20(:); end
    else
        cH20 = 101325./(8.314*T0(g.idx_H2));
    end

    if isfield(opt,'cO20')
        cO20 = opt.cO20;
        if isscalar(cO20), cO20 = cO20*ones(s.n.O2,1); else, cO20 = cO20(:); end
    else
        cO20 = yO2*101325./(8.314*T0(g.idx_O2));
    end

    if numel(cH20) ~= s.n.H2 || any(~isfinite(cH20)) || any(cH20 <= 0)
        error('pemfc_setup_ice:InvalidHydrogen', ...
            'cH20必须是正有限标量或长度为氢气网格数的向量。');
    end
    if numel(cO20) ~= s.n.O2 || any(~isfinite(cO20)) || any(cO20 <= 0)
        error('pemfc_setup_ice:InvalidOxygen', ...
            'cO20必须是正有限标量或长度为氧气网格数的向量。');
    end


    %% 5. 总水和冰初值

    if isfield(opt,'lambda0'), lambda0 = opt.lambda0; else, lambda0 = 3.0; end
    if ~isscalar(lambda0) || ~isfinite(lambda0) || lambda0 < 0
        error('pemfc_setup_ice:InvalidLambda', ...
            'lambda0必须是非负有限标量。');
    end

    if isfield(opt,'lambdaCCL0')
        lambdaCCL0 = opt.lambdaCCL0;
    elseif isfield(opt,'mw0')
        % 兼容用户直接给定总水场的情况：未指定时不额外凭空加入结合水。
        lambdaCCL0 = 0;
    else
        lambdaCCL0 = 2.0;
    end
    if ~isscalar(lambdaCCL0) || ~isfinite(lambdaCCL0) || ...
            lambdaCCL0 < 0 || lambdaCCL0 > 22
        error('pemfc_setup_ice:InvalidCCLLambda', ...
            'lambdaCCL0必须是[0,22]内的有限标量。');
    end

    % [HYD-1] 附件1给出CL离聚物体积分数0.3。单位控制体体积中，
    % cCL结合水质量浓度为0.3*rhoIonomer*Mw/EW*lambdaCCL。
    cCLIonomerWaterCoefficient = 0.3*2150*0.018/1.0;

    if isfield(opt,'mw0')
        mw0 = opt.mw0;
        if isscalar(mw0), mw0 = mw0*ones(g.N,1); else, mw0 = mw0(:); end
    else
        mw0 = zeros(g.N,1);
        mw0(g.idx_PEM) = lambda0*2150*0.018/1.0;
        mw0(g.idx_cCL) = cCLIonomerWaterCoefficient*lambdaCCL0;
    end
    if numel(mw0) ~= g.N || any(~isfinite(mw0)) || any(mw0 < 0)
        error('pemfc_setup_ice:InvalidInitialWater', ...
            'mw0必须是非负有限标量或长度为总网格数的向量。');
    end

    % [ICE-1] 吹扫基准：所有多孔层初始冰为0。
    if isfield(opt,'mi0')
        mi0 = opt.mi0;
        if isscalar(mi0), mi0 = mi0*ones(s.n.mi,1); else, mi0 = mi0(:); end
    else
        mi0 = zeros(s.n.mi,1);
    end
    if numel(mi0) ~= s.n.mi || any(~isfinite(mi0))
        error('pemfc_setup_ice:InvalidInitialIceSize', ...
            'mi0必须是有限标量或长度为多孔网格数的向量。');
    end

    ionomerWater0 = zeros(g.N,1);
    ionomerWater0(g.idx_cCL) = ...
        cCLIonomerWaterCoefficient*lambdaCCL0;
    poreWater0 = mw0-ionomerWater0;
    if any(poreWater0(g.idx_porous) < -1e-10)
        error('pemfc_setup_ice:InsufficientCCLWater', ...
            'mw0中的cCL总水小于lambdaCCL0对应的离聚物结合水。');
    end
    poreWater0 = max(poreWater0,0);

    eps0Layer = [0.8,0.3916,0,0.4207,0.8];
    eps0Porous = eps0Layer(g.layerId(g.idx_porous));
    eps0Porous = eps0Porous(:);
    if any(mi0 < 0) || any(mi0 > poreWater0(g.idx_porous)) || ...
            any(mi0/920 > eps0Porous)
        error('pemfc_setup_ice:InvalidInitialIce', ...
            '初始冰必须非负，且不能超过孔隙水或干孔隙容量。');
    end

    x0 = zeros(s.Nx,1);
    x0(s.idx.T) = T0;
    x0(s.idx.H2) = cH20;
    x0(s.idx.O2) = cO20;
    x0(s.idx.mw) = mw0;
    x0(s.idx.mi) = mi0;
    x0(s.idx.lambdaCCL) = lambdaCCL0;

    init.Tavg0 = sum(T0.*g.dx)/sum(g.dx);
    init.lambda0 = mean(mw0(g.idx_PEM)/(2150*0.018));
    init.lambdaCCL0 = lambdaCCL0;
    init.cCLIonomerWaterCoefficient = cCLIonomerWaterCoefficient;
    init.initialCCLIonomerWaterMassPerArea = ...
        sum(ionomerWater0(g.idx_cCL).*g.dx(g.idx_cCL));
    init.mi0 = mi0;
    init.initialIceMassPerArea = sum(mi0.*g.dx(g.idx_porous));
    init.porousLayersPurged = all(mi0 == 0) && ...
        all(poreWater0(g.idx_porous) <= 1e-12);
    init.membraneWaterRetained = any(mw0(g.idx_PEM) > 0);

end
