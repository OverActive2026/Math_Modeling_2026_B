function [Nw, out] = calc_water_flux(T, mw, j, g, p, bc)
%CALC_WATER_FLUX  分区域计算 PEMFC 一维水通量
%
%   [Nw,out] = calc_water_flux(T,mw,j,g,p)
%   [Nw,out] = calc_water_flux(T,mw,j,g,p,bc)
%
% 状态变量仍统一采用：
%
%       mw(x,t)   [kg/m^3]
%
% 但不同区域采用不同的实际输运变量：
%
%   aGDL + aCL : 孔隙水蒸气 mv 扩散
%   PEM        : 吸附态水 mw 扩散
%   cCL + cGDL : 孔隙水蒸气 mv 扩散
%
% 总水通量：
%
%       Nw = Ndiff + Neod
%
%       Neod = nd * Mw * im / F
%
% CL 中将总水拆分为：
%
%       mw = mIonCL + mPore
%
% 其中 mIonCL 为离聚物吸附水，mPore 再进行 vapor/liquid 相平衡。
%
% 当前 MVP 不直接用 CL 水蒸气浓度和 PEM 吸附水浓度做 Fick 扩散。
% 两个 CL/PEM 界面的扩散通量默认设为 0，只保留电渗拖曳通量。
%
% 后续若建立吸附/解吸界面模型，只需修改两个界面 face 的 Ndiff。
%
% -------------------------------------------------------------------------
% 输入
% -------------------------------------------------------------------------
% T   : 全局温度场 [K]
% mw  : 全局总水质量浓度 [kg/m^3]
% j   : 电流密度 [A/m^2]
% g   : build_grid() 输出
% p   : pemfc_params() 输出
%
% bc.left / bc.right:
%   多孔介质外部水蒸气边界条件
%
% bc.interface.aCL_PEM:
%   aCL/PEM 界面扩散水通量 [kg/(m^2 s)]
%
% bc.interface.PEM_cCL:
%   PEM/cCL 界面扩散水通量 [kg/(m^2 s)]
%
% 所有通量规定：
%
%       正方向 = +x = 阳极 -> 阴极
%
% -------------------------------------------------------------------------
% 输出
% -------------------------------------------------------------------------
% Nw            : 总水通量，global faces，长度 g.N+1
%
% out.Ndiff     : 扩散部分
% out.Neod      : 电渗拖曳部分
%
% out.mv        : 孔隙水蒸气浓度
% out.ml        : 孔隙液水浓度
% out.mIonCL    : CL 离聚物吸附水
% out.mPore     : CL/GDL 孔隙水总量
%
% out.Dw        : 各区域有效水扩散系数
% out.eps_g     : 剩余气相孔隙率
%
% out.lambdaPEM
% out.lambdaMean
% out.nd


%% 1. 输入检查

T  = T(:);
mw = mw(:);

if numel(T) ~= g.N || numel(mw) ~= g.N
    error('calc_water_flux:SizeMismatch', 'T and mw must both contain g.N elements.');
end

if any(~isfinite(T)) || any(T <= 0)
    error('calc_water_flux:InvalidTemperature', 'T must contain valid Kelvin temperatures.');
end

if any(~isfinite(mw)) || any(mw < 0)
    error('calc_water_flux:InvalidWater', 'mw must be finite and non-negative.');
end

if ~isscalar(j) || ~isfinite(j) || j < 0
    error('calc_water_flux:InvalidCurrent', 'j must be a finite non-negative scalar.');
end


%% 2. 默认外部边界：干氢气、干空气

if nargin < 6 || isempty(bc)
    bc = struct();
end

if ~isfield(bc,'left')
    bc.left.type  = 'dirichlet';
    bc.left.value = 0;
end

if ~isfield(bc,'right')
    bc.right.type  = 'dirichlet';
    bc.right.value = 0;
end


%% 3. 建立干孔隙率 eps0(x)

eps0 = zeros(g.N,1);

for k = 1:5
    idxLayer = find(g.layerId == k);
    eps0(idxLayer) = p.mat.eps0(k);
end


%% 4. PEM 膜含水量 lambda

lambdaPEM = membrane_lambda(mw(g.idx_PEM), p);

dxPEM = g.dx(g.idx_PEM);
lambdaMean = sum(lambdaPEM .* dxPEM) / sum(dxPEM);


%% 5. CL 中的离聚物水与孔隙水分离
%
% 附件1给定 CL 离聚物体积分数约为 0.3。
%
% MVP 假设：
% CL 离聚物的平衡 lambda 暂时取当前 PEM 平均 lambda。
%
% 离聚物本体中的水浓度：
%
%       mIonMaterial = lambda * rho_pem * Mw / EW
%
% 折算到整个 CL 总体积：
%
%       mIonTarget = phiIonCL * mIonMaterial

if isfield(p.water,'phiIonomerCL')
    phiIonCL = p.water.phiIonomerCL;
else
    phiIonCL = 0.3;
end

mIonTarget = phiIonCL * lambdaMean * p.water.rho_pem * p.const.Mw / p.water.EW;

mIonCL = zeros(g.N,1);
mPore  = zeros(g.N,1);

% GDL 中没有离聚物，总水全部属于孔隙水
mPore(g.idx_aGDL) = mw(g.idx_aGDL);
mPore(g.idx_cGDL) = mw(g.idx_cGDL);

% CL 中先分配离聚物水，剩余部分才属于孔隙水
idxCL = [g.idx_aCL; g.idx_cCL];

mIonCL(idxCL) = min(mw(idxCL), mIonTarget);
mPore(idxCL)  = max(mw(idxCL) - mIonCL(idxCL), 0);


%% 6. 多孔区域中：孔隙水进一步分成 vapor / liquid

mv    = zeros(g.N,1);
ml    = zeros(g.N,1);
eps_g = zeros(g.N,1);

idxA = [g.idx_aGDL; g.idx_aCL];
idxC = [g.idx_cCL; g.idx_cGDL];

% ---------- 阳极多孔区域 ----------
miA = zeros(numel(idxA),1);

[mvA, mlA, epsA] = water_phase_equilibrium(mPore(idxA), miA, T(idxA), eps0(idxA), p);

mv(idxA)    = mvA;
ml(idxA)    = mlA;
eps_g(idxA) = epsA;

% ---------- 阴极多孔区域 ----------
miC = zeros(numel(idxC),1);

[mvC, mlC, epsC] = water_phase_equilibrium(mPore(idxC), miC, T(idxC), eps0(idxC), p);

mv(idxC)    = mvC;
ml(idxC)    = mlC;
eps_g(idxC) = epsC;

% PEM 没有气相孔隙
eps_g(g.idx_PEM) = 0;


%% 7. 各区域水扩散系数

Dw = zeros(g.N,1);

% 阳极多孔区域：水蒸气扩散
Dw(idxA) = calc_gas_diffusivity(p.water.Dw_ref_anode, T(idxA), p.bc.pAnode, eps_g(idxA), p);

% PEM：吸附态水扩散
Dw(g.idx_PEM) = membrane_water_diffusivity(lambdaPEM, T(g.idx_PEM), p);

% 阴极多孔区域：水蒸气扩散
Dw(idxC) = calc_gas_diffusivity(p.water.Dw_ref_cathode, T(idxC), p.bc.pCathode, eps_g(idxC), p);


%% 8. 分区域计算扩散水通量
%
% 注意：
%
% 阳极多孔区扩散变量 = mv
% PEM 扩散变量        = mw
% 阴极多孔区扩散变量 = mv
%
% CL/PEM 界面不直接做普通 Fick 扩散。

Ndiff = zeros(g.N+1,1);


% ========================================================================
% 8.1 阳极多孔区域：aGDL + aCL
% ========================================================================

bcA.left = bc.left;
bcA.right.type = 'noflux';

[~, NdiffA, diffInfoA] = diffusion_fv(mv(idxA), Dw(idxA), g, idxA, bcA);

faceA = idxA(1):idxA(end)+1;
Ndiff(faceA) = NdiffA;


% ========================================================================
% 8.2 PEM
% ========================================================================

idxM = g.idx_PEM;

bcM.left.type  = 'noflux';
bcM.right.type = 'noflux';

[~, NdiffM, diffInfoM] = diffusion_fv(mw(idxM), Dw(idxM), g, idxM, bcM);

faceM = idxM(1):idxM(end)+1;
Ndiff(faceM) = NdiffM;


% ========================================================================
% 8.3 阴极多孔区域：cCL + cGDL
% ========================================================================

bcC.left.type = 'noflux';
bcC.right = bc.right;

[~, NdiffC, diffInfoC] = diffusion_fv(mv(idxC), Dw(idxC), g, idxC, bcC);

faceC = idxC(1):idxC(end)+1;
Ndiff(faceC) = NdiffC;


%% 9. CL / PEM 界面扩散通量
%
% 默认：
%
%       Ndiff,aCL/PEM = 0
%       Ndiff,PEM/cCL = 0
%
% 这样避免把：
%
%       水蒸气浓度 mv
%
% 和：
%
%       膜吸附水浓度 mw
%
% 直接做浓度梯度。
%
% 后续加入 sorption/desorption 模型时修改这里即可。

Nint_a = 0;
Nint_c = 0;

if isfield(bc,'interface')

    if isfield(bc.interface,'aCL_PEM')
        Nint_a = bc.interface.aCL_PEM;
    end

    if isfield(bc.interface,'PEM_cCL')
        Nint_c = bc.interface.PEM_cCL;
    end

end

Ndiff(g.face_aCL_PEM) = Nint_a;
Ndiff(g.face_PEM_cCL) = Nint_c;


%% 10. 质子电流分布

[~, imCell, ~, currentInfo] = calc_current_distribution(j, g, p);


%% 11. 电渗拖曳系数
%
%       nd = 2.5 * lambda / 22
%
% MVP 暂时采用 PEM 平均 lambda。

nd = 2.5 * lambdaMean / 22;


%% 12. 计算 face 上的质子电流

imFace = proton_current_faces(j, g, p);


%% 13. 电渗拖曳水通量
%
%       Neod = nd * Mw * im / F

Neod = nd * p.const.Mw .* imFace / p.const.F;


%% 14. 总水通量

Nw = Ndiff + Neod;


%% 15. 输出调试信息

out.Ndiff = Ndiff;
out.Neod  = Neod;

out.Dw = Dw;

out.mv = mv;
out.ml = ml;

out.mIonCL = mIonCL;
out.mPore  = mPore;

out.eps_g = eps_g;

out.lambdaPEM  = lambdaPEM;
out.lambdaMean = lambdaMean;
out.nd         = nd;

out.imCell = imCell;
out.imFace = imFace;

out.diffusionInfo.anode    = diffInfoA;
out.diffusionInfo.membrane = diffInfoM;
out.diffusionInfo.cathode  = diffInfoC;

out.interfaceFlux.aCL_PEM_diff = Nint_a;
out.interfaceFlux.PEM_cCL_diff = Nint_c;

out.interfaceFlux.aCL_PEM_total = Nw(g.face_aCL_PEM);
out.interfaceFlux.PEM_cCL_total = Nw(g.face_PEM_cCL);

out.interfaceModel = 'zero diffusive CL/PEM flux + EOD (MVP closure)';
out.currentInfo = currentInfo;

end


%% =========================================================================
% 辅助函数：直接在 face 上计算质子电流 im
% =========================================================================

function imFace = proton_current_faces(j, g, p)

x = g.xf;

imFace = zeros(size(x));

x1 = g.layerBoundary(2);   % aGDL / aCL
x2 = g.layerBoundary(3);   % aCL  / PEM
x3 = g.layerBoundary(4);   % PEM  / cCL
x4 = g.layerBoundary(5);   % cCL  / cGDL


%% aGDL：im = 0

imFace(x <= x1) = 0;


%% aCL：im 从 0 线性上升到 j

idx = x > x1 & x < x2;
xi = (x(idx) - x1) / p.geom.LaCL;

imFace(idx) = j .* xi;


%% aCL / PEM 界面

tol2 = 100 * eps(max(abs(x2),1));
imFace(abs(x-x2) <= tol2) = j;


%% PEM：im = j

idx = x > x2 & x < x3;
imFace(idx) = j;


%% PEM / cCL 界面

tol3 = 100 * eps(max(abs(x3),1));
imFace(abs(x-x3) <= tol3) = j;


%% cCL：im 从 j 线性下降到 0

idx = x > x3 & x < x4;
xi = (x(idx) - x3) / p.geom.LcCL;

imFace(idx) = j .* (1-xi);


%% cGDL：im = 0

imFace(x >= x4) = 0;

end