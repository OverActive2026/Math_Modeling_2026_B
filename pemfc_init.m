function [x0, out] = pemfc_init(g, s, p, opt)
%PEMFC_INIT  初始化一维 PEMFC 状态向量
%
%   [x0,out] = pemfc_init(g,s,p)
%   [x0,out] = pemfc_init(g,s,p,opt)
%
% 状态向量结构：
%
%   x = [ T ; cH2 ; cO2 ; mw ; mi ]
%
% 其中：
%
%   T    : 全区域温度 [K]
%   cH2  : aGDL + aCL 中 H2 浓度 [mol/m^3]
%   cO2  : cCL + cGDL 中 O2 浓度 [mol/m^3]
%   mw   : 全区域总水质量浓度 [kg/m^3]
%   mi   : 多孔区域冰质量浓度 [kg/m^3]
%
% -------------------------------------------------------------------------
% 默认初始条件
% -------------------------------------------------------------------------
%
% 1. 温度均匀：
%
%       T0 = p.ic.T0
%
% 2. 气体浓度由理想气体状态方程计算：
%
%       c = y * p / (R*T)
%
% 3. GDL 初始经过吹扫：
%
%       mw = 0
%
% 4. PEM 初始含水量：
%
%       lambda = lambda0
%
%       mw = lambda*rho_pem*Mw/EW
%
% 5. CL 中存在离聚物，因此不能简单令 mw=0。
%
%    MVP 假设 CL 离聚物初始 lambda 与 PEM 相同：
%
%       lambda_CL,0 = lambda0
%
%    折算到整个 CL 体积：
%
%       mw_CL = phiIonCL*lambda0*rho_pem*Mw/EW
%
%    因此初始 CL 孔隙仍然是干的，只包含离聚物吸附水。
%
% 6. 初始无冰：
%
%       mi = 0
%
%
% opt 可选字段：
%
%   opt.T0
%   opt.lambda0
%   opt.cH20
%   opt.cO20
%   opt.mw0
%   opt.mi0
%
% 标量会自动扩展；向量必须与对应状态区域大小一致。
% opt.lambda0 与 opt.mw0 分别表示“由 lambda 生成总水”和
% “直接指定总水”，二者不能同时给定。


%% 1. 可选参数结构

if nargin < 4 || isempty(opt)
    opt = struct();
end

if ~isstruct(opt) || ~isscalar(opt)
    error('pemfc_init:InvalidOptions', ...
        'opt must be a scalar structure.');
end

if isfield(opt,'lambda0') && isfield(opt,'mw0')
    error('pemfc_init:AmbiguousWaterInitialization', ...
        'Specify either opt.lambda0 or opt.mw0, not both.');
end


%% 2. 初始温度

if isfield(opt,'T0')
    T0 = expand_field(opt.T0, g.N, 'opt.T0');
else
    T0 = expand_field(p.ic.T0, g.N, 'p.ic.T0');
end

if any(T0 <= 0) || any(~isfinite(T0))
    error('pemfc_init:InvalidTemperature', ...
        'Initial temperature must be finite and greater than 0 K.');
end


%% 3. 初始 H2 浓度
%
% 阳极侧默认认为初始气体组成与入口一致：
%
%       cH2 = yH2*pAnode/(R*T)

if isfield(opt,'cH20')

    cH20 = expand_field(opt.cH20, numel(g.idx_H2), 'opt.cH20');

else

    TH2 = T0(g.idx_H2);
    cH20 = p.bc.yH2 * p.bc.pAnode ./ (p.const.R * TH2);

end

if any(cH20 < 0) || any(~isfinite(cH20))
    error('pemfc_init:InvalidH2', ...
        'Initial H2 concentration must be finite and non-negative.');
end


%% 4. 初始 O2 浓度
%
% 阴极侧同理：
%
%       cO2 = yO2*pCathode/(R*T)

if isfield(opt,'cO20')

    cO20 = expand_field(opt.cO20, numel(g.idx_O2), 'opt.cO20');

else

    TO2 = T0(g.idx_O2);
    cO20 = p.bc.yO2 * p.bc.pCathode ./ (p.const.R * TO2);

end

if any(cO20 < 0) || any(~isfinite(cO20))
    error('pemfc_init:InvalidO2', ...
        'Initial O2 concentration must be finite and non-negative.');
end


%% 5. 初始膜含水量 lambda0

if isfield(opt,'lambda0')
    lambda0 = opt.lambda0;
else
    lambda0 = p.water.lambda0;
end

if ~isscalar(lambda0) || ~isfinite(lambda0) || lambda0 < 0
    error('pemfc_init:InvalidLambda', ...
        'lambda0 must be a finite non-negative scalar.');
end


%% 6. 初始化总水 mw

if isfield(opt,'mw0')

    % 如果用户直接指定 mw0，则完全采用用户值
    mw0 = expand_field(opt.mw0, g.N, 'opt.mw0');

else

    mw0 = zeros(g.N,1);


    % ---------------------------------------------------------------------
    % 6.1 GDL：吹扫后初始孔隙水 = 0
    % ---------------------------------------------------------------------

    mw0(g.idx_aGDL) = 0;
    mw0(g.idx_cGDL) = 0;


    % ---------------------------------------------------------------------
    % 6.2 PEM：吸附态水
    %
    % lambda = EW*mw/(rho_pem*Mw)
    %
    % 所以：
    %
    % mw = lambda*rho_pem*Mw/EW
    % ---------------------------------------------------------------------

    mwPEM0 = lambda0 * p.water.rho_pem * p.const.Mw / p.water.EW;

    mw0(g.idx_PEM) = mwPEM0;


    % ---------------------------------------------------------------------
    % 6.3 CL：离聚物中的吸附水
    %
    % CL 内只有 phiIonCL 比例的总体积属于离聚物，因此：
    %
    % mwCL0 = phiIonCL * mwPEM0
    %
    % 这里假设：
    %
    % lambda_CL,0 = lambda_PEM,0
    %
    % 这是当前 MVP 的建模假设。
    % ---------------------------------------------------------------------

    phiIonCL = p.water.clIonomerVolumeFraction;

    if ~isscalar(phiIonCL) || ~isfinite(phiIonCL) || ...
            phiIonCL < 0 || phiIonCL > 1
        error('pemfc_init:InvalidCLIonomerFraction', ...
            ['p.water.clIonomerVolumeFraction must be a finite ', ...
             'scalar between 0 and 1.']);
    end

    mwCL0 = phiIonCL * mwPEM0;

    mw0(g.idx_aCL) = mwCL0;
    mw0(g.idx_cCL) = mwCL0;

end

if any(mw0 < 0) || any(~isfinite(mw0))
    error('pemfc_init:InvalidWater', ...
        'Initial total-water concentration must be finite and non-negative.');
end


%% 7. 初始冰
%
% MVP 和附件中的吹扫初值均采用：
%
%       mi = 0

if isfield(opt,'mi0')
    mi0 = expand_field(opt.mi0, numel(s.idx.mi), 'opt.mi0');
else
    mi0 = zeros(numel(s.idx.mi),1);
end

if any(mi0 < 0) || any(~isfinite(mi0))
    error('pemfc_init:InvalidIce', ...
        'Initial ice concentration must be finite and non-negative.');
end


%% 7.1 检查冰与总水、孔隙体积的一致性
%
% mw 是总水：
%
%       mw = mv + ml + mi
%
% 因此冰质量不能超过对应网格内的总水质量。

mwPorous0 = mw0(g.idx_porous);

if any(mi0 > mwPorous0)
    error('pemfc_init:IceExceedsTotalWater', ...
        ['Initial ice concentration cannot exceed total-water ', ...
         'concentration in a porous cell.']);
end


% 冰所占体积不能超过该材料的干孔隙体积。
% p.mat.eps0 按 [aGDL,aCL,PEM,cCL,cGDL] 的层序存储。

porousLayerId = g.layerId(g.idx_porous);
eps0Porous = p.mat.eps0(porousLayerId(:));
eps0Porous = eps0Porous(:);

rhoIce = p.ice.rho_ice;

if ~isscalar(rhoIce) || ~isfinite(rhoIce) || rhoIce <= 0
    error('pemfc_init:InvalidIceDensity', ...
        'p.ice.rho_ice must be a finite positive scalar.');
end

epsIce0 = mi0 ./ rhoIce;

if any(epsIce0 > eps0Porous)
    error('pemfc_init:IceExceedsPoreVolume', ...
        'Initial ice volume exceeds the available dry pore volume.');
end


%% 8. 打包状态向量

x0 = zeros(s.Nx,1);

x0(s.idx.T)   = T0;
x0(s.idx.H2)  = cH20;
x0(s.idx.O2)  = cO20;
x0(s.idx.mw)  = mw0;
x0(s.idx.mi)  = mi0;


%% 9. 初始化结果检查

if any(~isfinite(x0))
    error('pemfc_init:NonFiniteState', ...
        'Initial state contains NaN or Inf.');
end


%% 10. 调试输出

if nargout > 1

    out.T0   = T0;
    out.cH20 = cH20;
    out.cO20 = cO20;
    out.mw0  = mw0;
    out.mi0  = mi0;

    % PEM 初始吸附水
    out.mwPEM0 = mean(mw0(g.idx_PEM));


    % 两个 CL 初始总水
    out.mwACL0 = mean(mw0(g.idx_aCL));
    out.mwCCL0 = mean(mw0(g.idx_cCL));


    % GDL 初始水
    out.mwAGDL0 = mean(mw0(g.idx_aGDL));
    out.mwCGDL0 = mean(mw0(g.idx_cGDL));


    % 检查由 mw 反算出的 PEM lambda
    out.lambdaPEMCheck = membrane_lambda(mw0(g.idx_PEM), p);
    out.lambdaPEMMean  = mean(out.lambdaPEMCheck);

    % lambda0 始终表示最终初始状态的 PEM 平均含水量。
    % lambdaParameter 保留由 opt.lambda0 或 p.water.lambda0 读入的值；
    % 当直接指定 mw0 时，该参数并未参与初始水分布的构造。
    out.lambda0 = out.lambdaPEMMean;
    out.lambdaParameter = lambda0;
    out.usedDirectWaterInitialization = isfield(opt,'mw0');


    % 初始整体平均温度
    out.Tavg0 = sum(T0 .* g.dx) / sum(g.dx);

end

end


%% =========================================================================
% 辅助函数：标量自动扩展为列向量
% =========================================================================

function value = expand_field(value, N, name)

if ~isnumeric(value) || ~isreal(value)
    error('pemfc_init:InvalidFieldType', ...
        '%s must be a real numeric scalar or vector.', name);
end

if isscalar(value)

    value = value * ones(N,1);

else

    value = value(:);

    if numel(value) ~= N
        error('pemfc_init:SizeMismatch', ...
            '%s must be scalar or contain %d elements.', name, N);
    end

end

end
