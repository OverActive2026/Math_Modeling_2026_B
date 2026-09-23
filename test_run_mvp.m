%% test_run_mvp.m
% ================================================================
% PEMFC 一维模型 MVP 首次联调程序
%
% 目标：
%   1. 建立 grid / state index / parameters
%   2. 生成初始状态
%   3. 单独调用一次 pemfc_rhs_mvp，检查有没有 NaN / Inf
%   4. 使用 ode15s 进行短时间积分
%   5. 后处理并绘制 T、V、H2、O2、水含量
%
% 当前模型：
%   - 无冰生成
%   - mi 始终为 0
%   - 低电流测试
%
% ================================================================

clear;
clc;
close all;


%% ================================================================
% 1. 建立空间网格
% ================================================================

g = build_grid();

fprintf('Grid created.\n');
fprintf('Number of cells = %d\n', g.N);
fprintf('Total thickness = %.3f um\n', g.Ltotal*1e6);


%% ================================================================
% 2. 建立状态索引
% ================================================================

s = build_state_index(g);

fprintf('\nState index created.\n');
fprintf('Number of states = %d\n', s.Nx);


%% ================================================================
% 3. 读取模型参数
% ================================================================

p = pemfc_params();


%% ================================================================
% 4. MVP 联调所需参数检查 / 补充
% ================================================================
%
% 这一段以后等 pemfc_params.m 完全整理好以后可以删掉。
% 现在主要是为了避免因为字段没同步导致 RHS 直接报错。
% ================================================================


% ---------- 初始温度 ----------
if ~isfield(p,'op')
    p.op = struct();
end

p.op.T0 = 253.15;       % -20 degC


% ---------- 初始膜含水量 ----------
p.water.lambda0 = 3;

% CL 离聚物体积分数
p.water.phiIonomerCL = 0.3;


% ---------- 气体压力 ----------
p.bc.pAnode   = 101325;
p.bc.pCathode = 101325;


% ---------- 入口气体摩尔分数 ----------
p.bc.xH2 = 1.0;

% 附件给的是 O2/N2 质量分数 0.233 / 0.767
% 这里转换成摩尔分数

wO2 = 0.233;
wN2 = 0.767;

MO2 = 0.032;
MN2 = 0.028;

p.bc.xO2 = (wO2/MO2) / ((wO2/MO2) + (wN2/MN2));


% ---------- 阴极等效扩散距离 ----------
%
% MVP 暂取：
%
% Ldiff = LcCL + LcGDL

p.echem.Ldiff = p.geom.LcCL + p.geom.LcGDL;


% ---------- 单电池 MEA 厚度 ----------
%
% 电化学体积热源中需要：
%
% qgen = j*(Eth-V)/Lcell

p.geom.Lcell = sum(g.layerThickness);


% ---------- 热物性 ----------
%
% 顺序：
%
% [aGDL, aCL, PEM, cCL, cGDL]
%
% 这些值来自我们此前整理的附件参数：
%
% GDL:
%   rho = 185
%   cp  = 545
%   k   = 0.3
%
% CL:
%   rho = 970
%   cp  = 240
%   k   = 0.27
%
% PEM 使用离聚物物性：
%   rho = 2150
%   cp  = 1050
%   k   = 0.24

p.mat.rho = [185; 970; 2150; 970; 185];

p.mat.cp  = [545; 240; 1050; 240; 545];

p.mat.k   = [0.30; 0.27; 0.24; 0.27; 0.30];


%% ================================================================
% 5. 检查水通量函数需要的 interface face
% ================================================================
%
% 如果最新版 build_grid 已经定义这两个字段，这部分不会做任何事情。
%
% 对 cell-centered FV：
%
% cell i 与 cell i+1 之间的 face 编号为 i+1。
%

if ~isfield(g,'face_aCL_PEM')
    g.face_aCL_PEM = g.idx_aCL(end) + 1;
end

if ~isfield(g,'face_PEM_cCL')
    g.face_PEM_cCL = g.idx_PEM(end) + 1;
end


fprintf('\nInterface faces:\n');
fprintf('aCL / PEM face = %d\n', g.face_aCL_PEM);
fprintf('PEM / cCL face = %d\n', g.face_PEM_cCL);


%% ================================================================
% 6. 生成初始状态
% ================================================================

opt = struct();

opt.T0 = 253.15;
opt.lambda0 = 3;


[x0, ic] = pemfc_init(g, s, p, opt);


fprintf('\n========================================\n');
fprintf('Initial condition\n');
fprintf('========================================\n');

fprintf('Tavg0       = %.3f degC\n', ic.Tavg0 - 273.15);
fprintf('lambda PEM  = %.4f\n', ic.lambdaPEMMean);
fprintf('mw PEM      = %.4f kg/m^3\n', ic.mwPEM0);
fprintf('mw aCL      = %.4f kg/m^3\n', ic.mwACL0);
fprintf('mw cCL      = %.4f kg/m^3\n', ic.mwCCL0);
fprintf('mw aGDL     = %.4f kg/m^3\n', ic.mwAGDL0);
fprintf('mw cGDL     = %.4f kg/m^3\n', ic.mwCGDL0);


%% ================================================================
% 7. 基本初值检查
% ================================================================

assert(all(isfinite(x0)), ...
    'Initial state contains NaN or Inf.');

assert(abs(ic.lambdaPEMMean - 3) < 1e-10, ...
    'Initial PEM lambda is incorrect.');

assert(all(x0(s.idx.mi) == 0), ...
    'Initial ice state should be zero.');

fprintf('\nInitial-condition check PASSED.\n');


%% ================================================================
% 8. 设置外部输入
% ================================================================
%
% 1 A/cm^2 = 1e4 A/m^2
%
% 第一次联调建议：
%
% j = 0.01 A/cm^2
%
% 即：
%
% j = 100 A/m^2

u = struct();

u.j = 0.01e4;          % A/m^2
u.Tamb = 253.15;       % K
u.qaux = 0;            % W/m^3


fprintf('\nApplied current density = %.4f A/cm^2\n', ...
    u.j/1e4);


%% ================================================================
% 9. 第一次：不要直接 ode15s
%
% 先手动调用 RHS 一次
% ================================================================

fprintf('\n========================================\n');
fprintf('Testing RHS at t = 0\n');
fprintf('========================================\n');


[dx0, out0] = pemfc_rhs_mvp(0, x0, u, p, g, s);


%% ================================================================
% 10. RHS 健康检查
% ================================================================

if any(~isfinite(dx0))

    badIndex = find(~isfinite(dx0));

    fprintf('\nRHS FAILED.\n');
    fprintf('Non-finite states at indices:\n');
    disp(badIndex);

    error('dx0 contains NaN or Inf.');

end


fprintf('RHS finite check PASSED.\n');


fprintf('\nInitial model outputs:\n');

fprintf('Vcell       = %.6f V\n', out0.Vcell);
fprintf('Tavg        = %.4f degC\n', out0.Tavg - 273.15);
fprintf('lambdaMean  = %.6f\n', out0.lambdaMean);

fprintf('cH2,aCL     = %.6f mol/m^3\n', out0.cH2_aCL);
fprintf('cO2,cCL     = %.6f mol/m^3\n', out0.cO2_cCL);

fprintf('pH2         = %.3f kPa\n', out0.pH2/1000);
fprintf('pO2         = %.3f kPa\n', out0.pO2/1000);

fprintf('jLim        = %.4f A/cm^2\n', out0.jLim/1e4);

fprintf('qgen        = %.3e W/m^3\n', out0.qgen);

fprintf('max |dT/dt| = %.3e K/s\n', max(abs(dx0(s.idx.T))));
fprintf('max |dmw/dt|= %.3e kg/(m^3 s)\n', max(abs(dx0(s.idx.mw))));


%% ================================================================
% 11. 检查初始水的物理分配
% ================================================================

fprintf('\nInitial water distribution:\n');

fprintf('max vapor water      = %.6e kg/m^3\n', ...
    max(out0.water.mv));

fprintf('max liquid water     = %.6e kg/m^3\n', ...
    max(out0.water.ml));

fprintf('max CL ionomer water = %.6e kg/m^3\n', ...
    max(out0.water.mIonCL));

fprintf('max pore water       = %.6e kg/m^3\n', ...
    max(out0.water.mPore));


%% ================================================================
% 12. 检查界面水通量
% ================================================================

fprintf('\nWater interface flux:\n');

fprintf('aCL/PEM diffusive = %.6e kg/(m^2 s)\n', ...
    out0.water.interfaceFlux.aCL_PEM_diff);

fprintf('aCL/PEM total     = %.6e kg/(m^2 s)\n', ...
    out0.water.interfaceFlux.aCL_PEM_total);

fprintf('PEM/cCL diffusive = %.6e kg/(m^2 s)\n', ...
    out0.water.interfaceFlux.PEM_cCL_diff);

fprintf('PEM/cCL total     = %.6e kg/(m^2 s)\n', ...
    out0.water.interfaceFlux.PEM_cCL_total);


%% ================================================================
% 13. 建立 ode15s RHS
% ================================================================
%
% ode15s 要求函数形式：
%
%       dx = f(t,x)
%
% 我们自己的函数是：
%
%       dx = pemfc_rhs_mvp(t,x,u,p,g,s)
%
% 所以用匿名函数把其他参数“包起来”。

rhs = @(t,x) pemfc_rhs_mvp(t, x, u, p, g, s);


%% ================================================================
% 14. ODE 求解器设置
% ================================================================

nonNegativeStates = [ ...
    s.idx.H2(:);
    s.idx.O2(:);
    s.idx.mw(:);
    s.idx.mi(:) ];


odeOpt = odeset( ...
    'RelTol',1e-6, ...
    'AbsTol',1e-8, ...
    'MaxStep',0.01, ...
    'NonNegative',nonNegativeStates);


%% ================================================================
% 15. 第一次只跑 0.1 s
% ================================================================

tspan = [0 0.1];


fprintf('\n========================================\n');
fprintf('Starting ode15s: %.3f -> %.3f s\n', ...
    tspan(1), tspan(end));
fprintf('========================================\n');


tic;

[tSol, xSol] = ode15s(rhs, tspan, x0, odeOpt);

cpuTime = toc;


fprintf('\node15s finished successfully.\n');
fprintf('Number of output points = %d\n', numel(tSol));
fprintf('CPU time = %.3f s\n', cpuTime);


%% ================================================================
% 16. 后处理
%
% ode15s 只保存了 x(t)。
%
% 为了获得：
%
%   V(t)
%   Tavg(t)
%   lambda(t)
%   jLim(t)
%
% 我们在每个输出时刻重新调用一次 RHS。
%
% 注意：
% 这里不是重新积分，只是在计算代数量。
% ================================================================

Nt = numel(tSol);

Vcell      = zeros(Nt,1);
Tavg       = zeros(Nt,1);
lambdaMean = zeros(Nt,1);
jLim       = zeros(Nt,1);

minH2 = zeros(Nt,1);
minO2 = zeros(Nt,1);

for k = 1:Nt

    xk = xSol(k,:).';

    [~, outk] = pemfc_rhs_mvp(tSol(k), xk, u, p, g, s);

    Vcell(k)      = outk.Vcell;
    Tavg(k)       = outk.Tavg;
    lambdaMean(k) = outk.lambdaMean;
    jLim(k)       = outk.jLim;

    minH2(k) = outk.minH2;
    minO2(k) = outk.minO2;

end


%% ================================================================
% 17. 最终状态检查
% ================================================================

xEnd = xSol(end,:).';

fprintf('\n========================================\n');
fprintf('Final state at t = %.4f s\n', tSol(end));
fprintf('========================================\n');

fprintf('Vcell       = %.6f V\n', Vcell(end));
fprintf('Tavg        = %.6f degC\n', Tavg(end)-273.15);
fprintf('lambdaMean  = %.6f\n', lambdaMean(end));

fprintf('minimum H2  = %.6e mol/m^3\n', minH2(end));
fprintf('minimum O2  = %.6e mol/m^3\n', minO2(end));

fprintf('minimum mw  = %.6e kg/m^3\n', ...
    min(xEnd(s.idx.mw)));

fprintf('maximum ice = %.6e kg/m^3\n', ...
    max(xEnd(s.idx.mi)));


%% ================================================================
% 18. 最基本的数值健康检查
% ================================================================

assert(all(isfinite(xSol(:))), ...
    'ODE solution contains NaN or Inf.');

assert(all(minH2 >= -1e-10), ...
    'Negative H2 concentration detected.');

assert(all(minO2 >= -1e-10), ...
    'Negative O2 concentration detected.');

assert(all(isfinite(Vcell)), ...
    'Non-finite cell voltage detected.');

fprintf('\nSolution health check PASSED.\n');


%% ================================================================
% 19. 绘图：平均温度
% ================================================================

figure;

plot(tSol, Tavg-273.15, 'LineWidth',1.5);

xlabel('Time / s');
ylabel('Average temperature / ^\circC');

grid on;

title('PEMFC average temperature');


%% ================================================================
% 20. 绘图：电池电压
% ================================================================

figure;

plot(tSol, Vcell, 'LineWidth',1.5);

xlabel('Time / s');
ylabel('Cell voltage / V');

grid on;

title('PEMFC cell voltage');


%% ================================================================
% 21. 绘图：膜平均含水量
% ================================================================

figure;

plot(tSol, lambdaMean, 'LineWidth',1.5);

xlabel('Time / s');
ylabel('\lambda');

grid on;

title('Average membrane water content');


%% ================================================================
% 22. 绘图：最终 H2 浓度分布
% ================================================================

cH2End = xEnd(s.idx.H2);

figure;

plot(g.x(g.idx_H2)*1e6, cH2End, 'o-', 'LineWidth',1.5);

xlabel('x / \mum');
ylabel('H_2 concentration / mol m^{-3}');

grid on;

title('Final H_2 concentration');


%% ================================================================
% 23. 绘图：最终 O2 浓度分布
% ================================================================

cO2End = xEnd(s.idx.O2);

figure;

plot(g.x(g.idx_O2)*1e6, cO2End, 'o-', 'LineWidth',1.5);

xlabel('x / \mum');
ylabel('O_2 concentration / mol m^{-3}');

grid on;

title('Final O_2 concentration');


%% ================================================================
% 24. 绘图：最终总水分布
% ================================================================

mwEnd = xEnd(s.idx.mw);

figure;

plot(g.x*1e6, mwEnd, 'o-', 'LineWidth',1.5);

xlabel('x / \mum');
ylabel('m_w / kg m^{-3}');

grid on;

title('Final total-water distribution');


%% ================================================================
% 25. 完成
% ================================================================

fprintf('\n========================================\n');
fprintf('MVP test finished.\n');
fprintf('========================================\n');