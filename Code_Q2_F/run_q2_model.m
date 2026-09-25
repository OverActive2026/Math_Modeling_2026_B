%-------------------------------------------------------------------------------
% Code_Q2_F —— 问题2 总入口
%
% 电池模型：Code_ICE_F 的单电池物理内核（问题1 的标定结果）。
% 电堆层：五片热耦合 + 端板独立热节点 + 累积电荷 + 题目约束事件。
%
% 本脚本做两件事：
%   1. 用同一组初始温度把三种加载策略各跑一次，打印表3 需要的六项指标；
%   2. 给出每种策略的一维/二维粗扫描（可选），用于找"最优加载参数"。
%
% 三种工况（题目图3）：
%   恒流      j(t) = j_c                        见 run_q2_constant_current.m
%   线性升载  j(t) = min(r*t, j_p)              见 run_q2_linear_current.m
%   分段阶梯  j(t) = j1 / j2 / j3               见 run_q2_step_current.m
%
% 用法：
%   run_q2_model                 % 只跑三种策略的基准工况
%   run_q2_model('scan')         % 额外做粗扫描（耗时较长）
%-------------------------------------------------------------------------------

function run_q2_model(mode)

    if nargin < 1 || isempty(mode)
        mode = 'baseline';
    end

    codeDir = fileparts(mfilename('fullpath'));
    addpath(codeDir);

    %% 1. 公共设置

    T0 = -10;                 % 题目2(1) 固定初始温度 [degC]
    showPlot = true;
    tMax = 300;               % 最长积分时间 [s]

    simOpt = struct('tMax',tMax,'verbose',true,'plot',false, ...
        'outputStep',0.5,'MaxStep',0.05, ...
        'progressIntervalS',1);   % 每 1 s 仿真时间打印一次求解进度，0=关闭

    % 三种策略的已标定最优参数（粗网格扫描给出，见 README_Q2F.md 第3节）
    % [Q2F-13] 恒流必须带软启动：t=0 直接跳到 j_c 时 cCL 离聚物仍是干态，
    % j_c>=0.16 A/cm^2 会立刻把电压压到 0.30 V 以下；纯阶跃在整个
    % j_c = 0.16~0.50 区间全部失败。
    strategies = struct();
    strategies.constant = struct('type','constant','currentAcm2',0.26, ...
        'rampTimeS',5);
    strategies.linear = struct('type','linear','rampRateAcm2s',0.05, ...
        'plateauAcm2',0.50);
    strategies.step = struct('type','step', ...
        'levelsAcm2',[0.12 0.35 0.50],'switchTimesS',[10 20]);

    names = {'constant','linear','step'};
    labels = {'恒流','线性升载','分段阶梯'};

    fprintf('\n########################################################\n');
    fprintf('# 问题2：5片电堆自冷启动，T0 = %g degC\n',T0);
    fprintf('# 约束：q_max=20 C/cm^2, j_max=0.5 A/cm^2, V_min=0.30 V,\n');
    fprintf('#       eps_ice<0.99, 启动=五片平均温度全部达到0 degC\n');
    fprintf('########################################################\n');

    results = struct();
    for k = 1:numel(names)
        fprintf('\n---------- %s ----------\n',labels{k});
        results.(names{k}) = pemfc_stack5_simulate( ...
            strategies.(names{k}),T0,simOpt);
        print_table3_row(labels{k},results.(names{k}));
    end

    %% 2. 表3 汇总

    fprintf('\n========== 表3 不同加载策略下电堆冷启动性能对比 ==========\n');
    fprintf('%-10s | %-28s | %9s | %11s | %11s | %9s | %11s | %s\n', ...
        '加载策略','最优加载参数','启动时间/s','累计电荷/Ccm-2', ...
        '最大电流/Acm-2','最低电压/V','最大冰体积分','启动结果');
    for k = 1:numel(names)
        r = results.(names{k});
        fprintf('%-10s | %-28s | %9.2f | %11.4f | %11.4f | %9.5f | %11.4g | %s\n', ...
            labels{k},strategy_text(r.strategy),r.successTimeS, ...
            r.chargeUsedCcm2,r.maximumCurrentDensityAcm2, ...
            r.minimumCellVoltageV,r.maximumIceVolumeFraction, ...
            success_text(r.success));
    end

    %% 3. 可选：粗扫描找最优加载参数

    if strcmpi(mode,'scan')
        scanfcn = @(s,t0,o) pemfc_stack5_simulate(s,t0,o);
        scanOpt = simOpt;
        scanOpt.verbose = false;
        scanOpt.progressIntervalS = 0;
        scanOpt.tMax = 240;
        scanOpt.outputStep = 1.0;

        fprintf(['\n========== 恒流+软启动 二维扫描 ==========\n' ...
            't_ramp=0 即纯阶跃恒流，用于对照"一上来就失败"。\n' ...
            '注意 j_c*min(t/t_ramp,1) 与线性升载 min(r*t,j_p) 同族（r=j_c/t_ramp），\n' ...
            '故此处只对小斜坡（t_ramp<=5 s）寻优，以保持"恒流"的物理含义。\n']);
        fprintf('%8s %8s | %9s | %11s | %9s | %s\n', ...
            'j_c','t_ramp','t_s/s','q_use','V_min','结果');
        best = struct('score',inf);
        for tr = [0 4 5]
            for jc = 0.16:0.02:0.40
                s = struct('type','constant','currentAcm2',jc, ...
                    'rampTimeS',tr);
                r = scanfcn(s,T0,scanOpt);
                score = scan_score(r);
                fprintf('%8.3f %8.1f | %9.2f | %11.4f | %9.5f | %s\n', ...
                    jc,tr,r.successTimeS,r.chargeUsedCcm2, ...
                    r.minimumCellVoltageV,success_text(r.success));
                if score < best.score
                    best = struct('score',score,'strategy',s,'result',r);
                end
            end
        end
        fprintf('恒流(+短软启动)最优：j_c=%.3f A/cm^2，t_ramp=%g s，t_s=%.2f s\n', ...
            best.strategy.currentAcm2,best.strategy.rampTimeS, ...
            best.result.successTimeS);

        fprintf('\n========== 线性升载二维粗扫描 ==========\n');
        fprintf('%9s %9s | %9s | %11s | %9s | %s\n', ...
            'r','j_p','t_s/s','q_use','V_min','结果');
        bestLinear = struct('score',inf);
        for rr = [0.02 0.03 0.04 0.05 0.06]
            for jp = [0.35 0.40 0.45 0.50]
                s = struct('type','linear','rampRateAcm2s',rr, ...
                    'plateauAcm2',jp);
                r = scanfcn(s,T0,scanOpt);
                score = scan_score(r);
                fprintf('%9.4f %9.3f | %9.2f | %11.4f | %9.5f | %s\n', ...
                    rr,jp,r.successTimeS,r.chargeUsedCcm2, ...
                    r.minimumCellVoltageV,success_text(r.success));
                if score < bestLinear.score
                    bestLinear = struct('score',score,'strategy',s,'result',r);
                end
            end
        end
        fprintf('线性最优：r=%.4f A/(cm^2 s)，j_p=%.3f A/cm^2，t_s=%.2f s\n', ...
            bestLinear.strategy.rampRateAcm2s, ...
            bestLinear.strategy.plateauAcm2,bestLinear.result.successTimeS);
    end

    if showPlot
        for k = 1:numel(names)
            r = results.(names{k});
            plot_saved = r.simOpt.plot;
            if ~plot_saved
                plot_strategy_result(r);
            end
        end
    end

end


%% ========================================================================
% 辅助函数
% ========================================================================

function score = scan_score(r)

    % 可行解按启动时间排序；不可行解按违反约束的严重程度排在后面。
    if r.success
        score = r.successTimeS;
    else
        score = 1e6+1000*max(0,-r.margin.minimumVoltageV) ...
            +1000*max(0,-r.margin.qMaxCcm2)+1000*max(0,-r.margin.temperatureC);
    end

end


function text = strategy_text(strategy)

    switch lower(char(strategy.type))
        case 'constant'
            if isfield(strategy,'rampTimeS') && ~isempty(strategy.rampTimeS) ...
                    && strategy.rampTimeS > 0
                text = sprintf('j_c=%.3f A/cm^2, 软启动%.3g s', ...
                    strategy.currentAcm2,strategy.rampTimeS);
            else
                text = sprintf('j_c=%.3f A/cm^2',strategy.currentAcm2);
            end
        case 'linear'
            text = sprintf('r=%.4f, j_p=%.3f',strategy.rampRateAcm2s, ...
                strategy.plateauAcm2);
        case 'step'
            text = sprintf('[%s] @ [%s]s', ...
                num2str(strategy.levelsAcm2(:).'), ...
                num2str(strategy.switchTimesS(:).'));
        otherwise
            text = char(strategy.type);
    end

end


function text = success_text(success)

    if success
        text = '成功';
    else
        text = '失败';
    end

end


function print_table3_row(label,r)

    fprintf(['%s: t_s=%.3f s, q_use=%.4f C/cm^2, j_max=%.4f A/cm^2, ', ...
        'V_min=%.5f V, eps_ice,max=%.4g, success=%d\n'], ...
        label,r.successTimeS,r.chargeUsedCcm2, ...
        r.maximumCurrentDensityAcm2,r.minimumCellVoltageV, ...
        r.maximumIceVolumeFraction,r.success);

end


function plot_strategy_result(result)

    figure('Name',sprintf('Q2 %s',char(result.strategy.type)));
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

    nexttile;
    plot(result.t,result.currentDensityAcm2,'LineWidth',1.5);
    xlabel('Time / s'); ylabel('j / (A cm^{-2})');
    title('Current loading'); grid on;

    nexttile;
    plot(result.t,result.TavgC,'LineWidth',1.2); hold on;
    yline(0,'k--');
    xlabel('Time / s'); ylabel('T_{avg} / ^\circC');
    title('Five cell temperatures'); grid on;

    nexttile;
    plot(result.t,result.Vcell,'LineWidth',1.2); hold on;
    yline(0.30,'k--');
    xlabel('Time / s'); ylabel('V / V');
    title('Five cell voltages'); grid on;

    nexttile;
    plot(result.t,result.iceVolumeFractionMax,'LineWidth',1.2);
    xlabel('Time / s'); ylabel('\epsilon_{ice}');
    title('Maximum ice volume fraction'); grid on;

end
