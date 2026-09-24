%MAIN_Q2 问题2第一阶段唯一运行入口。
% 当前为 TEST PLACEHOLDER ONLY，不能据此报告启动时间或最优策略。
clear;
clc;
close all;
projectRoot = fileparts(mfilename('fullpath'));
addpath(genpath(projectRoot));
p = parameters_Q2();

% 正式参数仍为 NaN；仅在本次框架测试副本中填入虚拟导热参数。
% TEST PLACEHOLDER ONLY：不覆盖 parameters_Q2.m 中的未知物理量。
pRun = p;
testOnly = true; % 第一阶段入口显式开启；不得用于正式建模结果
if testOnly
    pRun.test.enabled = true;
    pRun.k_eff_stack = p.test.k_eff_stack;
    pRun.delta_stack = p.test.delta_stack;
    warning('Q2:TestOnly', ...
        '正在运行 TEST PLACEHOLDER ONLY 框架检查，结果不是冷启动物理预测。');
end

strategyType = "linear";
strategyPar.j0 = 0.05;  % TEST PLACEHOLDER ONLY [A/cm^2]
strategyPar.r = 0.005;  % TEST PLACEHOLDER ONLY [A/(cm^2 s)]
strategyPar.j_cap = 0.30; % TEST PLACEHOLDER ONLY [A/cm^2]

result = simulate_strategy(strategyType,strategyPar,pRun);
fprintf('状态=%s；结束时间=%.2f s；累计电荷=%.4f C/cm^2\n', ...
    result.failure_reason,result.time(end),result.q_use);
disp(result);
plot_results(result,pRun);
