%-------------------------------------------------------------------------------
% 问题4预冷温度场与冷启动初值传递测试
%-------------------------------------------------------------------------------

clearvars;
clc;

codeDir = fileparts(mfilename('fullpath'));
addpath(codeDir);

opt.nCellLayer = [2 3 4 4 2];
opt.initialTemperatureC = 25;
opt.ambientTemperatureC = -30;
opt.outputStepS = 0.25;
opt.MaxStep = 0.05;
opt.RelTol = 1e-8;
opt.AbsTol = 1e-10;
opt.plot = false;
opt.verbose = false;

[pre,preModel] = q4_precooling_temperature_field(0.05,opt);
assert(isequal(size(pre.initialState.temperatureCellK), ...
    [preModel.g.N,5]));
assert(all(pre.cellAverageTemperatureC(end,:) < ...
    opt.initialTemperatureC));
assert(all(pre.cellAverageTemperatureC(end,:) > ...
    opt.ambientTemperatureC));
assert(pre.symmetryErrorK < 1e-7);
assert(pre.energyBalanceRelativeError < 5e-5);

startupOpt.nCellLayer = opt.nCellLayer;
startupOpt.initialState = pre.initialState;
startupOpt.ambientTemperatureC = -30;
startupOpt.tMax = 0.01;
startupOpt.outputStep = 0.01;
startupOpt.MaxStep = 0.005;
startupOpt.RelTol = 1e-5;
startupOpt.AbsTol = 1e-8;
startupOpt.plot = false;
startupOpt.verbose = false;
startupOpt.useJacobianPattern = true;

[~,startupModel] = pemfc_stack5_simulate_q4( ...
    'preheat',zeros(1,5),0,startupOpt);
for k = 1:5
    temperatureIndex = startupModel.cellStateIndex{k}( ...
        startupModel.s.idx.T);
    transferredTemperatureK = startupModel.x0(temperatureIndex);
    assert(max(abs(transferredTemperatureK- ...
        pre.initialState.temperatureCellK(:,k))) < 1e-12);
end
assert(max(abs(startupModel.x0(startupModel.endPlateIndex)- ...
    pre.initialState.endPlateTemperatureK)) < 1e-12);

fprintf('test_q4_precooling: PASS\n');
