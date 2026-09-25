function scan = find_q2_optimized_min_temperature(resolutionC)
% Find the coldest successful initial temperature for the *fixed* optimized
% 2 s ramp-then-hold loading, using the same coarse-grid model and limits.
if nargin<1 || isempty(resolutionC), resolutionC = 0.1; end
assert(isscalar(resolutionC) && isfinite(resolutionC) && resolutionC>0);
codeDir = fileparts(mfilename('fullpath'));
addpath(codeDir,'-begin');
reportDir = fullfile(codeDir,'optimization_results', ...
    'constant_fixed_reference_kfreeze03');
files = dir(fullfile(reportDir,'*_report.mat'));
assert(numel(files)==1,'Expected one verified optimization report.');
saved = load(fullfile(files(1).folder,files(1).name),'report');
report = saved.report;
assert(report.verified && report.verifiedResult.success);
strategy = report.verifiedResult.strategy;
baseOpt = report.verifiedResult.options;
assert(strcmp(strategy.type,'linear') && ...
    abs(strategy.initialAcm2-0.02)<1e-12 && ...
    abs(strategy.plateauAcm2-0.2668808447608863)<1e-12);
assert(baseOpt.gammaIce==3.5 && baseOpt.kFreeze==0.3 && ...
    isequal(baseOpt.nCellLayer,[2 3 4 4 2]) && ...
    baseOpt.qMaxCcm2==20 && baseOpt.jMaxAcm2==0.5);
outDir = fullfile(codeDir,'optimization_results','fixed_optimized_temperature');
if ~exist(outDir,'dir'), mkdir(outDir); end
scanFile = fullfile(outDir,'scan.mat');
scan = struct('strategy',strategy,'baseOptions',baseOpt, ...
    'reportFile',fullfile(files(1).folder,files(1).name), ...
    'tests',struct('temperatureC',{},'success',{},'successTimeS',{}, ...
    'stopTimeS',{},'stopReason',{},'minimumCellVoltageV',{}, ...
    'minimumVoltageCellIndex',{},'maximumIceVolumeFraction',{}, ...
    'maximumIceCellIndex',{},'chargeUsedCcm2',{},'file',{}), ...
    'coldestSuccessC',NaN,'nearestFailureC',NaN, ...
    'resolutionC',resolutionC,'complete',false);
reference = report.verifiedResult;
if exist(scanFile,'file')
    previous = load(scanFile,'scan');
    assert(isequaln(previous.scan.strategy,strategy) && ...
        isequaln(previous.scan.baseOptions,baseOpt), ...
        'Saved temperature scan used different model settings.');
    scan = previous.scan;
    scan.resolutionC = resolutionC;
    % The copied scan may contain absolute paths from its source folder.
    scan.reportFile = fullfile(files(1).folder,files(1).name);
    for k = 1:numel(scan.tests)
        if abs(scan.tests(k).temperatureC+10)<1e-12
            scan.tests(k).file = scan.reportFile;
        else
            scan.tests(k).file = fullfile(outDir, ...
                sprintf('T_%+09.5f.mat',scan.tests(k).temperatureC));
        end
        assert(exist(scan.tests(k).file,'file')==2, ...
            'Missing saved trajectory: %s',scan.tests(k).file);
    end
    warm = scan.coldestSuccessC;
    cold = scan.nearestFailureC;
    fprintf('RESUME failed=%.6f C successful=%.6f C\n',cold,warm);
else
    scan.tests(end+1) = summarize(-10,reference,'verified optimization report');
    fprintf('KNOWN T=-10.000 C SUCCESS startup=%.6f s\n',reference.successTimeS);
    warm = -10;
    cold = NaN;
    for temperatureC = -11:-1:-30
        outcome = evaluate(temperatureC);
        if ~outcome.success
            cold = temperatureC;
            break
        end
        warm = temperatureC;
    end
end
if isnan(cold)
    scan.coldestSuccessC = warm;
    save(scanFile,'scan');
    fprintf('NO_FAILURE_FOUND: success through %.3f C; lower boundary unknown.\n',warm);
    return
end
while warm-cold > scan.resolutionC+1e-9
    temperatureC = (warm+cold)/2;
    outcome = evaluate(temperatureC);
    if outcome.success
        warm = temperatureC;
    else
        cold = temperatureC;
    end
end
scan.coldestSuccessC = warm;
scan.nearestFailureC = cold;
scan.complete = true;
save(scanFile,'scan');
fprintf('BOUNDARY failed=%.6f C successful=%.6f C width=%.6f C\n', ...
    cold,warm,warm-cold);

    function outcome = evaluate(temperatureC)
        opt = baseOpt;
        opt.ambientTemperatureC = temperatureC;
        opt.plot = false;
        opt.verbose = false;
        opt.progressIntervalS = 10;
        opt.progressWallIntervalS = 30;
        opt.maxWallTimeS = 300;
        fprintf('START T=%.6f C\n',temperatureC);
        [result,~] = pemfc_stack5_simulate(strategy,temperatureC,opt);
        file = fullfile(outDir,sprintf('T_%+09.5f.mat',temperatureC));
        save(file,'result','strategy','opt');
        outcome = summarize(temperatureC,result,file);
        scan.tests(end+1) = outcome;
        save(scanFile,'scan');
        fprintf('RESULT T=%.6f success=%d start=%.6f stop=%.6f ', ...
            temperatureC,outcome.success,outcome.successTimeS,outcome.stopTimeS);
        fprintf('Vmin=%.6f cell=%d ice=%.6f iceCell=%d charge=%.6f reason=%s\n', ...
            outcome.minimumCellVoltageV,outcome.minimumVoltageCellIndex, ...
            outcome.maximumIceVolumeFraction,outcome.maximumIceCellIndex, ...
            outcome.chargeUsedCcm2,outcome.stopReason);
    end
end

function row = summarize(temperatureC,result,file)
row = struct('temperatureC',temperatureC,'success',result.success, ...
    'successTimeS',result.successTimeS,'stopTimeS',result.stopTimeS, ...
    'stopReason',result.stopReason, ...
    'minimumCellVoltageV',result.minimumCellVoltageV, ...
    'minimumVoltageCellIndex',result.minimumVoltageCellIndex, ...
    'maximumIceVolumeFraction',result.maximumIceVolumeFraction, ...
    'maximumIceCellIndex',result.maximumIceCellIndex, ...
    'chargeUsedCcm2',result.chargeUsedCcm2,'file',file);
end
