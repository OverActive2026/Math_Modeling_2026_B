function summary=replay_q2_linear_temperature_checks(temperaturesC)
% Replay the fixed linear curve saved in this standalone paper package.
% The original -10 C optimum's 105% time limit is checked on every run.
if nargin<1 || isempty(temperaturesC)
    temperaturesC=[-10 -11 -12];
end
temperaturesC=double(temperaturesC(:).');
assert(all(isfinite(temperaturesC)) && ...
    all(temperaturesC>=-30 & temperaturesC<=-10), ...
    'Temperatures must be between -30 and -10 C.');
temperaturesC=unique([-10 temperaturesC],'stable');
codeDir=fileparts(mfilename('fullpath'));
addpath(codeDir,'-begin');
assert(strcmp(fileparts(which('pemfc_stack5_simulate')),codeDir));
saved=load(fullfile(codeDir,'optimization_results', ...
    'q2_linear_cap5_verified_report.mat'),'report');
report=saved.report;
assert(report.verified && report.best.success, ...
    'The saved linear curve was not verified.');
referenceTimeS=40.726546388;
maximumWarmTimeS=1.05*referenceTimeS;
assert(abs(report.maximumWarmStartupTimeS-maximumWarmTimeS)<1e-6);
x=report.bestX;
strategy=q2_decode_strategy(x,'linear');
simOpt=report.simOpt;
simOpt.progressIntervalS=5;
outputDir=fullfile(codeDir,'optimization_results','replayed_trajectories');
if ~isfolder(outputDir), mkdir(outputDir); end
n=numel(temperaturesC);
success=false(n,1);
startupTimeS=nan(n,1);
stopTimeS=nan(n,1);
minimumVoltageV=nan(n,1);
chargeCcm2=nan(n,1);
maximumPoreOccupancy=nan(n,1);
maximumIceFraction=nan(n,1);
criticalCell=nan(n,1);
reason=strings(n,1);
for k=1:n
    T0C=temperaturesC(k);
    localOpt=simOpt;
    localOpt.ambientTemperatureC=T0C;
    fprintf('REPLAY T0=%.2f C (%d/%d)\n',T0C,k,n);
    trial=q2_evaluate_real(x,'linear',T0C,localOpt);
    success(k)=trial.success;
    startupTimeS(k)=trial.startup_time;
    stopTimeS(k)=trial.termination_time;
    minimumVoltageV(k)=trial.V_min;
    chargeCcm2(k)=trial.q_use;
    maximumPoreOccupancy(k)=trial.raw.maximumPoreOccupancy;
    maximumIceFraction(k)=trial.ice_peak;
    criticalCell(k)=trial.critical_cell;
    reason(k)=trial.failure_reason;
    temperatureTag=strrep(sprintf('%+.2f',T0C),'.','p');
    save(fullfile(outputDir,['T_' temperatureTag '.mat']), ...
        'T0C','trial','strategy','localOpt','-v7.3');
    fprintf(['REPLAY_RESULT T0=%.2f C success=%d startup=%.6f s ', ...
        'stop=%.6f s Vmin=%.6f V pore=%.6f cell=%g\n'], ...
        T0C,trial.success,trial.startup_time,trial.termination_time, ...
        trial.V_min,trial.raw.maximumPoreOccupancy,trial.critical_cell);
    if T0C==-10
        assert(trial.success && trial.startup_time<=maximumWarmTimeS, ...
            'The -10 C startup did not satisfy the 5%% time limit.');
    end
end
summary=table(temperaturesC(:),success,startupTimeS,stopTimeS, ...
    minimumVoltageV,chargeCcm2,maximumPoreOccupancy, ...
    maximumIceFraction,criticalCell,reason, ...
    'VariableNames',{'temperatureC','success','startupTimeS', ...
    'stopTimeS','minimumVoltageV','chargeCcm2', ...
    'maximumPoreOccupancy','maximumIceFraction', ...
    'criticalCell','reason'});
writetable(summary,fullfile(outputDir,'temperature_checks.csv'));
save(fullfile(outputDir,'replay_summary.mat'), ...
    'summary','strategy','referenceTimeS','maximumWarmTimeS','-v7.3');
end
