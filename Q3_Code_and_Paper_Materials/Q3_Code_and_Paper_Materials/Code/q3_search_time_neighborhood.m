function logTable = q3_search_time_neighborhood()
% Check whether less than maximum power can beat the verified time bound.
% Uses the same Q3 model and constraints; each probe is a physical solve.
    codeDir = fileparts(mfilename('fullpath'));
    addpath(codeDir);
    outDir = fullfile(codeDir,'results');
    archive = load(fullfile(outDir,'time_priority_candidates.mat'), ...
        'pre','co','opt');
    opt = archive.opt;
    opt.plot = false;
    opt.verbose = false;
    opt.enforceChargeCap = true;
    baseQ = ones(1,5);
    step = 0.05;
    cutoffMarginS = 0.001;
    modeList = {'preheat','coheat'};
    rows = struct('mode',{},'changedCell',{},'q',{}, ...
        'plannedHeatingS',{},'searchCutoffS',{},'status',{}, ...
        'successTimeS',{},'minimumTemperatureC',{}, ...
        'minimumVoltageV',{},'maximumIceVolumeFraction',{}, ...
        'chargeCcm2',{},'elapsedS',{});
    n = 0;
    for m = 1:2
        mode = modeList{m};
        if strcmp(mode,'preheat')
            currentBest = archive.pre;
        else
            currentBest = archive.co;
        end
        for cellIndex = 1:5
            q = baseQ;
            q(cellIndex) = 1-step;
            cutoff = currentBest.successTimeS-cutoffMarginS;
            trialOpt = opt;
            trialOpt.tMax = cutoff;
            clock = tic;
            status = "numerical_failure";
            ts = NaN;
            minT = NaN;
            minV = NaN;
            maxIce = NaN;
            charge = NaN;
            try
                r = pemfc_stack5_simulate_q3(mode,q,cutoff,trialOpt);
                minT = min(r.TavgC(end,:));
                minV = r.minimumCellVoltageV;
                maxIce = r.maximumIceVolumeFraction;
                charge = r.chargeUsedCcm2;
                if ~r.solverTerminatedUnexpectedly
                    if r.success
                        status = "faster_success";
                        ts = r.successTimeS;
                        if strcmp(mode,'preheat')
                            archive.pre = r;
                        else
                            archive.co = r;
                        end
                        currentBest = r;
                    else
                        status = "not_faster";
                    end
                end
            catch ME
                status = string(ME.identifier);
            end
            n = n+1;
            rows(n) = struct('mode',mode,'changedCell',cellIndex, ...
                'q',q,'plannedHeatingS',cutoff, ...
                'searchCutoffS',cutoff,'status',char(status), ...
                'successTimeS',ts,'minimumTemperatureC',minT, ...
                'minimumVoltageV',minV, ...
                'maximumIceVolumeFraction',maxIce, ...
                'chargeCcm2',charge,'elapsedS',toc(clock));
            fprintf('%s cell%d q=%.2f: %s, minT=%.4f C, %.1f s\n', ...
                mode,cellIndex,q(cellIndex),status,minT,rows(n).elapsedS);
            write_checkpoint();
        end
    end
    logTable = flatten_rows(rows);
    save(fullfile(outDir,'time_neighborhood_search.mat'), ...
        'logTable','archive','step','cutoffMarginS','-v7.3');

    function write_checkpoint()
        partial = flatten_rows(rows);
        writetable(partial,fullfile(outDir,'time_neighborhood_search.csv'), ...
            'Encoding','UTF-8');
    end
end

function t = flatten_rows(rows)
    t = struct2table(rows);
    q = vertcat(rows.q);
    t.q = [];
    for k = 1:5
        t.(sprintf('q%d',k)) = q(:,k);
    end
end
