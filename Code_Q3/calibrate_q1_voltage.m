%-------------------------------------------------------------------------------
% @function: 用-20/-25 degC两组数据联合标定问题1电压参数
% @author:   PJ, GPT
% @date:     20260924
% @input:    maxFunctionEvaluations->最大目标函数调用次数，默认40
% @output:   fit->最优参数、粗网格RMSE及搜索记录
%
% [FIT-1] 标定时只改与电压形状直接相关的四个量：
%   j0Ref, fHydDry, ionomerBruggemanExponent, lambdaHydWet。
% Ea和水合曲线其余参数固定，避免小样本下不可辨识。
% [FIT-2] 搜索用粗网格和0.05 s最大步长，不改物理方程；
% 最终参数必须在pemfc_calculate_ice默认细网格上再复核。
%-------------------------------------------------------------------------------

function fit = calibrate_q1_voltage(maxFunctionEvaluations)

    if nargin < 1 || isempty(maxFunctionEvaluations)
        maxFunctionEvaluations = 40;
    end
    if ~isscalar(maxFunctionEvaluations) || ...
            ~isfinite(maxFunctionEvaluations) || ...
            maxFunctionEvaluations < 5
        error('calibrate_q1_voltage:InvalidBudget', ...
            'maxFunctionEvaluations必须是不小于5的有限标量。');
    end

    parameterNames = {'j0Ref','fHydDry', ...
        'ionomerBruggemanExponent','lambdaHydWet'};
    lowerBound = [0.08,0.001,1.20,6.0];
    upperBound = [1.50,0.050,2.30,10.0];
    % 以第一轮搜索+阈值一维扫描的最优点作为复标定起点。
    initialPhysical = [0.537556176898,0.00617369459314, ...
        1.77873441244,8.5];
    initialScaled = physical_to_scaled( ...
        initialPhysical,lowerBound,upperBound);

    baseOption = struct( ...
        'verbose',false, ...
        'nCellLayer',[2 3 4 4 2], ...
        'RelTol',1e-5, ...
        'AbsTol',1e-7, ...
        'MaxStep',0.05, ...
        'Ea',8000, ...
        'lambdaHydOn',3.5, ...
        'nHyd',3.0);

    evaluation = 0;
    history = struct('parameters',{},'score',{}, ...
        'rmse20',{},'rmse25',{});
    objective = @evaluate_candidate;
    searchOption = optimset('Display','iter','MaxFunEvals', ...
        floor(maxFunctionEvaluations),'MaxIter', ...
        floor(maxFunctionEvaluations),'TolX',2e-2,'TolFun',2e-4);
    [bestScaled,bestScore,exitFlag,output] = fminsearch( ...
        objective,initialScaled,searchOption);
    bestPhysical = scaled_to_physical( ...
        bestScaled,lowerBound,upperBound);

    [~,bestIndex] = min([history.score]);
    bestHistory = history(bestIndex);
    fit.parameterNames = parameterNames;
    fit.parameters = bestHistory.parameters;
    fit.score = bestHistory.score;
    fit.rmse20 = bestHistory.rmse20;
    fit.rmse25 = bestHistory.rmse25;
    fit.fminsearchReturnedParameters = bestPhysical;
    fit.fminsearchReturnedScore = bestScore;
    fit.exitFlag = exitFlag;
    fit.output = output;
    fit.fixedParameters = struct('Ea',baseOption.Ea, ...
        'lambdaHydOn',baseOption.lambdaHydOn,'nHyd',baseOption.nHyd);
    fit.searchGrid = baseOption.nCellLayer;
    fit.history = history;
    save('q1_voltage_calibration.mat','fit');

    fprintf('\n粗网格标定完成：score=%.8f V, ',fit.score);
    fprintf('RMSE20=%.8f V, RMSE25=%.8f V\n', ...
        fit.rmse20,fit.rmse25);
    for k = 1:numel(parameterNames)
        fprintf('  %s = %.12g\n',parameterNames{k},fit.parameters(k));
    end

    function score = evaluate_candidate(scaledParameter)

        physical = scaled_to_physical( ...
            scaledParameter,lowerBound,upperBound);
        option = baseOption;
        option.j0Ref = physical(1);
        option.fHydDry = physical(2);
        option.ionomerBruggemanExponent = physical(3);
        option.lambdaHydWet = physical(4);

        evaluation = evaluation+1;
        try
            result20 = pemfc_calculate_ice(-20,option);
            result25 = pemfc_calculate_ice(-25,option);
            rmse20 = result20.voltageRMSE;
            rmse25 = result25.voltageRMSE;
            score = sqrt(0.5*(rmse20^2+rmse25^2));
            if ~result20.health.allFinite || ...
                    ~result25.health.allFinite || ...
                    ~result20.health.integrationComplete || ...
                    ~result25.health.integrationComplete
                score = score+1;
            end
        catch exception
            rmse20 = Inf;
            rmse25 = Inf;
            score = 10;
            fprintf('标定候选求解失败：%s\n',exception.message);
        end

        history(end+1) = struct('parameters',physical, ...
            'score',score,'rmse20',rmse20,'rmse25',rmse25); %#ok<AGROW>
        fprintf(['FIT %02d score=%.7f  rmse=[%.7f %.7f]  ', ...
            'j0=%.5g fDry=%.5g bIon=%.5g lambdaWet=%.5g\n'], ...
            evaluation,score,rmse20,rmse25,physical);
    end

end


function scaled = physical_to_scaled(physical,lowerBound,upperBound)

    fraction = (physical-lowerBound)./(upperBound-lowerBound);
    fraction = min(max(fraction,1e-8),1-1e-8);
    scaled = log(fraction./(1-fraction));

end


function physical = scaled_to_physical(scaled,lowerBound,upperBound)

    fraction = 1./(1+exp(-scaled));
    physical = lowerBound+(upperBound-lowerBound).*fraction;

end
