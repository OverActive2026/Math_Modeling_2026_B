function strategy = q2_decode_parameters(x,cfg)
% 单位立方体 -> 合法、单调升载的策略；不是直接把归一化量当电流。
x = x(:).';
assert(all(isfinite(x)) && all(x>=0 & x<=1),'x must lie in [0,1].');
assert(cfg.currentMax>0 && cfg.currentMax<=0.5);
switch lower(cfg.strategyType)
    case 'constant'
        assert(isscalar(x));
        assert(cfg.constantInitialCurrentAcm2>=0 && ...
            cfg.constantInitialCurrentAcm2<=cfg.currentMax);
        assert(cfg.constantRampTimeS>0);
        bounds = cfg.constantPlateauBounds;
        assert(numel(bounds)==2 && bounds(1)>=cfg.constantInitialCurrentAcm2 && ...
            bounds(2)<=cfg.currentMax && bounds(2)>bounds(1));
        plateau = bounds(1)+diff(bounds)*x(1);
        % 与run_q2_constant_current完全相同：短时线性升载后恒流。
        strategy = struct('type','linear', ...
            'initialAcm2',cfg.constantInitialCurrentAcm2, ...
            'rampRateAcm2s', ...
            (plateau-cfg.constantInitialCurrentAcm2)/cfg.constantRampTimeS, ...
            'plateauAcm2',plateau);
    case 'linear'
        assert(numel(x)==3);
        assert(numel(cfg.rampBounds)==2 && cfg.rampBounds(1)>0 && ...
            cfg.rampBounds(2)>=cfg.rampBounds(1));
        plateau = cfg.currentMax*x(3);
        strategy = struct('type','linear','initialAcm2',plateau*x(1), ...
            'rampRateAcm2s',cfg.rampBounds(1)+diff(cfg.rampBounds)*x(2), ...
            'plateauAcm2',plateau);
    case 'step'
        nLevels = 3;
        if isfield(cfg,'stepStageCount'), nLevels = cfg.stepStageCount; end
        assert(nLevels>=2 && fix(nLevels)==nLevels);
        assert(numel(x)==2*nLevels-1);
        nSwitches = nLevels-1;
        assert(cfg.switchMinS>0 && cfg.switchGapS>0 && ...
            cfg.switchMaxS>cfg.switchMinS+(nSwitches-1)*cfg.switchGapS && ...
            cfg.switchMaxS<cfg.simOpt.tMax);
        levels = cfg.currentMax*sort(x(1:nLevels));
        u = sort(x(nLevels+1:end));
        available = cfg.switchMaxS-cfg.switchMinS-(nSwitches-1)*cfg.switchGapS;
        switchTimes = cfg.switchMinS + (0:nSwitches-1)*cfg.switchGapS + available*u;
        strategy = struct('type','step','levelsAcm2',levels, ...
            'switchTimesS',switchTimes);
    otherwise
        error('q2:UnknownStrategy','Unknown strategy type.');
end
end
