function strategy = q2_decode_strategy(x,kind)
%Q2_DECODE_STRATEGY Map a normalized search point to a valid current curve.
% x is in [0,1]^d; current density is in A/cm^2 and times are in s.
x = double(x(:).');
if any(~isfinite(x)) || any(x<0) || any(x>1)
    error('Q2:InvalidPoint','All normalized parameters must be in [0,1].');
end
switch lower(char(kind))
    case 'constant'
        % Matches run_q2_constant_current.m: a ramp followed by a plateau.
        assert(numel(x)==3,'Ramp-hold strategy needs three parameters.');
        % The problem permits 0 <= j <= 0.50 A/cm^2.  The former 0.10
        % upper bound on the initial current was a search-design choice.
        j0 = 0.50*x(1);
        jp = j0+(0.50-j0)*x(2);
        rampTime = 0.5+19.5*x(3);
        strategy = struct('type','linear','initialAcm2',j0, ...
            'rampRateAcm2s',(jp-j0)/rampTime,'plateauAcm2',jp);
    case 'strict_constant'
        assert(numel(x)==1,'Strict constant strategy needs one parameter.');
        strategy = struct('type','constant','currentAcm2',0.02+0.48*x(1));
    case 'linear'
        assert(numel(x)==3,'Linear strategy needs three parameters.');
        j0 = 0.50*x(1);
        strategy = struct('type','linear', ...
            'initialAcm2',j0, ...
            'rampRateAcm2s',0.001+0.049*x(2), ...
            'plateauAcm2',j0+(0.50-j0)*x(3));
    case 'step'
        % Matches run_q2_step_current.m: four independent levels and
        % three strictly increasing switch times.
        assert(numel(x)==7,'Step strategy needs seven parameters.');
        levels = 0.02+0.48*x(1:4);
        t1 = 2+28*x(5);
        t2 = t1+2+(60-t1-2)*x(6);
        t3 = t2+2+(90-t2-2)*x(7);
        strategy = struct('type','step', ...
            'levelsAcm2',levels, ...
            'switchTimesS',[t1,t2,t3]);
    otherwise
        error('Q2:UnknownStrategy','Unknown strategy type: %s',char(kind));
end
end
