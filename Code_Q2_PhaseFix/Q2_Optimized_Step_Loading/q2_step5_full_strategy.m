function strategy=q2_step5_full_strategy(x)
%Q2_STEP5_FULL_STRATEGY All five levels and four switches are free.
% Levels use the physical 0--0.5 A/cm^2 bound. Switch-time windows are
% broad but nonoverlapping so every normalized coordinate remains active.
% If the optimum touches a time-window boundary, expand the window and
% rerun; this is a search box, not a physical restriction from the paper.
x=double(x(:).');
assert(numel(x)==9 && all(isfinite(x)) && all(x>=0) && all(x<=1));
timeLower=[1 8.1 17 24];
timeUpper=[8 16.9 23.9 38];
strategy=struct('type','step','levelsAcm2',0.5*x(1:5), ...
    'switchTimesS',timeLower+(timeUpper-timeLower).*x(6:9));
end
