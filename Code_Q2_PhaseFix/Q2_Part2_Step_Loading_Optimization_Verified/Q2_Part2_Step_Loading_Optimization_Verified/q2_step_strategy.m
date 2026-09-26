function st=q2_step_strategy(x)
% Six/seven free levels and ordered switches; computational search bounds.
x=double(x(:).'); d=numel(x); assert(ismember(d,[11 13]));
assert(all(isfinite(x)) && all(x>=0 & x<=1)); n=(d+1)/2;
scale=[29.5 59.5 89.5 119.5 29.5 29.5];
st=struct('type','step','levelsAcm2',.5*x(1:n), ...
    'switchTimesS',cumsum(.5+scale(1:n-1).*x(n+1:end)));
end
