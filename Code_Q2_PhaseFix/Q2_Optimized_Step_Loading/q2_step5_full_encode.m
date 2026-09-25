function x=q2_step5_full_encode(strategy)
%Q2_STEP5_FULL_ENCODE Inverse of q2_step5_full_strategy for warm starts.
assert(strcmpi(char(strategy.type),'step'));
j=double(strategy.levelsAcm2(:).');
t=double(strategy.switchTimesS(:).');
assert(numel(j)==5 && numel(t)==4);
timeLower=[1 8.1 17 24];
timeUpper=[8 16.9 23.9 38];
x=[j/0.5,(t-timeLower)./(timeUpper-timeLower)];
assert(all(x>=0) && all(x<=1), ...
    'Prior five-step strategy lies outside the full-search box.');
end
