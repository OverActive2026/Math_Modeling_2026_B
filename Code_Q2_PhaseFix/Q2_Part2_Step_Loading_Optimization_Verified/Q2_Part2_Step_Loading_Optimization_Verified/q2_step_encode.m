function x=q2_step_encode(st)
n=numel(st.levelsAcm2); assert(ismember(n,[6 7]));
assert(numel(st.switchTimesS)==n-1);
scale=[29.5 59.5 89.5 119.5 29.5 29.5];
x=[st.levelsAcm2/.5,(diff([0 st.switchTimesS])-.5)./scale(1:n-1)];
assert(all(isfinite(x)) && all(x>=0 & x<=1),'Q2:OutsideSearchBox', ...
    'Strategy outside the documented computational search bounds.');
end
