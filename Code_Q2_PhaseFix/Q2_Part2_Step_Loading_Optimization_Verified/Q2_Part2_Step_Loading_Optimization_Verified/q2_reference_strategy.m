function st=q2_reference_strategy(stage)
if nargin<1, stage=7; end
assert(ismember(stage,[6 7]));
if stage==7
    st=struct('type','step', ...
        'levelsAcm2',[.175184 .291442 .242683 .261777 .5 .40 .43], ...
        'switchTimesS',[2.5713 12.1507 20.5979 27.5905 31 32.5]);
else
    a=load(fullfile(fileparts(mfilename('fullpath')),'reference_six.mat'),'strategy');
    st=a.strategy; % Exact audited six-stage curve, not rounded chat values.
end
end
