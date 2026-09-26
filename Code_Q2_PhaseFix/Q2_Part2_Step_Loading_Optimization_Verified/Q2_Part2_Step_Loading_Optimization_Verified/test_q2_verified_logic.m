function test_q2_verified_logic()
% Synthetic regression only. No mock startup number is a physical result.
folder=fileparts(mfilename('fullpath')); addpath(folder,'-begin'); rng(37);
for d=[11 13]
    for j=1:100
        x=rand(1,d); assert(max(abs(q2_step_encode(q2_step_strategy(x))-x))<1e-12);
    end
    for j=1:d
        x=.5*ones(1,d); a=q2_step_strategy(x); x(j)=.6;
        assert(~isequaln(a,q2_step_strategy(x)));
    end
end
% Later-than-warm-start switches are allowed, not falsely declared infeasible.
s=q2_reference_strategy(7); s.switchTimesS(end)=50;
assert(max(abs(q2_step_strategy(q2_step_encode(s)).switchTimesS-s.switchTimesS))<1e-12);
opt=q2_verified_options(); mode=1; seen=[];
st=q2_reference_strategy(7);
r=runMock(); assert(r.success && all(r.status==1) && r.startup_time==34);
assert(isequal(seen,[-10 -11.2]));
mode=2; r=runMock(); assert(~r.success && isnan(r.startup_time) && r.status(2)==-1);
mode=3; r=runMock(); assert(~r.success && r.time_censored && r.status(1)==-3);
mode=4; r=runMock(); assert(~r.success && r.time_censored && r.status(2)==-3);
mode=5; r=runMock(); assert(~r.success && r.solver_failed && r.status(2)==-2);
mode=6; r=runMock(); assert(~r.success && r.solver_failed && r.status(1)==-2);
mode=7; r=runMock(); assert(~r.success && r.time_censored && isnan(r.startup_time));
mode=1; overlay=opt; overlay.tMax=35; overlay.optimizationCutoffActive=true;
overlaySeen=[];
r=q2_evaluate_schedule(st,[-10 -11.2],overlay,600,@horizonMock);
assert(r.success && isequal(overlaySeen,[35 600]));
% Optimizer integration: custom cutoff physics carried through, censor excluded.
calls=0; mode=1;
o=struct('maxEvaluations',3,'initialSamples',3, ...
    'initialPoints',[.4*ones(1,13);.5*ones(1,13);.6*ones(1,13)], ...
    'simOpt',opt,'evaluator',@optMock,'useIncumbentCutoff',true);
report=q2_optimize_trbo('step7_full',o);
assert(report.nTimeCensored==1 && report.nSuccess==2);
assert(isnan(report.signedConstraint(3)) && isnan(report.startupTimeS(3)));
bad=o; bad.simOpt=struct(); rejected=false;
try, q2_optimize_trbo('step7_full',bad); catch ME, rejected=strcmp(ME.identifier,'Q2:IncompletePhysics'); end
assert(rejected);
fprintf('REGRESSION PASSED: roundtrips, all variables, joint failure, numerical unknown, both horizons, censor exclusion, incomplete overlay. SYNTHETIC ONLY.\n');
    function r=runMock()
        r=q2_evaluate_schedule(st,[-10 -11.2],opt,38,@mock);
    end
    function [r,m]=horizonMock(st,T,o)
        overlaySeen(end+1)=o.tMax;
        if T==-11.2, assert(~isfield(o,'optimizationCutoffActive')); end
        [r,m]=mock(st,T,o);
    end
    function r=optMock(x,o)
        if nargin<2, o=opt; end
        calls=calls+1;
        assert(o.gammaIce==3.5 && o.kFreeze==.4 && isequal(o.nCellLayer,[2 3 4 4 2]));
        if calls>=2, assert(o.tMax==34.5 && o.coldTMaxS==600); end
        mode=1; if calls==3,mode=3;end
        r=q2_evaluate_schedule(q2_step_strategy(x),[-10 -11.2],o,600,@mock);
    end
    function [r,m]=mock(~,T,o)
        seen(end+1)=T;
        if (mode==5 && T==-11.2) || mode==6, error('Q2:MockNumeric','Synthetic numerical failure'); end
        okay=true; event=1; t=34; if T==-11.2,t=50;end
        if mode==2 && T==-11.2, okay=false; event=3; end
        if (mode==3 && T==-10) || (mode==4 && T==-11.2)
            okay=false;event=[];
        end
        if mode==7 && T==-10,t=39;end
        if ~okay,t=NaN;end
        r=struct('solverTerminatedUnexpectedly',false,'success',okay, ...
            'successTimeS',t,'minimumCellVoltageV',.35,'chargeUsedCcm2',15, ...
            'maximumCurrentDensityAcm2',.4,'maximumIceVolumeFraction',.2, ...
            'maxIceVolumeFractionCell',.2*ones(1,5),'maximumPoreOccupancy',.4, ...
            'TavgC',zeros(1,5),'minimumVoltageCellIndex',1,'stopTimeS',o.tMax, ...
            'event',struct('index',event),'stopReason','synthetic test');m=[];
    end
end
