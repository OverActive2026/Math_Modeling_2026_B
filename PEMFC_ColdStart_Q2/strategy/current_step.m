function j = current_step(t,par)
%CURRENT_STEP 三段阶梯加载；t、t1、t2 [s]，电流密度 [A/cm^2]。
p = parameters_Q2();
if isfield(par,'j_max'), limit = par.j_max; else, limit = p.j_max; end
required = {'j1','j2','j3','t1','t2'};
if ~all(isfield(par,required))
    error('Q2:InvalidStepStrategy','缺少 j1、j2、j3、t1 或 t2。');
end
v = [par.j1,par.j2,par.j3,par.t1,par.t2];
if ~isscalar(t) || ~isfinite(t) || t < 0 || ...
        ~all(isfinite(v)) || par.j1 < 0 || ...
        par.j1 > par.j2 || par.j2 > par.j3 || ...
        par.j3 > limit || par.t1 <= 0 || par.t2 <= par.t1
    error('Q2:InvalidStepStrategy', ...
        '要求 0<=j1<=j2<=j3<=j_max 且 0<t1<t2。');
end
if t < par.t1
    j = par.j1;
elseif t < par.t2
    j = par.j2;
else
    j = par.j3;
end
end
