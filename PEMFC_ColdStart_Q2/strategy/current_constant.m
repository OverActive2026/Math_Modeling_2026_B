function j = current_constant(t,par)
%CURRENT_CONSTANT 恒流加载；t [s]，j [A/cm^2]。
p = parameters_Q2();
if isfield(par,'j_max'), limit = par.j_max; else, limit = p.j_max; end
if ~isscalar(t) || ~isfinite(t) || t < 0 || ...
        ~isfield(par,'J') || ~isscalar(par.J) || ...
        ~isfinite(par.J) || par.J < 0 || par.J > limit
    error('Q2:InvalidConstantStrategy','要求 t>=0 且 0<=J<=j_max。');
end
j = par.J;
end
