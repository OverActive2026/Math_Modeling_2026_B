function newState = stack_step(state,j,p)
%STACK_STEP 同一电流密度下推进五片状态；各片状态独立保存。
if ~isscalar(j) || ~isfinite(j) || j < 0 || j > p.j_max
    error('Q2:InvalidCurrent','电堆电流密度超出允许范围 [A/cm^2]。');
end
Qstack = stack_heat_transfer(state.T,p); % 五片净热功率 [W]
fields = {'T','V','ice_max','lambda_mem','m_v','m_l','m_i'};
for f = 1:numel(fields)
    values = state.(fields{f});
    if ~isequal(size(values),[p.Ncell,1])
        error('Q2:InvalidState','%s 必须为五片各自的列向量。',fields{f});
    end
end
newState = state;
for k = 1:p.Ncell
    oldCell = struct();
    for f = 1:numel(fields)
        values = state.(fields{f});
        oldCell.(fields{f}) = values(k);
    end
    newCell = single_cell_step(oldCell,j,Qstack(k),p);
    for f = 1:numel(fields)
        values = newState.(fields{f});
        values(k) = newCell.(fields{f});
        newState.(fields{f}) = values;
    end
end
end
