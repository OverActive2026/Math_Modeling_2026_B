function better = q3_compare_candidates(a,b)
% Feasibility-first strict ordering. A tie keeps the incumbent b.
    rank = @(s) 2*strcmp(s,'success') + strcmp(s,'infeasible');
    ra = rank(char(a.status));
    rb = rank(char(b.status));
    if ra ~= rb
        better = ra > rb;
        return;
    end
    if ra == 2
        scale = max([1,abs(a.energyJ),abs(b.energyJ)]);
        if abs(a.energyJ-b.energyJ) > 1e-6*scale
            better = a.energyJ < b.energyJ;
        elseif abs(a.safetyMargin-b.safetyMargin) > 1e-9
            better = a.safetyMargin > b.safetyMargin;
        else
            better = a.successTimeS < b.successTimeS-1e-7;
        end
    elseif ra == 1
        if abs(a.violation-b.violation) > 1e-9
            better = a.violation < b.violation;
        else
            better = a.safetyMargin > b.safetyMargin+1e-9;
        end
    else
        better = false;
    end
end
