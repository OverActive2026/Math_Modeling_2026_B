%------------------------------------------------------------------%
% @function: 按题目给定的等效模型计算电池的等效体积热容和面外导热系数
% @author:   PJ, GPT
% @date:     20260923
% @input:    g->空间网格与五层结构信息
%            p->参数列表
% @output:   rhoCp->等效体积热容网格场
%            kEff->等效面外导热系数网格场
%            out->分层热物性与等效参数诊断信息
%------------------------------------------------------------------%

function [rhoCp, kEff, out] = calc_thermal_properties(g, p)

    % 读取五层材料的密度、比热容、导热系数和厚度
    rho = p.mat.rho(:);
    cp  = p.mat.cp(:);
    k   = p.mat.k(:);

    delta = g.layerThickness(:);
    nLayer = numel(delta);

    % 每组材料参数都必须与层数一致
    if numel(rho) ~= nLayer || numel(cp) ~= nLayer || numel(k) ~= nLayer
        error('calc_thermal_properties:SizeMismatch', 'p.mat.rho, p.mat.cp and p.mat.k must each contain %d layer values.', nLayer);
    end

    if any(~isfinite(rho)) || any(~isfinite(cp)) || any(~isfinite(k))
        error('calc_thermal_properties:MissingParameter', 'Thermal material properties contain NaN/Inf. Populate p.mat.rho, p.mat.cp and p.mat.k first.');
    end

    if any(rho <= 0)
        error('All layer densities must be positive.');
    end

    if any(cp <= 0)
        error('All layer specific heats must be positive.');
    end

    if any(k <= 0)
        error('All layer thermal conductivities must be positive.');
    end

    if any(delta <= 0)
        error('All layer thicknesses must be positive.');
    end

    % 按各层厚度加权平均计算等效体积热容
    layerRhoCp = rho .* cp;
    Ltotal = sum(delta);
    rhoCpScalar = sum(layerRhoCp .* delta) ./ Ltotal;

    % 五层材料在面外方向上串联导热，因此先求单位面积总热阻
    RthPerArea = sum(delta ./ k);
    kEffScalar = Ltotal ./ RthPerArea;

    % MVP将等效热物性扩展到所有网格，当前不考虑液态水和冰对热物性的影响
    rhoCp = rhoCpScalar .* ones(g.N,1);
    kEff = kEffScalar .* ones(g.N,1);

    % 输出分层与等效热物性，用于调试和结果检查
    out.rhoCpScalar = rhoCpScalar;
    out.kEffScalar  = kEffScalar;

    out.layerRhoCp = layerRhoCp;
    out.layerThermalResistance = delta ./ k;
    out.thermalResistancePerArea = RthPerArea;
    out.Ltotal = Ltotal;

end
