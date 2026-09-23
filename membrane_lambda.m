%--------------------------------------------------------------%
% @function: 由膜内水质量浓度计算膜含水量Lambda，即题目中式(21)
% @author:   PJ, GPT
% @date:     20260923
% @input:    mw->PEM中水质量浓度
%            p->参数列表
% @output:   Lambda->膜含水量
%--------------------------------------------------------------%

function lambda = membrane_lambda(mw, p)

    % 膜内水质量浓度不能为负
    if any(mw(:) < 0)
        error('Membrane water concentration cannot be negative.');
    end
    
    % EW需使用SI单位kg/mol，附件1中1000g/mol即1kg/mol
    lambda = p.water.EW .* mw ./ (p.water.rho_pem .* p.const.Mw);

end
