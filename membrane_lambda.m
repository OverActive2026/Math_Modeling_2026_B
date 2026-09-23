function lambda = membrane_lambda(mw, p)

    %   Calculate membrane water content lambda
    %
    %   lambda = EW * mw / (rho_pem * Mw)
    %
    %   mw : water concentration in PEM [kg/m^3]
    %
    %   NOTE:
    %   p.water.EW must use SI units [kg/mol].
    %
    %   Attachment 1:
    %       EW = 1000 g/mol = 1 kg/mol
    
    if any(mw(:) < 0)
        error('Membrane water concentration cannot be negative.');
    end
    
    lambda = ...
        p.water.EW .* mw ./ ...
        (p.water.rho_pem .* p.const.Mw);

end