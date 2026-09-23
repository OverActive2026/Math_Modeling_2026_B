function Deff = calc_gas_diffusivity(Dref, T, pGas, eps_g, p)

    %   Effective gas diffusivity in porous media
    %
    %   Deff = calc_gas_diffusivity(Dref,T,pGas,eps_g,p)
    %
    %   Formula:
    %
    %   Deff = Dref * (T/Tref)^1.75 ...
    %                * (p0/pGas) ...
    %                * eps_g^1.5
    %
    %   INPUT
    %   -----
    %   Dref  : reference diffusivity [m^2/s]
    %   T     : local temperature [K]
    %   pGas  : local absolute pressure [Pa]
    %   eps_g : available gas porosity [-]
    
    if any(T(:) <= 0)
        error('Temperature must be greater than zero.');
    end
    
    if any(pGas(:) <= 0)
        error('Gas pressure must be greater than zero.');
    end
    
    if any(eps_g(:) < 0)
        error('Negative gas porosity detected.');
    end
    
    Deff = Dref .* ...
        (T ./ p.const.Tref).^p.trans.Texp .* ...
        (p.const.p0 ./ pGas) .* ...
        eps_g.^p.trans.epsExp;

end