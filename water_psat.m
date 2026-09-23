function psat = water_psat(T)

    %   Saturated water-vapor pressure using Buck correlation
    %
    %   psat = water_psat(T)
    %
    %   INPUT
    %   -----
    %   T : temperature [K], scalar/vector/matrix
    %
    %   OUTPUT
    %   ------
    %   psat : saturation vapor pressure [Pa]

    if any(T(:) <= 0)
        error('water_psat:InvalidTemperature', ...
            'Temperature must be given in Kelvin and greater than zero.');
    end
    
    Tc = T - 273.15;
    
    psat = zeros(size(T));
    
    idxWarm = Tc >= 0;
    idxCold = ~idxWarm;
    
    %% Liquid-water branch
    Tw = Tc(idxWarm);
    
    psat(idxWarm) = 611.21 .* exp( ...
        (18.678 - Tw./234.5) .* ...
        Tw ./ (257.14 + Tw) );
    
    %% Ice branch
    Tcold = Tc(idxCold);
    
    psat(idxCold) = 611.15 .* exp( ...
        (23.036 - Tcold./333.7) .* ...
        Tcold ./ (279.82 + Tcold) );

end