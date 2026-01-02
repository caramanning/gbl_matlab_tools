function [Ko] = co2_solubility(T,S,varargin),
%- Calculate solubility of CO2 as a non-ideal gas in sea water.  This
%- normally requires correction using the virial correction coefficient (B),
%- but we can avoid doing that if we're close to 1 atm.
%- OUTPUT:  Ko, where C_liquid = Ko*Pp,  Pp = gas partial pressure
%-
%- INPUTS:  temp (degrees kelvin).  Salinity (psu)
%
%- SOURCE:  Weiss, R.F, (1974).  Carbon Dioxide in water and seawater: The
%           Solubility of a non-ideal gas.  Marine Chemistry V.2, pp.203-15
%
%-          We are required to use the fugacity eqn.
%-
%- OUTPUT:  mols/L/atm

T = T + 273.15;

if nargin>2,
    mass_or_vol = varargin{1};
else
    mass_or_vol = 'vol';
end

if mass_or_vol == 'vol', %- this gives mols/L/atm

    a1 = -58.0931;
    a2 = 90.5069;
    a3 = 22.2940;

    b1 = 0.027766;
    b2 = -0.025888;
    b3 = 0.0050578;

elseif mass_or_vol == 'mas',  %- this gives mols/kg/atm
    a1 = -60.2409;
    a2 = 93.4517;
    a3 = 23.3585;
    
    b1 = 0.023517;
    b2 = -0.023656;
    b3 = 0.0047036;
end

lnKo = a1 + a2*100./T + a3*log(T/100) + S.*[b1 + b2*T/100 + b3*(T/100).^2];
Ko = exp(lnKo);
    

%- Gas volume and fugacity pressure corrections based on virial eqn of
%- state
%B = -1636.75 + 12.0408*T - 3.27957e-2*T.^2 + 3.16528e-5*T.^3;
%R = 8.3141;
%fco2 = P.*exp(B.*P./R./T);