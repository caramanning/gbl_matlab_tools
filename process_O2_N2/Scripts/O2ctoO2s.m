function O2sat=O2ctoO2s(O2conc,T,S,P,p_atm)
%function O2sat=O2ctoO2s(O2conc,T,S,P,p_atm)
%
% convert molar oxygen concentration to oxygen saturation
%
% inputs:
%   O2conc - oxygen concentration in umol L-1
%   T      - temperature in �C
%   S      - salinity (PSS-78)
%   P      - hydrostatic pressure in dbar (default: 0 dbar)
%   p_atm  - atmospheric (air) pressure in mbar (default: 1013.25 mbar)
%
% output:
%   O2sat  - oxygen saturation in %
%
% according to recommendations by SCOR WG 142 "Quality Control Procedures
% for Oxygen and Other Biogeochemical Sensors on Floats and Gliders"
%
% Henry Bittig
% Laboratoire d'Oc�anographie de Villefranche-sur-Mer, France
% bittig@obs-vlfr.fr
% 28.10.2015
% 19.04.2018, v1.1, fixed typo in B2 exponent

% set a few input defaults
if nargin<5, p_atm = 1013.25; end
if nargin<4, P     = 0;       end

pH2Osat = 1013.25.*(exp(24.4543-(67.4509*(100./(T+273.15)))-(4.8489*log(((273.15+T)./100)))-0.000544.*S)); % saturated water vapor in mbar
sca_T   = log((298.15-T)./(273.15+T)); % scaled temperature for use in TCorr and SCorr
TCorr   = 44.6596.*exp(2.00907+3.22014.*sca_T+4.05010.*sca_T.^2+4.94457.*sca_T.^3-2.56847e-1.*sca_T.^4+3.88767.*sca_T.^5); % temperature correction part from Garcia and Gordon (1992), Benson and Krause (1984) refit mL(STP) L-1; and conversion from mL(STP) L-1 to umol L-1
Scorr   = exp(S.*(-6.24523e-3-7.37614e-3.*sca_T-1.03410e-2.*sca_T.^2-8.17083e-3.*sca_T.^3)-4.88682e-7.*S.^2); % salinity correction part from Garcia and Gordon (1992), Benson and Krause (1984) refit ml(STP) L-1
Vm      = 0.317; % molar volume of O2 in m3 mol-1 Pa dbar-1 (Enns et al. 1965)
R       = 8.314; % universal gas constant in J mol-1 K-1

O2sat=O2conc.*100./(TCorr.*Scorr)./(p_atm-pH2Osat).*(1013.25-pH2Osat).*exp(Vm.*P./(R.*(T+273.15)));
