%%
load LIS_gas_flux_Aug_avg.mat;
%%
% CH4 flux to CO2 equivalents
mol_per_umol = 1e-6;

% At EXCR
k_N2O =  [1.28 1.56 1.20];
k_CH4 = [1.32 1.60 1.23];

CH4_flux_umol = [236 236 80];
N2O_flux_umol = [2.96 4.85 5.81];

% At MID4
% Aug Oct May]
k_N2O =  [1.32 1.45 1.24];
k_CH4 = [1.39 1.49 1.28];

CH4_flux_umol = [69 72 59];
N2O_flux_umol = [1.89 1.82 4.90];


% CH4 flux 
CH4_g_per_mol = 16;
CH4_GWP_100 = 80; % over 100 years
CH4_flux_ug = CH4_flux_umol .* CH4_g_per_mol;
CH4_CO2_eq_100 = CH4_flux_ug .* CH4_GWP_100 .* mol_per_umol % g/m2/d

mean(CH4_CO2_eq_100)

% N2O flux 
N2O_g_per_mol = 44;
N2O_GWP_100 = 273; % over 100 years
N2O_flux_ug = N2O_flux_umol .* N2O_g_per_mol;
N2O_CO2_eq_100 = N2O_flux_ug .* N2O_GWP_100 .* mol_per_umol % g/m2/d

mean(N2O_CO2_eq_100)
%%
% CO2 flux
% 1500-400
S = 26;
T = 22;
[~,Sc] = gasmoldiff(S,T,'CO2');
CO2_g_per_mol = 44;
pCO2_atm = 400e-6;
pCO2_surf = 1300e-6;
sol_CO2 = co2_solubility(T,S,'vol').*1000; % inputs T,S,'vol', output mol/L/atm
CO2_flux = 1.4 .* sol_CO2.* (pCO2_surf - pCO2_atm); % m/d * mol/m3/atm .* atm = mol/m2/d
CO2_flux_g = CO2_flux .* CO2_g_per_mol

% flux at MID4

% 562 at k=1.4 = 0.3
% 672 at k=1.4 = 0.5
% 586 at k=1.5 = 0.4
% 1333 at k=1.4 = 1.84

% flux at EXCR
% 800 at k=1.4 
% 1300 at k=1.4
% 1200 at k=1.5 
% 1300 at k=1.4

