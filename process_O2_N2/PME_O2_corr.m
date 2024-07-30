% PME oxygen calibration

% O2 is reported in mg/L at freshwater salinity, need to convert to umol/kg
% at the true salinity

% ---------
% method 1: mg/L freshwater --> % sat --> umol/kg at true salinity

O2_mgL_PME = 8.3;
temp_PME = 25;
sal_PME = 0.*O2_mgL_PME; % set to 0 salinity for all data points
press_PME = 0.*O2_mgL_PME; % set to 0 pressure for all data points

sal_meas = 30; %salinity measured with external sensor

MW_O2 = 31.998; % g/mol
g_per_mg = 1/1000;
umol_per_mol = 1e6;
m3_per_kg = 1/1000;
mgL_to_umolL = 1/MW_O2 .* g_per_mg .* umol_per_mol;
PME_O2_umolL = O2_mgL_PME .* mgL_to_umolL;
PME_O2_umolkg = PME_O2_umolL .* sw_dens(temp_PME,sal_PME,press_PME) .* m3_per_kg;
PME_O2_sat_ratio = PME_O2_umolkg ./ O2sol(sal_PME,temp_PME);

% salinity corrected O2 in umolkg
% note this is based on in situ temperature not conservative or potential
% temperature
O2_umolkg_sal_corr = PME_O2_sat_ratio .* O2sol(sal_meas,temp_PME);

% -------------
% method 2: 

