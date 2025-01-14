% correction for CH4 and N2O to match with updated standard concentrations
load('LISAug23_CH4N2O.mat')
LISAug23_CH4N2O_uncal = LISAug23_CH4N2O;
save LISAug23_CH4N2O_uncal.mat LISAug23_CH4N2O_uncal;

load('LISAug23_CH4N2O_CTD.mat')
LISAug23_CH4N2O_CTD_uncal = LISAug23_CH4N2O_CTD;
save LISAug23_CH4N2O_CTD_uncal.mat LISAug23_CH4N2O_CTD_uncal;


%%
CH4_cf = 0.9687;
N2O_cf = 0.9898;
LISAug23_CH4N2O.mean_CH4_nM = LISAug23_CH4N2O.mean_CH4_nM .* CH4_cf;
LISAug23_CH4N2O.std_CH4_nM = LISAug23_CH4N2O.std_CH4_nM .* CH4_cf;
LISAug23_CH4N2O.mean_N2O_nM = LISAug23_CH4N2O.mean_N2O_nM .* N2O_cf;
LISAug23_CH4N2O.std_N2O_nM = LISAug23_CH4N2O.std_N2O_nM .* N2O_cf;

save LISAug23_CH4N2O.mat LISAug23_CH4N2O;


%%
LISAug23_CH4N2O_CTD.mean_CH4_nM = LISAug23_CH4N2O_CTD.mean_CH4_nM .* CH4_cf;
LISAug23_CH4N2O_CTD.std_CH4_nM = LISAug23_CH4N2O_CTD.std_CH4_nM .* CH4_cf;
LISAug23_CH4N2O_CTD.mean_N2O_nM = LISAug23_CH4N2O_CTD.mean_N2O_nM .* N2O_cf;
LISAug23_CH4N2O_CTD.std_N2O_nM = LISAug23_CH4N2O_CTD.std_N2O_nM .* N2O_cf;

save LISAug23_CH4N2O_CTD.mat LISAug23_CH4N2O_CTD;

%%
%now correct the O2 data based on winkler titrations.
% for August 2024
% Difficult to calibrate due to the very strong vertical gradients
% therefore only using near-surface measurements as the niskins and sensors
% are not perfectly aligned
% From LIS August 2023 oxygen results.xlsx
% m = 1.0387
% b = -7.7438
% O2cal = O2uncal.* m + b

% For October 2023, again focusing on surface samples due to uncertainties
% in subsurface
% From LIS Oct 2023 titrations.xlsx
% m = 1.0283; 
% b = -3.1327;
%O2cal = O2uncal.*m + b

% For May 2024, results are from lis Cruise May 2024 Titration.xlsx
% m = 1.0438
% b = -6.0476
% O2cal = O2uncal.*m + b


%%
% correction for O2 to match with updated standard concentrations
load('LISAug23_CH4N2O_CTD.mat')

mAug = 1.0387;
bAug = -7.7438;
LISAug23_CH4N2O_CTD.O2_umolkg = LISAug23_CH4N2O_CTD.O2_umolkg .* mAug + bAug;

save LISAug23_CH4N2O_CTD.mat LISAug23_CH4N2O_CTD;

%%
load('LISOct23_CH4N2O_CTD.mat')

mOct = 1.0283;
bOct = -3.1327;
LISOct23_CH4N2O_CTD.O2_umolkg = LISOct23_CH4N2O_CTD.O2_umolkg .* mOct + bOct;

save LISOct23_CH4N2O_CTD.mat LISOct23_CH4N2O_CTD;

%%
load('LISMay24_CH4N2O_CTD.mat')

mMay = 1.0438;
bMay = -6.0476;
LISMay24_CH4N2O_CTD.O2_umolkg = LISMay24_CH4N2O_CTD.O2_umolkg .* mMay + bMay;

save LISMay24_CH4N2O_CTD.mat LISMay24_CH4N2O_CTD;