load exrx2023_met.mat

load LISAug23_CH4N2O_CTD;
LIS = LISAug23_CH4N2O_CTD;

load LISOct23_CH4N2O_CTD;
LIS = LISOct23_CH4N2O_CTD;

h_WSPD = 3.5; % wind speed height in meters
ms_per_kts = 0.514444; % conversion factor to go from kts to m/s
exrx.u10_ms = exrx.windSpd_Kts .* ms_per_kts .*(10/h_WSPD)^0.11;


% cruise dates
%dc = datetime(2023,8,2); % Aug cruise date
dc = datetime(2023,10,19); % Oct cruise date
ds = dc-15;
de = dc + 15;

%dsdt = datetime(ds);
%dsde = datetime(de);

A = find(exrx.TIMESTAMP>=ds & exrx.TIMESTAMP<=de); % 15 days before and after cruise

u10_sq = exrx.u10_ms(A).^2;
avg_u10_sq = mean(u10_sq);
u10_monthly = sqrt(avg_u10_sq);



LIS.n2o_nmolkg = LIS.mean_N2O_nM./sw_dens(LIS.S,LIS.T,LIS.P).*1000;
LIS.ch4_nmolkg = LIS.mean_CH4_nM./sw_dens(LIS.S,LIS.T,LIS.P).*1000;


% need to use the exact atmospheric concentration corrected for water vapor here, just approximating
% for now
N2Oatm = 330e-9;
LIS.n2o_eq_nmolkg = N2Osol(LIS.S,LIS.T,N2Oatm).*1000;

CH4atm = 1920e-9;
LIS.ch4_eq_nmolkg = CH4sol(LIS.S,LIS.T,CH4atm)'.*1000;

LIS.Dch4 = (LIS.ch4_nmolkg - LIS.ch4_eq_nmolkg)./LIS.ch4_eq_nmolkg.*100;
LIS.Dn2o = (LIS.n2o_nmolkg - LIS.n2o_eq_nmolkg)./LIS.n2o_eq_nmolkg.*100;


si = find(LIS.Depth<=3); % surface indices
u10_monthly_rep = repmat(u10_monthly,numel(si),1);


[~,LIS.Sc_CH4] = gasmoldiff(LIS.S,LIS.T,'CH4');
[~,LIS.Sc_N2O] = gasmoldiff(LIS.S,LIS.T,'N2O');

LIS.k_CH4 = nan.*LIS.S;
LIS.k_N2O = nan.*LIS.S;
spd = 60*60*24;
LIS.k_CH4_md(si) = kgas(u10_monthly_rep,LIS.Sc_CH4(si),'W14').*spd; %k in m/d
LIS.k_N2O_md(si) = kgas(u10_monthly_rep,LIS.Sc_N2O(si),'W14').*spd; %k in m/d

LIS.F_CH4_umolm2d = nan.*LIS.S;
LIS.F_N2O_umolm2d = nan.*LIS.S;
% Flux = k .* (CH4meas - CH4eq)
% Flux in umol m-2 d-1 = m d-1 .* nmol/kg .* umol/nmol .* kg/m3 

umol_per_nmol = 1e-3;
LIS.F_CH4_umolm2d(si) = LIS.k_CH4_md(si) .* (LIS.ch4_nmolkg(si) - LIS.ch4_eq_nmolkg(si)) .* umol_per_nmol .* LIS.Dens(si);
LIS.F_N2O_umolm2d(si) = LIS.k_N2O_md(si) .* (LIS.n2o_nmolkg(si) - LIS.n2o_eq_nmolkg(si)) .* umol_per_nmol .* LIS.Dens(si);



