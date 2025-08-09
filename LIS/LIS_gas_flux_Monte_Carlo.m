load LIS_gas_flux_May_avg.mat;
gf = gfM_avg;
k_rsd = 0.2;
wind_rsd = 0.03;
n2o_rsd = 0.04;
ch4_rsd = 0.12;

dens = 1017;

sd_F_N2O = nan(length(gf.station),1);
sd_F_CH4 = nan(length(gf.station),1);

% diffusivity, 3% for CH4, 5% for N2O?
% solubility, 

n=1000;
for i = 1:length(gf.station)
    mc_k_wt_15_N2O = gf.k_wt_15_N2O(i) + k_rsd.*gf.k_wt_15_N2O(i).*randn(n,1); % uncertainty in kn2o linear factor
    mc_k_wt_15_CH4 = gf.k_wt_15_CH4(i) + k_rsd .* gf.k_wt_15_CH4(i) .* randn(n, 1); % uncertainty in kch4 linear factor
    mc_wind = (1 + wind_rsd .* randn(n, 1)); % additional factor to add based on wind speed squared
    mc_n2o = gf.n2o_nmolkg(i) + n2o_rsd .* gf.n2o_nmolkg(i) .* randn(n, 1);
    mc_ch4 = gf.ch4_nmolkg(i) + ch4_rsd .* gf.ch4_nmolkg(i) .* randn(n, 1);

    mc_F_N2O = mc_k_wt_15_N2O .* mc_wind .* (mc_n2o - repmat(gf.n2o_eq_nmolkg(i),n,1)) .* dens ./1000;
    mc_F_CH4 = mc_k_wt_15_CH4 .* mc_wind .* (mc_ch4 - repmat(gf.ch4_eq_nmolkg(i),n,1)) .* dens ./1000;

    sd_F_N2O(i) = std(mc_F_N2O);
    sd_F_CH4(i) = std(mc_F_CH4);

end;

rsd_F_N2O = sd_F_N2O./gf.F_N2O_15 .*100;
rsd_F_CH4 = sd_F_CH4./gf.F_CH4_15 .*100;


[mean(rsd_F_N2O) mean(rsd_F_CH4)]