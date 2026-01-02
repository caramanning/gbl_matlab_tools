load LIS_gas_flux_Aug.mat;
load LIS_gas_flux_Oct.mat;
load LIS_gas_flux_May.mat;

matches = contains(string(LIS_gas_flux_Aug.station), "EXRX");  % logical vector of matching categories
CH4_MID4_ms = [mean(LIS_gas_flux_Aug.ch4_nmolkg(matches)) std(LIS_gas_flux_Aug.ch4_nmolkg(matches)) std(LIS_gas_flux_Aug.ch4_nmolkg(matches))./mean(LIS_gas_flux_Aug.ch4_nmolkg(matches)).*100]
FCH4_MID4_ms = [mean(LIS_gas_flux_Aug.F_CH4_15(matches)) std(LIS_gas_flux_Aug.F_CH4_15(matches)) std(LIS_gas_flux_Aug.F_CH4_15(matches))./mean(LIS_gas_flux_Aug.F_CH4_15(matches)).*100]

N2O_MID4_ms = [mean(LIS_gas_flux_Aug.n2o_nmolkg(matches)) std(LIS_gas_flux_Aug.n2o_nmolkg(matches)) std(LIS_gas_flux_Aug.n2o_nmolkg(matches))./mean(LIS_gas_flux_Aug.n2o_nmolkg(matches)).*100]
FN2O_MID4_ms = [mean(LIS_gas_flux_Aug.F_N2O_15(matches)) std(LIS_gas_flux_Aug.F_N2O_15(matches)) std(LIS_gas_flux_Aug.F_N2O_15(matches))./mean(LIS_gas_flux_Aug.F_N2O_15(matches)).*100]

matches = contains(string(LIS_gas_flux_Oct.station), "EXRX");  % logical vector of matching categories
CH4_MID4_ms = [mean(LIS_gas_flux_Oct.ch4_nmolkg(matches)) std(LIS_gas_flux_Oct.ch4_nmolkg(matches)) std(LIS_gas_flux_Oct.ch4_nmolkg(matches))./mean(LIS_gas_flux_Oct.ch4_nmolkg(matches)).*100]
FCH4_MID4_ms = [mean(LIS_gas_flux_Oct.F_CH4_15(matches)) std(LIS_gas_flux_Oct.F_CH4_15(matches)) std(LIS_gas_flux_Oct.F_CH4_15(matches))./mean(LIS_gas_flux_Oct.F_CH4_15(matches)).*100]

N2O_MID4_ms = [mean(LIS_gas_flux_Oct.n2o_nmolkg(matches)) std(LIS_gas_flux_Oct.n2o_nmolkg(matches)) std(LIS_gas_flux_Oct.n2o_nmolkg(matches))./mean(LIS_gas_flux_Oct.n2o_nmolkg(matches)).*100]
FN2O_MID4_ms = [mean(LIS_gas_flux_Oct.F_N2O_15(matches)) std(LIS_gas_flux_Oct.F_N2O_15(matches)) std(LIS_gas_flux_Oct.F_N2O_15(matches))./mean(LIS_gas_flux_Oct.F_N2O_15(matches)).*100]

matches = contains(string(LIS_gas_flux_May.station), "EXRX");  % logical vector of matching categories
CH4_MID4_ms = [mean(LIS_gas_flux_May.ch4_nmolkg(matches)) std(LIS_gas_flux_May.ch4_nmolkg(matches)) std(LIS_gas_flux_May.ch4_nmolkg(matches))./mean(LIS_gas_flux_May.ch4_nmolkg(matches)).*100]
FCH4_MID4_ms = [mean(LIS_gas_flux_May.F_CH4_15(matches)) std(LIS_gas_flux_May.F_CH4_15(matches)) std(LIS_gas_flux_May.F_CH4_15(matches))./mean(LIS_gas_flux_May.F_CH4_15(matches)).*100]

N2O_MID4_ms = [mean(LIS_gas_flux_May.n2o_nmolkg(matches)) std(LIS_gas_flux_May.n2o_nmolkg(matches)) std(LIS_gas_flux_May.n2o_nmolkg(matches))./mean(LIS_gas_flux_May.n2o_nmolkg(matches)).*100]
FN2O_MID4_ms = [mean(LIS_gas_flux_May.F_N2O_15(matches)) std(LIS_gas_flux_May.F_N2O_15(matches)) std(LIS_gas_flux_May.F_N2O_15(matches))./mean(LIS_gas_flux_May.F_N2O_15(matches)).*100]

%%
matches = contains(string(LIS_gas_flux_Aug.station), "MID4");  % logical vector of matching categories
FCH4_MID4_ms = [mean(LIS_gas_flux_Aug.F_CH4_inst(matches)) std(LIS_gas_flux_Aug.F_CH4_inst(matches)) std(LIS_gas_flux_Aug.F_CH4_inst(matches))./mean(LIS_gas_flux_Aug.F_CH4_inst(matches)).*100]

FN2O_MID4_ms = [mean(LIS_gas_flux_Aug.F_N2O_inst(matches)) std(LIS_gas_flux_Aug.F_N2O_inst(matches)) std(LIS_gas_flux_Aug.F_N2O_inst(matches))./mean(LIS_gas_flux_Aug.F_N2O_inst(matches)).*100]

matches = contains(string(LIS_gas_flux_Oct.station), "MID4");  % logical vector of matching categories
FCH4_MID4_ms = [mean(LIS_gas_flux_Oct.F_CH4_inst(matches)) std(LIS_gas_flux_Oct.F_CH4_inst(matches)) std(LIS_gas_flux_Oct.F_CH4_inst(matches))./mean(LIS_gas_flux_Oct.F_CH4_inst(matches)).*100]

FN2O_MID4_ms = [mean(LIS_gas_flux_Oct.F_N2O_inst(matches)) std(LIS_gas_flux_Oct.F_N2O_inst(matches)) std(LIS_gas_flux_Oct.F_N2O_inst(matches))./mean(LIS_gas_flux_Oct.F_N2O_inst(matches)).*100]

matches = contains(string(LIS_gas_flux_May.station), "MID4");  % logical vector of matching categories
FCH4_MID4_ms = [mean(LIS_gas_flux_May.F_CH4_inst(matches)) std(LIS_gas_flux_May.F_CH4_inst(matches)) std(LIS_gas_flux_May.F_CH4_inst(matches))./mean(LIS_gas_flux_May.F_CH4_inst(matches)).*100]

FN2O_MID4_ms = [mean(LIS_gas_flux_May.F_N2O_inst(matches)) std(LIS_gas_flux_May.F_N2O_inst(matches)) std(LIS_gas_flux_May.F_N2O_inst(matches))./mean(LIS_gas_flux_May.F_N2O_inst(matches)).*100]

%%
matches = contains(string(LIS_gas_flux_Aug.station), "MID4");  % logical vector of matching categories
k_wt_15_CH4_ms = [mean(LIS_gas_flux_Aug.k_wt_15_CH4(matches)) std(LIS_gas_flux_Aug.k_wt_15_CH4(matches)) std(LIS_gas_flux_Aug.k_wt_15_CH4(matches))./mean(LIS_gas_flux_Aug.k_wt_15_CH4(matches)).*100]

FN2O_MID4_ms = [mean(LIS_gas_flux_Aug.k_wt_15_N2O(matches)) std(LIS_gas_flux_Aug.k_wt_15_N2O(matches)) std(LIS_gas_flux_Aug.k_wt_15_N2O(matches))./mean(LIS_gas_flux_Aug.k_wt_15_N2O(matches)).*100]

matches = contains(string(LIS_gas_flux_Oct.station), "MID4");  % logical vector of matching categories
FCH4_MID4_ms = [mean(LIS_gas_flux_Oct.k_wt_15_CH4(matches)) std(LIS_gas_flux_Oct.k_wt_15_CH4(matches)) std(LIS_gas_flux_Oct.k_wt_15_CH4(matches))./mean(LIS_gas_flux_Oct.k_wt_15_CH4(matches)).*100]

FN2O_MID4_ms = [mean(LIS_gas_flux_Oct.k_wt_15_N2O(matches)) std(LIS_gas_flux_Oct.k_wt_15_N2O(matches)) std(LIS_gas_flux_Oct.k_wt_15_N2O(matches))./mean(LIS_gas_flux_Oct.k_wt_15_N2O(matches)).*100]

matches = contains(string(LIS_gas_flux_May.station), "MID4");  % logical vector of matching categories
FCH4_MID4_ms = [mean(LIS_gas_flux_May.k_wt_15_CH4(matches)) std(LIS_gas_flux_May.k_wt_15_CH4(matches)) std(LIS_gas_flux_May.k_wt_15_CH4(matches))./mean(LIS_gas_flux_May.k_wt_15_CH4(matches)).*100]

FN2O_MID4_ms = [mean(LIS_gas_flux_May.k_wt_15_N2O(matches)) std(LIS_gas_flux_May.k_wt_15_N2O(matches)) std(LIS_gas_flux_May.k_wt_15_N2O(matches))./mean(LIS_gas_flux_May.k_wt_15_N2O(matches)).*100]

