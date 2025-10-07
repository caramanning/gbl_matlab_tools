% load LIS_gas_flux_Aug.mat;
% LISo = LIS_gas_flux_Aug;
% cr = 'RVCT_20230802_';
% fn = 'PANGAEA_A23_flux.csv';

% load LIS_gas_flux_Oct.mat
% cr = 'RVCT_20231019_';
% LISo = LIS_gas_flux_Oct;
% fn = 'PANGAEA_O23_flux.csv';


load LIS_gas_flux_May.mat
cr = 'RVCT_20240522_';
LISo = LIS_gas_flux_May;
fn = 'PANGAEA_M24_flux.csv';


Event = repmat(string(cr),length(LISo.station),1);

LISn = table(Event);
LISn.Station = LISo.station;
LISn.dt_ISO = datestr(LISo.dt, 'yyyy-mm-ddTHH:MM:SS');
LISn.Latitude = LISo.lat;
LISn.Longitude = LISo.lon;
LISn.Depth = LISo.depth;
LISn.Pressure = LISo.P;
LISn.Salinity = LISo.S;
LISn.Temperature = LISo.T;
LISn.MLD = LISo.mld;
LISn.CH4_diss = LISo.ch4_nmolkg;
LISn.CH4_std_dev = LISo.ch4_std_nmolkg;
LISn.CH4_eq = LISo.ch4_eq_nmolkg;
LISn.N2O_diss = LISo.n2o_nmolkg;
LISn.N2O_std_dev = LISo.n2o_std_nmolkg;
LISn.N2O_eq = LISo.n2o_eq_nmolkg;
LISn.slp_15 = LISo.slp_15;
LISn.k_wt_15_CH4 = LISo.k_wt_15_CH4;
LISn.F_CH4_15 = LISo.F_CH4_15;
LISn.k_wt_15_N2O = LISo.k_wt_15_N2O;
LISn.F_N2O_15 = LISo.F_N2O_15;

for i=1:length(LISn.Event)
    LISn.Event(i) = strcat(cr,string(LISn.Station(i)));
end;

writetable(LISn,fn);

%%

LISn.O2 = LISo.O2_umolkg;
LISn.StationDepth = LISo.StationDepth;

LISn.Cast = LISo.CastNum;
LISn.Bottle = LISo.NiskinNum; 




