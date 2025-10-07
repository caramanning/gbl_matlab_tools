% load LISAug23_CH4N2O_CTD.mat;
% cr = 'RVCT_20230802_';
% LISo = LISAug23_CH4N2O_CTD;

% load LISOct23_CH4N2O_CTD.mat;
% cr = 'RVCT_20231019_';
% LISo = LISOct23_CH4N2O_CTD;

load LISMay24_CH4N2O_CTD.mat;
cr = 'RVCT_20240522_';
LISo = LISMay24_CH4N2O_CTD;

Event = LISo.Cruise;
LISn = table(Event);
LISn.Station = LISo.Station;
LISn.dt_ISO = datestr(LISo.datetime, 'yyyy-mm-ddTHH:MM:SS');
LISn.Latitude = LISo.Lat;
LISn.Longitude = LISo.Lon;
LISn.Depth = LISo.Depth;
LISn.Pressure = LISo.P;
LISn.Salinity = LISo.S;
LISn.Temperature = LISo.T;
LISn.CH4_diss = LISo.CH4_mean_nmolkg;
LISn.CH4_std_dev = LISo.CH4_std_nmolkg;
LISn.CH4_eq = LISo.CH4_eq_nmolkg;
LISn.N2O_diss = LISo.N2O_mean_nmolkg;
LISn.N2O_std_dev = LISo.N2O_std_nmolkg;
LISn.N2O_eq = LISo.N2O_eq_nmolkg;
LISn.O2 = LISo.O2_umolkg;
LISn.StationDepth = LISo.StationDepth;
LISn.MLD = LISo.mld;
LISn.Cast = LISo.CastNum;
LISn.Bottle = LISo.NiskinNum; 


for i=1:length(LISn.Event)
    LISn.Event(i) = strcat(cr,string(LISn.Station(i)));
end;

writetable(LISn,'PANGAEA_M24.csv');