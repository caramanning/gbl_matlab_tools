% this code calculates the gas flux based on NDBC mooring data collected at
% one location (mooring 44022, EXRX)
% data from mooring 44040, WLIS, was not available for the desired time
% period

clear all; clc;

load LISAug23_CH4N2O_CTD.mat;

surf_i = find(LISAug23_CH4N2O_CTD.Depth<3);

LIS_s = LISAug23_CH4N2O_CTD(surf_i,:);

dl = datetime(2023,7,1,0,0,0):15/60/24:datetime(2023,8,31,0,0,0);

%%
load LISOct23_CH4N2O_CTD.mat;

surf_i = find(LISOct23_CH4N2O_CTD.Depth<3);

LIS_s = LISOct23_CH4N2O_CTD(surf_i,:);

dl = datetime(2023,9,1,0,0,0):15/60/24:datetime(2023,10,31,0,0,0);

%%
load LISMay24_CH4N2O_CTD.mat;

surf_i = find(LISMay24_CH4N2O_CTD.Depth<3);

LIS_s = LISMay24_CH4N2O_CTD(surf_i,:);

dl = datetime(2024,4,1,0,0,0):15/60/24:datetime(2024,5,31,0,0,0);


%%
% atmospheric concentrations averaged over the 3 cruises
N2Oatm_dry = 338e-9; %dry mole fraction
CH4atm_dry = 2020e-9; %dry mole fraction



% load in wind data
%load NDBC_44040h2023.mat;
% winds are collected at 3.5 m height and need to be adjusted to u10
load exrxMet_adj.mat;
% 
% mbar_to_atm = 1./1013.25; % convert pressure data from mbar to atm for use in gas sol eqns
% 
% exrxMet_adj.baroPress_atm = exrxMet_adj.baroPress_Avg .* mbar_to_atm;


wd = exrxMet_adj;
%%

% find the unique values in exrx
[unique_times_exrx, ~, idx] = unique(wd.TIMESTAMP);  % idx gives the group index for each key
num_keys = numel(unique_times_exrx);
avg_wspd_exrx = accumarray(idx(:), wd.windSpd_MS(:), [], @mean); % take average if there are repeats
avg_slp_exrx = accumarray(idx(:), wd.baroPress_atm(:), [], @mean); % take average if there are repeats

i_wspd_exrx = interp1(unique_times_exrx,avg_wspd_exrx,dl); % interpolate to fill in the gaps
i_slp_exrx = interp1(unique_times_exrx,avg_slp_exrx,dl); % interpolate to fill in the gaps

i_wspd_exrx = naninterp(i_wspd_exrx);
i_slp_exrx = naninterp(i_slp_exrx);

wspd_exrx_dl = nan(size(dl));


% if there are repeats, take the average of all values, otherwise just
% insert existing

for i = 1:length(dl)
    A = find(unique_times_exrx == dl(i));
    if ~isempty(A)
        wspd_exrx_dl(i) = avg_wspd_exrx(A);
    end;
end;

% get a list of the nan values
for i = 1:length(i_wspd_exrx)
    A = find(unique_times_exrx == dl(i));
    if ~isempty(A)
        wspd_exrx_dl(i) = avg_wspd_exrx(A);
    end;
end;



%%
%b = wd;

b.datenum = datenum(dl);
b.slp = i_slp_exrx; % slp
time_int_wind_hr = 0.25; % interval between wind speed observations = 0.25 hr

b.WSPD = i_wspd_exrx;
h_WSPD = 3.5; % wind speed height in meters
b.u10 = b.WSPD .*(10/h_WSPD)^0.11; % correcting to 10 m height

figure(1)
clf; hold on;
plot(datenum(wd.TIMESTAMP),wd.windSpd_MS,'o');
plot(b.datenum,b.WSPD);
datetick;
xlim([datenum(2024,4,1),datenum(2024,6,1)]);

%%

% slp_mbar = ERA5.slp_Pa./100;
% mbar_to_atm = 1./1013.25; % convert pressure data from mbar to atm for use in gas sol eqns
% slp_atm = slp_mbar .* mbar_to_atm;

%clear slp_mbar;

% rename the imported files
%raice = demo_AMSR2; 
%rawind = ERA5.wspd;
%raslp = ERA5.slp_atm;
%gd = demo_gas_data;

wt_t_30 = 30; % weighting time of 30 days
wt_t_15 = 15; % weighting time of 15 days

% we want the wind to have dimesions time x lat x lon
% but ERA5.wspd has dimensions lon x lat x time so we need to permute

% rawind.wind = permute(ERA5.wspd,[3,2,1]);
% rawind.lat = ERA5.lat;
% rawind.lon = ERA5.lon;
% rawind.datet = ERA5.datetime;
% rawind.daten = datenum(ERA5.datetime);

% we want the slp to have dimesions time x lat x lon
% but ERA5.slp has dimensions lon x lat x time so we need to permute
% raslp.slp = permute(slp_atm,[3,2,1]);
% raslp.lat = ERA5.lat;
% raslp.lon = ERA5.lon;
% raslp.datet = ERA5.datetime;
% raslp.daten = datenum(ERA5.datetime);

%clear slp_atm;


% datetimes
gd = table;
gd.station = LIS_s.Station;

gd.dt = LIS_s.datetime;

gd.mld = LIS_s.mld;

gd.lat = LIS_s.Lat;

gd.lon = LIS_s.Lon;

gd.P = LIS_s.P;

gd.S = LIS_s.S;

gd.T = LIS_s.T;

gd.depth = LIS_s.Depth;



gd.n2o_nmolkg = LIS_s.N2O_mean_nmolkg;
gd.ch4_nmolkg = LIS_s.CH4_mean_nmolkg;

gd.n2o_std_nmolkg = LIS_s.N2O_std_nmolkg;
gd.ch4_std_nmolkg = LIS_s.CH4_std_nmolkg;



% calculate n2o concentration with water vapor pressure at saturation
% c_H2Osat = c_dry * (1 - H2Opress)
N2Oatm_H2Osat = N2Oatm_dry .* (1 - vpress(gd.S,gd.T));
gd.n2o_eq_nmolkg = N2Osol(gd.S,gd.T,N2Oatm_H2Osat).*1000;

CH4atm_H2Osat = CH4atm_dry .* (1 - vpress(gd.S,gd.T));
gd.ch4_eq_nmolkg = CH4sol(gd.S,gd.T,CH4atm_H2Osat)'.*1000;


% save the gas data to new variable names
%s_time = datenum(gd.yyyy,gd.mm,gd.dd,gd.HH,gd.MM,gd.SS);
s_time = datenum(gd.dt);
%s_lat = gd.lat;
%s_lon = gd.lon;
s_mld = gd.mld;
s_T = gd.T;
s_S = gd.S;
s_P = gd.P;
s_ch4 = gd.ch4_nmolkg;
s_n2o = gd.n2o_nmolkg;
s_ch4_std = gd.ch4_std_nmolkg;
s_n2o_std = gd.n2o_std_nmolkg;
s_ch4_eq = gd.ch4_eq_nmolkg;
s_n2o_eq = gd.n2o_eq_nmolkg;
s_station = gd.station;
s_depth = gd.depth;


%%

%und.lat = s_lat; 
%und.lon= s_lon;
% %%
% %--- Wind matrix
% % Interpolate 3D wind matrix to the observations
% % need to use numeric date in the ndgrid function
% % here X_wind is date, Y_wind is lat, Z_wind is lon
% %und.wind = interp1(und)
% 
% [X_wind,Y_wind,Z_wind] = ndgrid(rawind.daten,rawind.lat,rawind.lon);
% V_wind = rawind.wind;
% Y_wind = double(Y_wind);
% Z_wind = double(Z_wind);
% % calculate instantaneous wind speeds for time of each sample
% und.wind = interpn(X_wind,Y_wind,Z_wind,V_wind,und.time,und.lat,und.lon); % need to turn the Y and Z into doubles
% 
% %--- SLP matrix
% % Interpolate 3D SLP matrix to the observations
% [X_slp,Y_slp,Z_slp] = ndgrid(raslp.daten,raslp.lat,raslp.lon);
% V_slp = raslp.slp;
% Y_slp = double(Y_slp);
% Z_slp = double(Z_slp);
% und.slp = interpn(X_slp,Y_slp,Z_slp,V_slp,und.time,und.lat,und.lon); %instantaneous slp for time of cruise

%%

% 15 day weighting

% create structure to save the interpolated data
und.time= s_time; 

% calculate instantaneous wind speeds and slp for time of each sample
und.wind = interpn(b.datenum,b.u10,und.time); 
und.slp = interpn(b.datenum,b.slp,und.time);

%-------------- instantaneous and 15-day weighted fluxes -----------------
%--- Create historical wind and slp matrices for the samples using the
% weighting period wt_t_15 (15 days) 
wt_t = wt_t_15;
int = time_int_wind_hr; % time interval for wind speed and slp observations
lag = 1:(wt_t *(24/int)); % # observations back in time (1 : # days * # observations per day)
windmat = nan(length(und.time),length(lag)+1); % empty matrix; note: the +1 is to hold the current (i.e. 0 hr lag) observation
t_back = nan(length(und.time),length(lag)+1);

slp_t = wt_t_15;
lag_slp = 1:(slp_t *(24/int));
t_back_slp = nan(length(und.time),length(lag)+1);
slpmat = nan(length(und.time),length(lag_slp)+1); %empty matrix; note: the +1 is to hold the current (i.e. 0 hr lag) observation


for kk=fliplr(lag)
    t_back(:,kk) = und.time - datenum(0,0,0,kk*int,0,0); %lag sample time backwards by kk*int hours (e.g. if kk = 120 and int = 6 hrs ==> lag backwards 720 hrs (30 days)
    wind_back = interpn(b.datenum,b.u10,t_back(:,kk)); %interpolate to lat/lon/time grid - wind speed lagged backwards, interpolated onto cruise track grid

    windmat(:,lag(end)-lag(kk)+1) = wind_back; %fill matrix
end

t_back(:,end) = und.time;
windmat(:,end) = und.wind; %add instantaneous wind
und.windmat = windmat;

for kk=fliplr(lag_slp)
   t_back_slp(:,kk) = und.time - datenum(0,0,0,kk*int,0,0); %lag sample time backwards by kk*int hours (e.g. if kk = 120 and int = 6 hrs ==> lag backwards 720 hrs (30 days)
   slp_back = interpn(b.datenum,b.slp,t_back(:,kk)); %interpolate to lat/lon/time grid - wind speed lagged backwards, interpolated onto cruise track grid   
   slpmat(:,lag(end)-lag(kk)+1) = slp_back; %fill matrix
end

slpmat(:,end) = und.slp; %add instantaneous slp
und.slpmat = slpmat;

% clear variables no longer needed
clear  X_wind Y_wind Z_wind V_wind ...
    X_slp Y_slp Z_slp V_slp ...
    kk t_back wind_back uw_wind

% SKIPPING THE ICE FRACTION CORRECTION SINCE THIS IS ICE-FREE
% using kw_weighting rather than kw_weighting_ice

%--- make matrices of T, S, MLD, and Schmidt number
% Assume constant backwards in time
Tmat = repmat(s_T,1,length(lag)+1);
Smat = repmat(s_S,1,length(lag)+1);
zMLmat = repmat(s_mld,1,length(lag)+1);

[~,Scmat_CH4] = gasmoldiff(Smat,Tmat,'CH4');
[~,Scmat_N2O] = gasmoldiff(Smat,Tmat,'N2O');

% calculate open water gas transfer velocity following Wanninkhof et al.
% (2014)
spd = 60*60*24; % seconds per day
k.CH4 = kgas(windmat,Scmat_CH4,'W14').* spd; % m/d
k.N2O = kgas(windmat,Scmat_N2O,'W14') .* spd; % m/d
%%
%ice = und.icemat;

slp_inst = slpmat(:,end); % instantaneous slp

% calculate instantaneous gas transfer velocity for ice-free condition
gd.k_inst_CH4 = k.CH4(:,1); % instantaneous CH4 gas tranfer velocity in m/d
gd.k_inst_N2O = k.N2O(:,1); % instantaneous N2O gas tranfer velocity in m/d

% instantenous sea-air gas flux corrected for slp
gd.F_CH4_inst = gd.k_inst_CH4.*(s_ch4-s_ch4_eq.*slp_inst)./1000.*sw_dens(s_S,s_T,s_P); % umol/m2/d
gd.F_N2O_inst = gd.k_inst_N2O.*(s_n2o-s_n2o_eq.*slp_inst)./1000.*sw_dens(s_S,s_T,s_P); % umol/m2/d

% calculate 15-day weighted gas transfer velocity corrected for ice cover
gd.k_wt_15_CH4 = nan.*s_time; % initialize variable
gd.k_wt_15_N2O = nan.*s_time; % initialize variable
%gd.wt_15 = nan.*s_time;

for kk = 1:length(Scmat_CH4(:,1))
    gd.k_wt_15_CH4(kk) = kw_weighting(k.CH4(kk,:), int/24, wt_t, zMLmat(kk,:)); % m/d
    gd.k_wt_15_N2O(kk) = kw_weighting(k.N2O(kk,:), int/24, wt_t, zMLmat(kk,:)); % m/d
end

clear kk

slp_15 = mean(slpmat,2); % 15-day average slp

%ice_inst = ice(:,end); % instantaneous ice fraction
%ice_15 = mean(ice,2); % 30-day average ice fraction

% calculate 30-day weighted sea-air flux of gases
% Units are umol/m2/d = m/d * umol/kg * kg/m3
gd.F_CH4_15 = gd.k_wt_15_CH4.*(s_ch4-s_ch4_eq.*slp_15)./1000.*sw_dens(s_S,s_T,s_P); % umol/m2/d
gd.F_N2O_15 = gd.k_wt_15_N2O.*(s_n2o-s_n2o_eq.*slp_15)./1000.*sw_dens(s_S,s_T,s_P); % umol/m2/d

gd.slp_inst = slp_inst;
gd.slp_15 = slp_15;

%gd.tres_15_N2O = gd.mld./gd.k_wt_15_N2O;
%gd.tres_15_CH4 = gd.mld./gd.k_wt_15_CH4;


disp('processed 15-day weighted gas fluxes')

%%
% 30 day weighting
% create structure to save the interpolated data
und.time= s_time; 

% calculate instantaneous wind speed and slp for time of each sample
und.wind = interpn(b.datenum,b.u10,und.time); 
und.slp = interpn(b.datenum,b.slp,und.time);


%-------------- instantaneous and 30-day weighted fluxes -----------------
%--- Create historical wind and slp matrices for the samples using the
% weighting period wt_t_30 (30 days) 
wt_t = wt_t_30;
int = time_int_wind_hr; % time interval for wind speed and slp observations
lag = 1:(wt_t *(24/int)); % # observations back in time (1 : # days * # observations per day)
windmat = nan(length(und.time),length(lag)+1); % empty matrix; note: the +1 is to hold the current (i.e. 0 hr lag) observation
t_back = nan(length(und.time),length(lag)+1);

slp_t = wt_t_30;
lag_slp = 1:(slp_t *(24/int));
t_back_slp = nan(length(und.time),length(lag)+1);
slpmat = nan(length(und.time),length(lag_slp)+1); %empty matrix; note: the +1 is to hold the current (i.e. 0 hr lag) observation


for kk=fliplr(lag)
    t_back(:,kk) = und.time - datenum(0,0,0,kk*int,0,0); %lag sample time backwards by kk*int hours (e.g. if kk = 120 and int = 6 hrs ==> lag backwards 720 hrs (30 days)
    wind_back = interpn(b.datenum,b.u10,t_back(:,kk)); %interpolate to lat/lon/time grid - wind speed lagged backwards, interpolated onto cruise track grid

    windmat(:,lag(end)-lag(kk)+1) = wind_back; %fill matrix
end

t_back(:,end) = und.time;
windmat(:,end) = und.wind; %add instantaneous wind
und.windmat = windmat;


% now calculate slp
und.slp = interpn(b.datenum,b.slp,und.time); 

for kk=fliplr(lag_slp)
    t_back_slp(:,kk) = und.time - datenum(0,0,0,kk*int,0,0); %lag sample time backwards by kk*int hours (e.g. if kk = 120 and int = 6 hrs ==> lag backwards 720 hrs (30 days)
    slp_back = interpn(b.datenum,b.slp,t_back(:,kk));
    slpmat(:,lag(end)-lag(kk)+1) = slp_back; %fill matrix
end

slpmat(:,end) = und.slp; %add instantaneous slp
und.slpmat = slpmat;


% slpmat = ones(size(und.windmat));
% und.slpmat = ones(size(und.windmat));

% clear variables no longer needed
clear  X_wind Y_wind Z_wind V_wind ...
    X_slp Y_slp Z_slp V_slp ...
    kk t_back wind_back uw_wind

% SKIPPING THE ICE FRACTION CORRECTION SINCE THIS IS ICE-FREE
% using kw_weighting rather than kw_weighting_ice

%--- make matrices of T, S, MLD, and Schmidt number
% Assume constant backwards in time
Tmat = repmat(s_T,1,length(lag)+1);
Smat = repmat(s_S,1,length(lag)+1);
zMLmat = repmat(s_mld,1,length(lag)+1);

[~,Scmat_CH4] = gasmoldiff(Smat,Tmat,'CH4');
[~,Scmat_N2O] = gasmoldiff(Smat,Tmat,'N2O');

% calculate open water gas transfer velocity following Wanninkhof et al.
% (2014)
spd = 60*60*24; % seconds per day
k.CH4 = kgas(windmat,Scmat_CH4,'W14').* spd; % m/d
k.N2O = kgas(windmat,Scmat_N2O,'W14') .* spd; % m/d

%ice = und.icemat;

slp_inst = slpmat(:,end); % instantaneous slp

% calculate instantaneous gas transfer velocity for ice-free condition
gd.k_inst_CH4 = k.CH4(:,1); % instantaneous CH4 gas tranfer velocity in m/d
gd.k_inst_N2O = k.N2O(:,1); % instantaneous N2O gas tranfer velocity in m/d

% instantenous sea-air gas flux corrected for slp
gd.F_CH4_inst = gd.k_inst_CH4.*(s_ch4-s_ch4_eq.*slp_inst)./1000.*sw_dens(s_S,s_T,s_P); % umol/m2/d
gd.F_N2O_inst = gd.k_inst_N2O.*(s_n2o-s_n2o_eq.*slp_inst)./1000.*sw_dens(s_S,s_T,s_P); % umol/m2/d

% calculate 30-day weighted gas transfer velocity corrected for ice cover
gd.k_wt_30_CH4 = nan.*s_time; % initialize variable
gd.k_wt_30_N2O = nan.*s_time; % initialize variable
%gd.wt_30 = nan.*s_time;


for kk = 1:length(Scmat_CH4(:,1))
    gd.k_wt_30_CH4(kk) = kw_weighting(k.CH4(kk,:), int/24, wt_t, zMLmat(kk,:)); % m/d
    gd.k_wt_30_N2O(kk) = kw_weighting(k.N2O(kk,:), int/24, wt_t, zMLmat(kk,:)); % m/d
end

% for kk = 1:length(Scmat_CH4(:,1))
%     [gd.k_wt_30_CH4(kk) gd.wt_30(kk,:)] = kw_weighting(k.CH4(kk,:), int/24, wt_t, zMLmat(kk,:)); % m/d
%     gd.k_wt_30_N2O(kk) = kw_weighting(k.N2O(kk,:), int/24, wt_t, zMLmat(kk,:)); % m/d
% end

clear kk

slp_30 = mean(slpmat,2); % 30-day average slp

%ice_inst = ice(:,end); % instantaneous ice fraction
%ice_30 = mean(ice,2); % 30-day average ice fraction

% calculate 30-day weighted sea-air flux of gases
% Units are umol/m2/d = m/d * umol/kg * kg/m3
gd.F_CH4_30 = gd.k_wt_30_CH4.*(s_ch4-s_ch4_eq.*slp_30)./1000.*sw_dens(s_S,s_T,s_P); % umol/m2/d
gd.F_N2O_30 = gd.k_wt_30_N2O.*(s_n2o-s_n2o_eq.*slp_30)./1000.*sw_dens(s_S,s_T,s_P); % umol/m2/d

disp('processed 30-day weighted gas fluxes')

% for Oct only
%gd.k_wt_30_CH4 = nan.*gd.k_wt_30_CH4;
%gd.k_wt_30_N2O = nan.*gd.k_wt_30_N2O;
%gd.F_CH4_30 = nan.*gd.k_wt_30_CH4;
%gd.F_N2O_30 = nan.*gd.k_wt_30_CH4;

%gd.tres_30_N2O = gd.mld./gd.k_wt_30_N2O;
%gd.tres_30_CH4 = gd.mld./gd.k_wt_30_CH4;


%%

LIS_gas_flux_Oct = gd; %rename variable for saving
save LIS_gas_flux_Oct.mat LIS_gas_flux_Oct;



%%
A=find(b.datenum>=datenum(2023,7,15,0,0,0)&b.datenum<datenum(2023,8,15,0,0,0));
k_JA = kgas(b.u10(A),Scmat_CH4(1),'W14').* spd; % m/d
k_avg.CH4 = mean(k_JA);
F_CH4_avg = k_avg.CH4.*(s_ch4-s_ch4_eq.*slp_30)./1000.*sw_dens(s_S,s_T,s_P); 

%%

%save LIS_N2O_flux.mat SBUS_N2O_flux; % save data

% CH4 flux is averaging 83 umol/m2/d
% N2O flux is averaging 1.41 umol/m2/d

% We know that CO2 in WLIS ranges from 750-1500 uatm in summer
% CO2 flux is around 13 mmol m-2 d-1
% based on:
% [F_CO2, dpCO2]=FCO2(1000,420,22,28,3)
% CO2 flux is 13 000 umol m-2 d-1
% CH4 flux is ~100 umol m-2 d-1 --> multiply by 10 to get CO2 equivalents
% for GWP = 1000 umol m-2 d-1
% at max, GWP = 3130


