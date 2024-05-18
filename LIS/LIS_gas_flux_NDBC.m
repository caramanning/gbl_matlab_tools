% this code calculates the gas flux based on NDBC mooring data collected at
% one location (mooring 44022, EXRX)
% data from mooring 44040, WLIS, was not available for the desired time
% period

load LISOct23_CH4N2O_CTD.mat;

surf_i = find(LISOct23_CH4N2O_CTD.Depth<3);

LIS_s = LISOct23_CH4N2O_CTD(surf_i,:);

%%

% load in wind data
%load NDBC_44040h2023.mat;
% winds are collected at 3.5 m height and need to be adjusted to u10
load NDBC_44022h2023.mat;
b = NDBC_44022h2023;

b.datenum = datenum(b.YY,b.MM,b.DD,b.hh,b.mm,0.*b.mm);

b.slp = ones(size(b.WSPD));
time_int_wind_hr = 0.25; % interval between wind speed observations = 0.25 hr

b.WSPD(b.WSPD>99)=nan;
h_WSPD = 3.5; % wind speed height in meters
b.u10 = b.WSPD .*(10/3.5)^0.11; % correcting to 10 m height

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
wt_t_60 = 60; % weighting time of 60 days

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

clear slp_atm;


% fake date times, the date is correct but times are fake
gd.dt = LIS_s.datetime;
gd.dn = datenum(gd.dt);
[gd.Y, gd.M, gd.D] = ymd(gd.dt);
[gd.H, gd.MI, gd.S] = hms(gd.dt);

% fake mld
gd.mld = 10.*ones(size(gd.dn));

gd.lat = LIS_s.Lat;

gd.lon = LIS_s.Lon;

gd.P = LIS_s.P;

gd.S = LIS_s.S;

gd.T = LIS_s.T;

gd.depth = LIS_s.Depth;

gd.station = LIS_s.Station;

gd.n2o_nmolkg = LIS_s.mean_N2O_nM./sw_dens(gd.S,gd.T,gd.P).*1000;
gd.ch4_nmolkg = LIS_s.mean_CH4_nM./sw_dens(gd.S,gd.T,gd.P).*1000;


% need to use the exact atmospheric concentration corrected for water vapor here, just approximating
% for now
N2Oatm = 330e-9;
gd.n2o_eq_nmolkg = N2Osol(gd.S,gd.T,N2Oatm).*1000;

CH4atm = 1920e-9;
gd.ch4_eq_nmolkg = CH4sol(gd.S,gd.T,CH4atm)'.*1000;


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
% create structure to save the interpolated data
und.time= s_time; 

% calculate instantaneous wind speeds for time of each sample
und.wind = interpn(b.datenum,b.u10,und.time); % need to turn the Y and Z into doubles

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
    wind_back = interpn(b.datenum,b.WSPD,t_back(:,kk)); %interpolate to lat/lon/time grid - wind speed lagged backwards, interpolated onto cruise track grid
    %wind_back = interpn(X_wind,Y_wind,Z_wind,V_wind,t_back(:,kk),und.lat,und.lon); %interpolate to lat/lon/time grid - wind speed lagged backwards, interpolated onto cruise track grid
    
    windmat(:,lag(end)-lag(kk)+1) = wind_back; %fill matrix
end
t_back(:,end) = und.time;
windmat(:,end) = und.wind; %add instantaneous wind
und.windmat = windmat;
%%
%for kk=fliplr(lag_slp)
%    t_back_slp(:,kk) = und.time - datenum(0,0,0,kk*int,0,0); %lag sample time backwards by kk*int hours (e.g. if kk = 120 and int = 6 hrs ==> lag backwards 720 hrs (30 days)
%    slp_back = interpn(X_slp,Y_slp,Z_slp,V_slp,t_back_slp(:,kk),und.lat,und.lon); %interpolate to lat/lon/time grid - wind speed lagged backwards, interpolated onto cruise track grid
%    slpmat(:,lag(end)-lag(kk)+1) = slp_back; %fill matrix
%end

%slpmat(:,end) = und.slp; %add instantaneous slp
%und.slpmat = slpmat;

% temporarily putting in 1 values for pressure
slpmat = ones(size(und.windmat));
und.slpmat = ones(size(und.windmat));

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

for kk = 1:length(Scmat_CH4(:,1))
    gd.k_wt_30_CH4(kk) = kw_weighting(k.CH4(kk,:), int/24, wt_t, zMLmat(kk,:)); % m/d
    gd.k_wt_30_N2O(kk) = kw_weighting(k.N2O(kk,:), int/24, wt_t, zMLmat(kk,:)); % m/d
end

clear kk

slp_30 = mean(slpmat,2); % 30-day average slp

%ice_inst = ice(:,end); % instantaneous ice fraction
%ice_30 = mean(ice,2); % 30-day average ice fraction

% calculate 30-day weighted sea-air flux of gases
% Units are umol/m2/d = m/d * umol/kg * kg/m3
gd.F_CH4_30 = gd.k_wt_30_CH4.*(s_ch4-s_ch4_eq.*slp_30)./1000.*sw_dens(s_S,s_T,s_P); % umol/m2/d
gd.F_N2O_30 = gd.k_wt_30_N2O.*(s_n2o-s_n2o_eq.*slp_30)./1000.*sw_dens(s_S,s_T,s_P); % umol/m2/d
disp('processed 30-day weighted gas fluxes')

LIS_gas_flux = gd; %rename variable for saving


%%
A=find(b.datenum>=datenum(2023,7,15,0,0,0)&b.datenum<datenum(2023,8,15,0,0,0));
k_JA = kgas(b.WSPD(A),Scmat_CH4(1),'W14').* spd; % m/d
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


