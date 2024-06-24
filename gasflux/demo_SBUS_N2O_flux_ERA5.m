
cd ../data
load IEPAug23_CH4N2O_CTD.mat
IEP = IEPAug23_CH4N2O_CTD;
%% THINGS THAT STILL NEED TO BE DONE:
%Shallow depth findong
%MLD finding
%exact N2O atmosphere and then N2Oeq
%Flux


% for this code I have done a demo with just three measurements and only
% done a 30-day weighting to demonstrate how to run it

% you need to extract the shallowest measurement from each profile and only
% calculate fluxes for those (probably depth < 6 m)
% save these to a new table that contains all of the associated CTD data
% for the surface samples only

% If code is above 6 m, then we will not get KML02. So, I took shallowest
% depth

load('IEPAug23_CH4N2O_CTD.mat');

IEPAug23_CH4N2O_CTD.Station= string(IEPAug23_CH4N2O_CTD.Station);

% Find indices of stations
nmlIndices = startsWith(IEPAug23_CH4N2O_CTD.Station, 'NML');
kmlIndices = startsWith(IEPAug23_CH4N2O_CTD.Station, 'KML');
sbmlIndices = startsWith(IEPAug23_CH4N2O_CTD.Station, 'SHBML');

% Combine indices
validIndices = nmlIndices | kmlIndices | sbmlIndices;

% Filter the data for only the valid indices
filteredData = IEPAug23_CH4N2O_CTD(validIndices, :);

% Find the shallowest depth for each unique station
uniqueStations = unique(filteredData.Station);
shallowestData = table();

for i = 1:numel(uniqueStations)
    stationData = filteredData(filteredData.Station == uniqueStations(i), :);
    [~, shallowestIndex] = min(stationData.Depth);
    shallowestData= [shallowestData; stationData(shallowestIndex, :)];
    
    shallowest_measurements = shallowestData; 
   
end
%%
cd ..\data

save ('shallowest_CTD_data.mat', 'shallowest_measurements');


%%
%Need to insert CombinedData_MLD to this shallowest_CTD_data.mat 

load shallowest_CTD_data.mat;
%%

load combinedData_MLD.mat;
Shallow = shallowest_measurements;

Shallow.MLD = nan.*Shallow.mean_CH4_NM;

Shallow.MLD = combinedData_MLD.MLD;
 
shallowest_measurements = Shallow; %Back to the original name

cd ../data 

save shallowest_CTD_data.mat shallowest_measurements;

%%

% you need to go back and get the exact date and time for each cast (year
% month day hours mins seconds)

% you need to calculate mixed layer depth using calcmld.m and save that
% with every profile. I suggest using the density based criterion but we
% should compare both temperature and density and see if they agree.

%calculating mld seperate for kml, nml, and shbml

%%
% you need to calculate the exact equilibrium N2O concentration based on
% atmospheric data aligned with cruise and water vapor pressure from the
% SST and SSS
% you can look at the historical record from Cape Town and Namibia and plot it and we
% can interpolate to estimate the concentration in August 2023
% https://gml.noaa.gov/dv/site/?stacode=CPT
% https://gml.noaa.gov/dv/site/?stacode=NMB
cd ..\gasflux
load AtmN2OCPT.mat
figure(100)
clf;
hold on;
h1=plot(AtmN2OCPT.YearMonth,AtmN2OCPT.N2O, 'LineWidth',2); % original data
xlabel('Year.Month');
ylabel('Atm N2O');

% Extrapolate to December 2023
% Create new years vector including 2023
new_years = [AtmN2OCPT.YearMonth; 2023.8];

% Fit a linear model to the data
p = polyfit(AtmN2OCPT.YearMonth, AtmN2OCPT.N2O, 1);

% Predict N2O levels for the new year
new_n2o_levels = polyval(p, new_years);

% Plot the extrapolated data
h2=plot(new_years, new_n2o_levels, '--r', 'LineWidth', 2); % extrapolated data

% Mark the extrapolated point
h3=plot(2023, new_n2o_levels(end), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

legend([h1,h2,h3],'Original Data', 'Extrapolated Data', 'Extrapolated Point');
hold off;

title('Atmospheric N2O Extrapolation to December 2023');
grid on;

%%
%Resulted concentration 335.558
%Ready to calculate the eq N2O and delta N2O
cd ..\data

load shallowest_CTD_data.mat
IEP = shallowest_measurements;
IEP.PTemp = sw_ptmp(IEP.S,IEP.T,IEP.P,0); %Made a column in IEP and calcualted the Ptemp

N2Oatm_dry = 335.558e-9;

% calculate n2o concentration with water vapor pressure at saturation
% N2O_H2Osat = N2O_dry * (1 - H2Opress)
IEP.n2o_atm_h2osat = N2Oatm_dry .* (1 - vpress(IEP.S,IEP.T));
 
% calculate equilibrium concentration
IEP.n2o_eq_nmolkg = N2Osol(IEP.S,IEP.PTemp,IEP.n2o_atm_h2osat).*1000;

IEP.Dn2o = (IEP.mean_N2O_nmolperkg - IEP.n2o_eq_nmolkg)./IEP.n2o_eq_nmolkg.*100;

shallowest_measurements = IEP;

save shallowest_CTD_data.mat shallowest_measurements;


%%

cd ..\gasflux

% you need to add in the code for 60 day weighting (available in the
% gas_toolbox but will need to be modified - see demo_weighted_gas_fluxes.m
% as well as the code below 


% load in ERA5 data
% I have loaded hourly data with every 0.25 deg resolution
info = ncinfo("ERA5_SBUS_2023.nc");

ERA5.timeEMCWF = ncread('ERA5_SBUS_2023.nc','time');

ERA5.datestr = datestr(double(ERA5.timeEMCWF)./24 + datetime('1900-01-01 0:0:0'));
ERA5.datetime = datetime(ERA5.datestr);
u10 = ncread('ERA5_SBUS_2023.nc','u10'); % lon x lat x time
v10 = ncread('ERA5_SBUS_2023.nc','v10'); % lon x lat x time
ERA5.u10 = u10;
ERA5.v10 = v10;
ERA5.wspd = sqrt(u10.^2+v10.^2);
ERA5.lon = ncread('ERA5_SBUS_2023.nc','longitude');
ERA5.lat = ncread('ERA5_SBUS_2023.nc','latitude');
ERA5.slp_Pa = ncread('ERA5_SBUS_2023.nc','msl'); %(1 hPa = 1 mb = 100 Pa).
ERA5.sst = ncread('ERA5_SBUS_2023.nc','sst');

time_int_wind_hr = 1; % interval between wind speed observations = 1 hr for wind and slp; in the gas toolbox demo int = 6 because the data source is different
%%
% to check the time interval can use
%time_int = datenum(ERA5.datetime(2) - ERA5.datetime(1)).*24;

clear u10 v10;


slp_mbar = ERA5.slp_Pa./100;
mbar_to_atm = 1./1013.25; % convert pressure data from mbar to atm for use in gas sol eqns
slp_atm = slp_mbar .* mbar_to_atm;

clear slp_mbar;

% rename the imported files
%raice = demo_AMSR2; 
%rawind = ERA5.wspd;
%raslp = ERA5.slp_atm;
%gd = demo_gas_data;

wt_t_30 = 30; % weighting time of 30 days
wt_t_60 = 60; % weighting time of 60 days

% we want the wind to have dimesions time x lat x lon
% but ERA5.wspd has dimensions lon x lat x time so we need to permute

rawind.wind = permute(ERA5.wspd,[3,2,1]);
rawind.lat = ERA5.lat;
rawind.lon = ERA5.lon;
rawind.datet = ERA5.datetime;
rawind.daten = datenum(ERA5.datetime);

% we want the slp to have dimesions time x lat x lon
% but ERA5.slp has dimensions lon x lat x time so we need to permute
raslp.slp = permute(slp_atm,[3,2,1]);
raslp.lat = ERA5.lat;
raslp.lon = ERA5.lon;
raslp.datet = ERA5.datetime;
raslp.daten = datenum(ERA5.datetime);

clear slp_atm;


%%
% fake date times, the date is correct but times are fake
%KML

cd ..\data
%Using Shallowest_measurements sheet to do further codes
Shallow = shallowest_measurements;

gd.dt = Shallow.dt;
[gd.Y, gd.M, gd.D] = ymd(gd.dt);
[gd.H, gd.MI, gd.S] = hms(gd.dt);

% fake mld
%gd.mld = [10 10 10]';
gd.mld = Shallow.MLD;

gd.lat = Shallow.Lat;

gd.lon = Shallow.Lon;

gd.n2o_nmolkg = Shallow.mean_N2O_nmolperkg;

gd.P = Shallow.P;

gd.S = Shallow.S;

gd.T = Shallow.T;

gd.depth = Shallow.Depth;

gd.station = Shallow.Station;

% need to use the exact atmospheric concentration corrected for water vapor here, just approximating
% for now
%N2Oatm = 335e-9;
gd.n2o_eq_nmolkg = Shallow.n2o_eq_nmolkg;

% save the gas data to new variable names
%s_time = datenum(gd.yyyy,gd.mm,gd.dd,gd.HH,gd.MM,gd.SS);
s_time = datenum(gd.dt);
s_lat = gd.lat;
s_lon = gd.lon;
s_mld = gd.mld;
s_T = gd.T;
s_S = gd.S;
s_P = gd.P;
%s_ch4 = gd.ch4_nmolkg;
s_n2o = gd.n2o_nmolkg;
%s_ch4_eq = gd.ch4_eq_nmolkg;
s_n2o_eq = gd.n2o_eq_nmolkg;
s_station = gd.station;
s_depth = gd.depth;

%%

% create structure to save the interpolated data
und.time=s_time; 
und.lat = s_lat; 
und.lon= s_lon;

%--- Wind matrix
% Interpolate 3D wind matrix to the observations
% need to use numeric date in the ndgrid function
% here X_wind is date, Y_wind is lat, Z_wind is lon
[X_wind,Y_wind,Z_wind] = ndgrid(rawind.daten,rawind.lat,rawind.lon);
V_wind = rawind.wind;
Y_wind = double(Y_wind);
Z_wind = double(Z_wind);
% calculate instantaneous wind speeds for time of each sample
und.wind = interpn(X_wind,Y_wind,Z_wind,V_wind,und.time,und.lat,und.lon); % need to turn the Y and Z into doubles

%--- SLP matrix
% Interpolate 3D SLP matrix to the observations
[X_slp,Y_slp,Z_slp] = ndgrid(raslp.daten,raslp.lat,raslp.lon);
V_slp = raslp.slp;
Y_slp = double(Y_slp);
Z_slp = double(Z_slp);
und.slp = interpn(X_slp,Y_slp,Z_slp,V_slp,und.time,und.lat,und.lon); %instantaneous slp for time of cruise
%%
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
    wind_back = interpn(X_wind,Y_wind,Z_wind,V_wind,t_back(:,kk),und.lat,und.lon); %interpolate to lat/lon/time grid - wind speed lagged backwards, interpolated onto cruise track grid
    windmat(:,lag(end)-lag(kk)+1) = wind_back; %fill matrix
end
t_back(:,end) = und.time;
und.windmat = windmat;
windmat(:,end) = und.wind; %add instantaneous wind

for kk=fliplr(lag_slp)
    t_back_slp(:,kk) = und.time - datenum(0,0,0,kk*int,0,0); %lag sample time backwards by kk*int hours (e.g. if kk = 120 and int = 6 hrs ==> lag backwards 720 hrs (30 days)
    slp_back = interpn(X_slp,Y_slp,Z_slp,V_slp,t_back_slp(:,kk),und.lat,und.lon); %interpolate to lat/lon/time grid - wind speed lagged backwards, interpolated onto cruise track grid
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
%k.CH4 = kgas(windmat,Scmat_CH4,'W14').* spd; % m/d
k.N2O = kgas(windmat,Scmat_N2O,'W14') .* spd; % m/d

%ice = und.icemat;

slp_inst = slpmat(:,end); % instantaneous slp

% calculate instantaneous gas transfer velocity for ice-free condition
%gd.k_inst_CH4 = k.CH4(:,1); % instantaneous CH4 gas tranfer velocity in m/d
gd.k_inst_N2O = k.N2O(:,1); % instantaneous N2O gas tranfer velocity in m/d

% instantenous sea-air gas flux corrected for slp
%gd.F_CH4_inst = gd.k_inst_CH4.*(s_ch4-s_ch4_eq.*slp_inst)./1000.*sw_dens(s_S,s_T,s_P); % umol/m2/d
gd.F_N2O_inst = gd.k_inst_N2O.*(s_n2o-s_n2o_eq.*slp_inst)./1000.*sw_dens(s_S,s_T,s_P); % umol/m2/d

% calculate 30-day weighted gas transfer velocity corrected for ice cover
%gd.k_wt_30_CH4 = nan.*s_time; % initialize variable
gd.k_wt_30_N2O = nan.*s_time; % initialize variable

for kk = 1:length(Scmat_CH4(:,1))
    %gd.k_wt_30_CH4(kk) = kw_weighting(k.CH4(kk,:), int/24, wt_t, zMLmat(kk,:)); % m/d
    gd.k_wt_30_N2O(kk) = kw_weighting(k.N2O(kk,:), int/24, wt_t, zMLmat(kk,:)); % m/d
end

clear kk

slp_30 = mean(slpmat,2); % 30-day average slp

%ice_inst = ice(:,end); % instantaneous ice fraction
%ice_30 = mean(ice,2); % 30-day average ice fraction

% calculate 30-day weighted sea-air flux of gases
% Units are umol/m2/d = m/d * umol/kg * kg/m3
%gd.F_CH4_30 = gd.k_wt_30_CH4.*(s_ch4-s_ch4_eq.*slp_30)./1000.*sw_dens(s_S,s_T,s_P); % umol/m2/d
gd.F_N2O_30 = gd.k_wt_30_N2O.*(s_n2o-s_n2o_eq.*slp_30)./1000.*sw_dens(s_S,s_T,s_P); % umol/m2/d
disp('processed 30-day weighted gas fluxes')

SBUS_N2O_flux = gd; %rename variable for saving

save SBUS_N2O_flux.mat SBUS_N2O_flux; % save data

%%


sortedData = sort(SBUS_N2O_flux.F_N2O_30);
Q1 = prctile(sortedData, 25);
Q3 = prctile(sortedData, 75);
median (sortedData)

disp(['The first quartile (Q1) is: ', num2str(Q1)])
disp(['The third quartile (Q3) is: ', num2str(Q3)])


