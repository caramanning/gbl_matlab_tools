
stnlist = ["MID4-cast01"
    "MID4-cast03"
    "MID4-cast05"
    "MID4-cast07"
    "MID4-cast09"
    "MID4-cast11"
    "MID4-cast13"
    "MID4-cast14"];


%%
%load LISAug23_CH4N2O_CTD.mat;
%LIS = LISAug23_CH4N2O_CTD;
%LIS.Depth(isnan(LIS.Depth)) = 20;

% 
% load LISOct23_CH4N2O_CTD.mat;
% LIS = LISOct23_CH4N2O_CTD;
% 
% load LISMay24_CH4N2O_CTD.mat;
% LIS = LISMay24_CH4N2O_CTD;
% 
% 
% UTC_to_local = -4/24;
% LIS.datetime_local = LIS.datetime + UTC_to_local;
% 
% 
% %LIS = LISO;
%% AUGUST
UTC_to_local = -4/24;
CH4airA = 1920e-9;
N2OairA = 335e-9;

load LISAug23_CH4N2O_CTD.mat
LISA = LISAug23_CH4N2O_CTD; % August
LISA.datetime_local = LISA.datetime + UTC_to_local;
LISA.dn_local = datenum(LISA.datetime_local);

LISA.CH4_mean_nmolkg = LISA.mean_CH4_nM./(1000+LISA.PDen).*1000;
LISA.N2O_mean_nmolkg = LISA.mean_N2O_nM./(1000+LISA.PDen).*1000;
LISA.CH4_std_nmolkg = LISA.std_CH4_nM./(1000+LISA.PDen).*1000;
LISA.N2O_std_nmolkg = LISA.std_N2O_nM./(1000+LISA.PDen).*1000;

LISA.DCH4_nmolkg = LISA.CH4_mean_nmolkg - CH4sol(LISA.S,LISA.T,CH4airA)'.*1000;
LISA.DN2O_nmolkg = LISA.N2O_mean_nmolkg - N2Osol(LISA.S,LISA.T,N2OairA).*1000;
LISA.DO2_umolkg = LISA.O2_umolkg - O2sol(LISA.S,LISA.T);

LISA.DCH4 = (LISA.DCH4_nmolkg./CH4sol(LISA.S,LISA.T,CH4airA)'.*1000 - 1).*100;
LISA.DN2O = (LISA.DN2O_nmolkg./N2Osol(LISA.S,LISA.T,N2OairA).*1000 - 1).* 100;
LISA.DO2 = (LISA.O2_umolkg./O2sol(LISA.S,LISA.T) - 1).*100;


load LISAug2023CastData.mat;
LISCDA.Cast = LISAug2023CastData.Cast;
LISCDA.Lat = LISAug2023CastData.Lat;
LISCDA.Lon = LISAug2023CastData.Lon;
LISCDA.Station = LISAug2023CastData.Station;
LISCDA.Dmax = sw_dpth(LISAug2023CastData.Pmax_dbar,LISCDA.Lat);
LISCDA.datetime_UTC = datetime(LISAug2023CastData.DateTime_UTC,'InputFormat',"MMM dd yyyy HH:mm:ss");
LISCDA.datetime_local = LISCDA.datetime_UTC + UTC_to_local;
LISCDA.dn_local = datenum(LISCDA.datetime_local);


% add in the station depth for all casts
LISA.StationDepth = NaN.*LISA.Lat;

for i = 1:numel(LISA.StationDepth)    
    A = find(LISCDA.Cast==LISA.CastNum(i));
    LISA.StationDepth(i) = LISCDA.Dmax(A);
end;

% get the cast numbers and station depth for MID4 casts
MID4castA = [];

for i = 1:numel(stnlist)
     q = find(LISCDA.Station==stnlist(i));
     MID4castA = [MID4castA; q];
end;

DmaxA = LISCDA.Dmax(MID4castA); % max depth for station in transect based on deep cast
dnlocalA = LISCDA.dn_local(MID4castA); % local datetime for MID4 casts

%% OCTOBER
UTC_to_local = -4/24;
CH4airO = 1920e-9;
N2OairO = 335e-9;

load LISOct23_CH4N2O_CTD.mat
LISO = LISOct23_CH4N2O_CTD; % Octust
LISO.datetime_local = LISO.datetime + UTC_to_local;
LISO.dn_local = datenum(LISO.datetime_local);

LISO.CH4_mean_nmolkg = LISO.mean_CH4_nM./(1000+LISO.PDen).*1000;
LISO.N2O_mean_nmolkg = LISO.mean_N2O_nM./(1000+LISO.PDen).*1000;
LISO.CH4_std_nmolkg = LISO.std_CH4_nM./(1000+LISO.PDen).*1000;
LISO.N2O_std_nmolkg = LISO.std_N2O_nM./(1000+LISO.PDen).*1000;

LISO.DCH4_nmolkg = LISO.CH4_mean_nmolkg - CH4sol(LISO.S,LISO.T,CH4airA)'.*1000;
LISO.DN2O_nmolkg = LISO.N2O_mean_nmolkg - N2Osol(LISO.S,LISO.T,N2OairA).*1000;
LISO.DO2_umolkg = LISO.O2_umolkg - O2sol(LISO.S,LISO.T);

LISO.DCH4 = (LISO.DCH4_nmolkg./CH4sol(LISO.S,LISO.T,CH4airA)'.*1000 - 1).*100;
LISO.DN2O = (LISO.DN2O_nmolkg./N2Osol(LISO.S,LISO.T,N2OairA).*1000 - 1).* 100;
LISO.DO2 = (LISO.O2_umolkg./O2sol(LISO.S,LISO.T) - 1).*100;


load LISOct2023CastData.mat;
LISCDO.Cast = LISOct2023CastData.Cast;
LISCDO.Lat = LISOct2023CastData.Lat;
LISCDO.Lon = LISOct2023CastData.Lon;
LISCDO.Station = LISOct2023CastData.Station;
LISCDO.Dmax = sw_dpth(LISOct2023CastData.Pmax_dbar,LISCDO.Lat);
LISCDO.datetime_UTC = datetime(LISOct2023CastData.DateTime_UTC,'InputFormat',"MMM dd yyyy HH:mm:ss");
LISCDO.datetime_local = LISCDO.datetime_UTC + UTC_to_local;
LISCDO.dn_local = datenum(LISCDO.datetime_local);


% add in the station depth for all casts
LISO.StationDepth = NaN.*LISO.Lat;

for i = 1:numel(LISO.StationDepth)    
    q = find(LISCDO.Cast==LISO.CastNum(i));
    LISO.StationDepth(i) = LISCDO.Dmax(q);
end;

% get the cast numbers and station depth for MID4 casts
MID4castO = [];

for i = 1:numel(stnlist)
     q = find(LISCDO.Station==stnlist(i));
     MID4castO = [MID4castO; q];
end;

DmaxO = LISCDO.Dmax(MID4castA); % max depth for station in transect based on deep cast
dnlocalO = LISCDO.dn_local(MID4castA); % local datetime for MID4 casts


%% MAY
UTC_to_local = -4/24;
CH4airM = 1920e-9;
N2OairM = 335e-9;

load LISMay24_CH4N2O_CTD.mat
LISM = LISMay24_CH4N2O_CTD; % August
LISM.datetime_local = LISM.datetime + UTC_to_local;
LISM.dn_local = datenum(LISM.datetime_local);

LISM.CH4_mean_nmolkg = LISM.mean_CH4_nM./(1000+LISM.PDen).*1000;
LISM.N2O_mean_nmolkg = LISM.mean_N2O_nM./(1000+LISM.PDen).*1000;
LISM.CH4_std_nmolkg = LISM.std_CH4_nM./(1000+LISM.PDen).*1000;
LISM.N2O_std_nmolkg = LISM.std_N2O_nM./(1000+LISM.PDen).*1000;

LISM.DCH4_nmolkg = LISM.CH4_mean_nmolkg - CH4sol(LISM.S,LISM.T,CH4airA)'.*1000;
LISM.DN2O_nmolkg = LISM.N2O_mean_nmolkg - N2Osol(LISM.S,LISM.T,N2OairA).*1000;
LISM.DO2_umolkg = LISM.O2_umolkg - O2sol(LISM.S,LISM.T);

LISM.DCH4 = (LISM.DCH4_nmolkg./CH4sol(LISM.S,LISM.T,CH4airA)'.*1000 - 1).*100;
LISM.DN2O = (LISM.DN2O_nmolkg./N2Osol(LISM.S,LISM.T,N2OairA).*1000 - 1).* 100;
LISM.DO2 = (LISM.O2_umolkg./O2sol(LISM.S,LISM.T) - 1).*100;

load LISMay2024CastData.mat;
LISCDM.Cast = LISMay2024CastData.Cast;
LISCDM.Lat = LISMay2024CastData.Lat;
LISCDM.Lon = LISMay2024CastData.Lon;
LISCDM.Station = LISMay2024CastData.Station;
LISCDM.Dmax = sw_dpth(LISMay2024CastData.Pmax_dbar,LISCDM.Lat);
LISCDM.datetime_UTC = datetime(LISMay2024CastData.DateTime_UTC,'InputFormat',"MMM dd yyyy HH:mm:ss");
LISCDM.datetime_local = LISCDM.datetime_UTC + UTC_to_local;
LISCDM.dn_local = datenum(LISCDM.datetime_local);


% add in the station depth for all casts
LISM.StationDepth = NaN.*LISM.Lat;

for i = 1:numel(LISM.StationDepth)    
    q = find(LISCDM.Cast==LISM.CastNum(i));
    LISM.StationDepth(i) = LISCDM.Dmax(q);
end;

% get the cast numbers and station depth for MID4 casts
MID4castM = [];

for i = 1:numel(stnlist)
     q = find(LISCDM.Station==stnlist(i));
     MID4castM = [MID4castM; q];
end;

DmaxM = LISCDM.Dmax(MID4castM); % max depth for station in transect based on deep cast
dnlocalM = LISCDM.dn_local(MID4castM); % local datetime for MID4 casts
%%
% AUGUST INTERPOLATION
% we need to make a grid that is evenly spaced so that all the casts are
% interpolated onto the same spacing
% additionally, we add on time at the start and the end to make the first
% and last profiles a bit easier to look at

dl = [0:0.1:20]'; % depth spacing for interpolation
x = [dnlocalA(1)-0.5/24 dnlocalA' dnlocalA(end)+0.5/24]; % times for x axis

%t_grid = repmat(ti_EDT,length(dl),1);
%t_grid = datenum(t_grid);

n_stn = numel(x);

dl_grid = repmat(dl,1,n_stn);
x_gridA = repmat(x,length(dl),1);

% AUGUST 2023
CH4iA = nan(numel(dl),n_stn); % make a blank grid for storing CH4
N2OiA = nan(numel(dl),n_stn); % make a blank grid for storing N2O
O2iA = nan(numel(dl),n_stn); % make a blank grid for storing O2
SiA = nan(numel(dl),n_stn); % make a blank grid for storing S
PDeniA = nan(numel(dl),n_stn); % make a blank grid for storing PDen
DCH4iA = nan(numel(dl),n_stn); % make a blank grid for storing CH4
DN2OiA = nan(numel(dl),n_stn); % make a blank grid for storing N2O
DO2iA = nan(numel(dl),n_stn); % make a blank grid for storing O2
tiA = repmat(datetime(0,0,0), 1, n_stn);
dl_gridA = repmat(dl,1,n_stn);
x_gridA = repmat(x,length(dl),1);


asA = []; % indices containing samples we want to use

for i = 1:numel(stnlist)    
    q = find(LISA.Station==stnlist(i));
    asA = [asA; q];

    tiA(i+1) = LISA.datetime(q(1));
    CH4iA(:,i+1) = interp1(LISA.Depth(q),LISA.mean_CH4_nM(q),dl);
    N2OiA(:,i+1) = interp1(LISA.Depth(q),LISA.mean_N2O_nM(q),dl);
    SiA(:,i+1) = interp1(LISA.Depth(q),LISA.S(q),dl);
    O2iA(:,i+1) = interp1(LISA.Depth(q),LISA.O2_umolkg(q),dl);
    PDeniA(:,i+1) = interp1(LISA.Depth(q),LISA.PDen(q),dl); 
    DCH4iA(:,i+1) = interp1(LISA.Depth(q),LISA.DCH4(q),dl);
    DN2OiA(:,i+1) = interp1(LISA.Depth(q),LISA.DN2O(q),dl);
    DO2iA(:,i+1) = interp1(LISA.Depth(q),LISA.DO2(q),dl);    

    % fill in the values at start and end with the closest non-NaN value
    nnC = find(~isnan(CH4iA(:,i+1)));
    CH4iA(1:min(nnC),i+1) = CH4iA(min(nnC),i+1);
    CH4iA(max(nnC):end,i+1) = CH4iA(max(nnC),i+1);

    nnC = find(~isnan(N2OiA(:,i+1)));
    N2OiA(1:min(nnC),i+1) = N2OiA(min(nnC),i+1);
    N2OiA(max(nnC):end,i+1) = N2OiA(max(nnC),i+1);

    nnC = find(~isnan(SiA(:,i+1)));
    SiA(1:min(nnC),i+1) = SiA(min(nnC),i+1);
    SiA(max(nnC):end,i+1) = SiA(max(nnC),i+1);

    nnC = find(~isnan(O2iA(:,i+1)));
    O2iA(1:min(nnC),i+1) = O2iA(min(nnC),i+1);
    O2iA(max(nnC):end,i+1) = O2iA(max(nnC),i+1);    

    nnC = find(~isnan(PDeniA(:,i+1)));
    PDeniA(1:min(nnC),i+1) = PDeniA(min(nnC),i+1);
    PDeniA(max(nnC):end,i+1) = PDeniA(max(nnC),i+1); 

    nnC = find(~isnan(DCH4iA(:,i+1)));
    DCH4iA(1:min(nnC),i+1) = DCH4iA(min(nnC),i+1);
    DCH4iA(max(nnC):end,i+1) = DCH4iA(max(nnC),i+1);

    nnC = find(~isnan(DN2OiA(:,i+1)));
    DN2OiA(1:min(nnC),i+1) = DN2OiA(min(nnC),i+1);
    DN2OiA(max(nnC):end,i+1) = DN2OiA(max(nnC),i+1);

    nnC = find(~isnan(DO2iA(:,i+1)));
    DO2iA(1:min(nnC),i+1) = DO2iA(min(nnC),i+1);
    DO2iA(max(nnC):end,i+1) = DO2iA(max(nnC),i+1);     
end;

dmin = 1.4; % minimum depth to plot 
CH4iA(dl_gridA<dmin) = NaN;
N2OiA(dl_gridA<dmin) = NaN;
SiA(dl_gridA<dmin) = NaN;
O2iA(dl_gridA<dmin) = NaN;
PDeniA(dl_gridA<dmin) = NaN;
DCH4iA(dl_gridA<dmin) = NaN;
DN2OiA(dl_gridA<dmin) = NaN;
DO2iA(dl_gridA<dmin) = NaN;


% now add in data at start and end 
tiA(1) = tiA(2);
tiA(end) = tiA(end-1);

CH4iA(:,1) = CH4iA(:,2);
CH4iA(:,end) = CH4iA(:,end-1);

N2OiA(:,1) = N2OiA(:,2);
N2OiA(:,end) = N2OiA(:,end-1);

SiA(:,1) = SiA(:,2);
SiA(:,end) = SiA(:,end-1);

O2iA(:,1) = O2iA(:,2);
O2iA(:,end) = O2iA(:,end-1);

PDeniA(:,1) = PDeniA(:,2);
PDeniA(:,end) = PDeniA(:,end-1);

DCH4iA(:,1) = DCH4iA(:,2);
DCH4iA(:,end) = DCH4iA(:,end-1);

DN2OiA(:,1) = DN2OiA(:,2);
DN2OiA(:,end) = DN2OiA(:,end-1);

DO2iA(:,1) = DO2iA(:,2);
DO2iA(:,end) = DO2iA(:,end-1);

%%
% OCTOBER INTERPOLATION
% we need to make a grid that is evenly spaced so that all the casts are
% interpolated onto the same spacing
% additionally, we add on time at the start and the end to make the first
% and last profiles a bit easier to look at

dl = [0:0.1:20]'; % depth spacing for interpolation
x = [dnlocalO(1)-0.5/24 dnlocalO' dnlocalO(end)+0.5/24]; % times for x axis

n_stn = numel(x);

dl_grid = repmat(dl,1,n_stn);
x_gridO = repmat(x,length(dl),1);

% OCTOBER 2023
CH4iO = nan(numel(dl),n_stn); % make a blank grid for storing CH4
N2OiO = nan(numel(dl),n_stn); % make a blank grid for storing N2O
O2iO = nan(numel(dl),n_stn); % make a blank grid for storing O2
SiO = nan(numel(dl),n_stn); % make a blank grid for storing S
PDeniO = nan(numel(dl),n_stn); % make a blank grid for storing PDen
DCH4iO = nan(numel(dl),n_stn); % make a blank grid for storing CH4
DN2OiO = nan(numel(dl),n_stn); % make a blank grid for storing N2O
DO2iO = nan(numel(dl),n_stn); % make a blank grid for storing O
tiO = repmat(datetime(0,0,0), 1, n_stn);
dl_gridO = repmat(dl,1,n_stn);
x_gridO = repmat(x,length(dl),1);
asO = []; % indices containing samples we want to use

for i = 1:numel(stnlist)    
    q = find(LISO.Station==stnlist(i));
    asO = [asO; q];
    tiO(i+1) = LISO.datetime(q(1));
    CH4iO(:,i+1) = interp1(LISO.Depth(q),LISO.mean_CH4_nM(q),dl);
    N2OiO(:,i+1) = interp1(LISO.Depth(q),LISO.mean_N2O_nM(q),dl);
    SiO(:,i+1) = interp1(LISO.Depth(q),LISO.S(q),dl);
    O2iO(:,i+1) = interp1(LISO.Depth(q),LISO.O2_umolkg(q),dl);
    PDeniO(:,i+1) = interp1(LISO.Depth(q),LISO.PDen(q),dl); 
    DCH4iO(:,i+1) = interp1(LISO.Depth(q),LISO.DCH4(q),dl);
    DN2OiO(:,i+1) = interp1(LISO.Depth(q),LISO.DN2O(q),dl);
    DO2iO(:,i+1) = interp1(LISO.Depth(q),LISO.DO2(q),dl);        

    % fill in the values at start and end with the closest non-NaN value
    nnC = find(~isnan(CH4iO(:,i+1)));
    CH4iO(1:min(nnC),i+1) = CH4iO(min(nnC),i+1);
    CH4iO(max(nnC):end,i+1) = CH4iO(max(nnC),i+1);

    nnC = find(~isnan(N2OiO(:,i+1)));
    N2OiO(1:min(nnC),i+1) = N2OiO(min(nnC),i+1);
    N2OiO(max(nnC):end,i+1) = N2OiO(max(nnC),i+1);

    nnC = find(~isnan(SiO(:,i+1)));
    SiO(1:min(nnC),i+1) = SiO(min(nnC),i+1);
    SiO(max(nnC):end,i+1) = SiO(max(nnC),i+1);

    nnC = find(~isnan(O2iO(:,i+1)));
    O2iO(1:min(nnC),i+1) = O2iO(min(nnC),i+1);
    O2iO(max(nnC):end,i+1) = O2iO(max(nnC),i+1);    

    nnC = find(~isnan(PDeniO(:,i+1)));
    PDeniO(1:min(nnC),i+1) = PDeniO(min(nnC),i+1);
    PDeniO(max(nnC):end,i+1) = PDeniO(max(nnC),i+1); 

    nnC = find(~isnan(DCH4iO(:,i+1)));
    DCH4iO(1:min(nnC),i+1) = DCH4iO(min(nnC),i+1);
    DCH4iO(max(nnC):end,i+1) = DCH4iO(max(nnC),i+1);

    nnC = find(~isnan(DN2OiO(:,i+1)));
    DN2OiO(1:min(nnC),i+1) = DN2OiO(min(nnC),i+1);
    DN2OiO(max(nnC):end,i+1) = DN2OiO(max(nnC),i+1);

    nnC = find(~isnan(DO2iO(:,i+1)));
    DO2iO(1:min(nnC),i+1) = DO2iO(min(nnC),i+1);
    DO2iO(max(nnC):end,i+1) = DO2iO(max(nnC),i+1);     
end;

dmin = 1.4; % minimum depth to plot 
CH4iO(dl_gridO<dmin) = NaN;
N2OiO(dl_gridO<dmin) = NaN;
SiO(dl_gridO<dmin) = NaN;
O2iO(dl_gridO<dmin) = NaN;
PDeniO(dl_gridO<dmin) = NaN;
DCH4iO(dl_gridO<dmin) = NaN;
DN2OiO(dl_gridO<dmin) = NaN;
DO2iO(dl_gridO<dmin) = NaN;


% now add in data at start and end 
tiO(1) = tiO(2);
tiO(end) = tiO(end-1);

CH4iO(:,1) = CH4iO(:,2);
CH4iO(:,end) = CH4iO(:,end-1);

N2OiO(:,1) = N2OiO(:,2);
N2OiO(:,end) = N2OiO(:,end-1);

SiO(:,1) = SiO(:,2);
SiO(:,end) = SiO(:,end-1);

O2iO(:,1) = O2iO(:,2);
O2iO(:,end) = O2iO(:,end-1);

PDeniO(:,1) = PDeniO(:,2);
PDeniO(:,end) = PDeniO(:,end-1);

DCH4iO(:,1) = DCH4iO(:,2);
DCH4iO(:,end) = DCH4iO(:,end-1);

DN2OiO(:,1) = DN2OiO(:,2);
DN2OiO(:,end) = DN2OiO(:,end-1);

DO2iO(:,1) = DO2iO(:,2);
DO2iO(:,end) = DO2iO(:,end-1);

%%
% MAY INTERPOLATION
% we need to make a grid that is evenly spaced so that all the casts are
% interpolated onto the same spacing
% additionally, we add on time at the start and the end to make the first
% and last profiles a bit easier to look at

dl = [0:0.1:20]'; % depth spacing for interpolation
x = [dnlocalM(1)-0.5/24 dnlocalM' dnlocalM(end)+0.5/24]; % times for x axis

n_stn = numel(x);

dl_grid = repmat(dl,1,n_stn);
x_gridM = repmat(x,length(dl),1);

% MAY 2024
CH4iM = nan(numel(dl),n_stn); % make a blank grid for storing CH4
N2OiM = nan(numel(dl),n_stn); % make a blank grid for storing N2M
O2iM = nan(numel(dl),n_stn); % make a blank grid for storing O2
SiM = nan(numel(dl),n_stn); % make a blank grid for storing S
PDeniM = nan(numel(dl),n_stn); % make a blank grid for storing PDen
DCH4iM = nan(numel(dl),n_stn); % make a blank grid for storing CH4
DN2OiM = nan(numel(dl),n_stn); % make a blank grid for storing N2O
DO2iM = nan(numel(dl),n_stn); % make a blank grid for storing O
tiM = repmat(datetime(0,0,0), 1, n_stn);
dl_gridM = repmat(dl,1,n_stn);
x_gridM = repmat(x,length(dl),1);
asM = []; % indices containing samples we want to use

dl_gridM = repmat(dl,1,n_stn);
x_gridM = repmat(x,length(dl),1);

for i = 1:numel(stnlist)    
    q = find(LISM.Station==stnlist(i));
    asM = [asM; q];

    tiM(i+1) = LISM.datetime(q(1));
    CH4iM(:,i+1) = interp1(LISM.Depth(q),LISM.mean_CH4_nM(q),dl);
    N2OiM(:,i+1) = interp1(LISM.Depth(q),LISM.mean_N2O_nM(q),dl);
    SiM(:,i+1) = interp1(LISM.Depth(q),LISM.S(q),dl);
    O2iM(:,i+1) = interp1(LISM.Depth(q),LISM.O2_umolkg(q),dl);
    PDeniM(:,i+1) = interp1(LISM.Depth(q),LISM.PDen(q),dl); 
    DCH4iM(:,i+1) = interp1(LISM.Depth(q),LISM.DCH4(q),dl);
    DN2OiM(:,i+1) = interp1(LISM.Depth(q),LISM.DN2O(q),dl);
    DO2iM(:,i+1) = interp1(LISM.Depth(q),LISM.DO2(q),dl);        

    % fill in the values at start and end with the closest non-NaN value
    nnC = find(~isnan(CH4iM(:,i+1)));
    CH4iM(1:min(nnC),i+1) = CH4iM(min(nnC),i+1);
    CH4iM(max(nnC):end,i+1) = CH4iM(max(nnC),i+1);

    nnC = find(~isnan(N2OiM(:,i+1)));
    N2OiM(1:min(nnC),i+1) = N2OiM(min(nnC),i+1);
    N2OiM(max(nnC):end,i+1) = N2OiM(max(nnC),i+1);

    nnC = find(~isnan(SiM(:,i+1)));
    SiM(1:min(nnC),i+1) = SiM(min(nnC),i+1);
    SiM(max(nnC):end,i+1) = SiM(max(nnC),i+1);

    nnC = find(~isnan(O2iM(:,i+1)));
    O2iM(1:min(nnC),i+1) = O2iM(min(nnC),i+1);
    O2iM(max(nnC):end,i+1) = O2iM(max(nnC),i+1);    

    nnC = find(~isnan(PDeniM(:,i+1)));
    PDeniM(1:min(nnC),i+1) = PDeniM(min(nnC),i+1);
    PDeniM(max(nnC):end,i+1) = PDeniM(max(nnC),i+1);

    nnC = find(~isnan(DCH4iM(:,i+1)));
    DCH4iM(1:min(nnC),i+1) = DCH4iM(min(nnC),i+1);
    DCH4iM(max(nnC):end,i+1) = DCH4iM(max(nnC),i+1);

    nnC = find(~isnan(DN2OiM(:,i+1)));
    DN2OiM(1:min(nnC),i+1) = DN2OiM(min(nnC),i+1);
    DN2OiM(max(nnC):end,i+1) = DN2OiM(max(nnC),i+1);

    nnC = find(~isnan(DO2iM(:,i+1)));
    DO2iM(1:min(nnC),i+1) = DO2iM(min(nnC),i+1);
    DO2iM(max(nnC):end,i+1) = DO2iM(max(nnC),i+1);     
end;


dmin = 1.4; % minimum depth to plot 
CH4iM(dl_gridM<dmin) = NaN;
N2OiM(dl_gridM<dmin) = NaN;
SiM(dl_gridM<dmin) = NaN;
O2iM(dl_gridM<dmin) = NaN;
PDeniM(dl_gridM<dmin) = NaN;
DCH4iM(dl_gridM<dmin) = NaN;
DN2OiM(dl_gridM<dmin) = NaN;
DO2iM(dl_gridM<dmin) = NaN;

% now add in data at start and end 
tiM(1) = tiM(2);
tiM(end) = tiM(end-1);

CH4iM(:,1) = CH4iM(:,2);
CH4iM(:,end) = CH4iM(:,end-1);

N2OiM(:,1) = N2OiM(:,2);
N2OiM(:,end) = N2OiM(:,end-1);

SiM(:,1) = SiM(:,2);
SiM(:,end) = SiM(:,end-1);

O2iM(:,1) = O2iM(:,2);
O2iM(:,end) = O2iM(:,end-1);

PDeniM(:,1) = PDeniM(:,2);
PDeniM(:,end) = PDeniM(:,end-1);

DCH4iM(:,1) = DCH4iM(:,2);
DCH4iM(:,end) = DCH4iM(:,end-1);

DN2OiM(:,1) = DN2OiM(:,2);
DN2OiM(:,end) = DN2OiM(:,end-1);

DO2iM(:,1) = DO2iM(:,2);
DO2iM(:,end) = DO2iM(:,end-1);

%%
% set values above the first sample to the value of the top sample
% sdl = find(dl>1.5,1);
% for i = 1:length(CH4i(1,:))
%     % find first non-NaN value
%     fnn = find(~isnan(CH4i(:,i)),1);
%     CH4i(sdl:fnn-1,i) = CH4i(fnn,i);
%     N2Oi(sdl:fnn-1,i) = N2Oi(fnn,i);
%     Si(sdl:fnn-1,i) = Si(fnn,i);
%     O2i(sdl:fnn-1,i) = O2i(fnn,i);
%     PDeni(sdl:fnn-1,i) = PDeni(fnn,i);
% end;
%%

%% PLOT THE DATA
% help for getting ranges
% Density
[min(LISO.PDen) max(LISO.PDen)]

[min(LISM.mean_N2O_nM) max(LISM.mean_N2O_nM)]

[min(min(N2OiM)) max(max(N2OiM))]

[min(min(CH4iM)) max(max(CH4iM))]

[min(min(O2iM)) max(max(O2iM))]
%%
disp('DN2O ranges A O M')
[min(min(DN2OiA)) max(max(DN2OiA))
min(min(DN2OiO)) max(max(DN2OiO))
min(min(DN2OiM)) max(max(DN2OiM))]

%%
disp('DCH4 ranges A O M')
[min(min(DCH4iA)) max(max(DCH4iA))
min(min(DCH4iO)) max(max(DCH4iO))
min(min(DCH4iM)) max(max(DCH4iM))]

%%
disp('DO2 ranges A O M')
[min(min(DO2iA)) max(max(DO2iA))
min(min(DO2iO)) max(max(DO2iO))
min(min(DO2iM)) max(max(DO2iM))]

%%
disp('PDen ranges A O M')
[min(min(PDeniA)) max(max(PDeniA))
min(min(PDeniO)) max(max(PDeniO))
min(min(PDeniM)) max(max(PDeniM))]

%%
% get mean, median, min, max for stations MID4
    asA = [];
    asO = [];
    asM = [];
    for i = 1:length(stnlist)
        q = find(LISA.Station==stnlist(i));
        asA = [asA; q];

        q = find(LISO.Station==stnlist(i));
        asO = [asO; q];

        q = find(LISM.Station==stnlist(i));
        asM = [asM; q];
    end;

disp('N2O ranges A O M')
[mean(LISA.N2O_mean_nmolkg(asA)) median(LISA.N2O_mean_nmolkg(asA)) min(LISA.N2O_mean_nmolkg(asA)) max(LISA.N2O_mean_nmolkg(asA))
mean(LISO.N2O_mean_nmolkg(asO)) median(LISO.N2O_mean_nmolkg(asO)) min(LISO.N2O_mean_nmolkg(asO)) max(LISO.N2O_mean_nmolkg(asO))
mean(LISM.N2O_mean_nmolkg(asM)) median(LISM.N2O_mean_nmolkg(asM)) min(LISM.N2O_mean_nmolkg(asM)) max(LISM.N2O_mean_nmolkg(asM))]


%%
disp('CH4 ranges A O M')
[mean(LISA.CH4_mean_nmolkg(asA)) median(LISA.CH4_mean_nmolkg(asA)) min(LISA.CH4_mean_nmolkg(asA)) max(LISA.CH4_mean_nmolkg(asA))
mean(LISO.CH4_mean_nmolkg(asO)) median(LISO.CH4_mean_nmolkg(asO)) min(LISO.CH4_mean_nmolkg(asO)) max(LISO.CH4_mean_nmolkg(asO))
mean(LISM.CH4_mean_nmolkg(asM)) median(LISM.CH4_mean_nmolkg(asM)) min(LISM.CH4_mean_nmolkg(asM)) max(LISM.CH4_mean_nmolkg(asM))]


%%
[mean(LISM.DN2O) median(LISM.DN2O) min(LISM.DN2O(asM)) max(LISM.DN2O)]

%%
[min(min(DN2OiO)) max(max(DN2OiO))
min(min(DN2OiM)) max(max(DN2OiM))]


%% SUBPLOT WITH
% ROW 1: DCH4
% ROW 2: DN2O
% ROW 3: DO2
% ROW 4: PDEN

nr = 4; % number of rows
nc = 3; % number of columns
lw = 2; % default line width
fs = 10; % default font size
ms = 3; % default marker size
yt = [0 10 20 30]; %y-axis ticks;
xtA = [datenum(2023,08,02,08,0,0), datenum(2023,08,02,12,0,0) datenum(2023,08,02,16,0,0) datenum(2023,08,02,20,0,0), datenum(2023,08,03,0,0,0), datenum(2023,08,03,4,0,0), datenum(2023,08,03,8,0,0)];
xtlA = ['08:00';'12:00';'16:00';'20:00';'00:00';'04:00';'08:00'];

xtO = [datenum(2023,10,19,08,0,0), datenum(2023,10,19,12,0,0) datenum(2023,10,19,16,0,0) datenum(2023,10,19,20,0,0), datenum(2023,10,20,0,0,0), datenum(2023,10,20,4,0,0), datenum(2023,10,20,8,0,0)];
xtlO = ['08:00';'12:00';'16:00';'20:00';'00:00';'04:00';'08:00'];

xtM = [datenum(2024,05,22,08,0,0), datenum(2024,05,22,12,0,0) datenum(2024,05,22,16,0,0) datenum(2024,05,22,20,0,0), datenum(2024,05,23,0,0,0), datenum(2024,05,23,4,0,0), datenum(2024,05,23,8,0,0)];
xtlM = ['08:00';'12:00';'16:00';'20:00';'00:00';'04:00';'08:00'];


yl = [0 20]; % y-axis limit in m
xlA = [min(min(x_gridA)) max(max(x_gridA))]; %x-axis limit in time
xlO = [min(min(x_gridO)) max(max(x_gridO))]; %x-axis limit in time
xlM = [min(min(x_gridM)) max(max(x_gridM))]; %x-axis limit in time
caDCH4 = [20 90]; % CH4 axis limits
clevelDCH4 = [caDCH4(1):1:caDCH4(2)]; %CH4 colorbar levels

caDN2O = [0 5]; % N2O axis limits
clevelDN2O = [caDN2O(1):0.02:caDN2O(2)]; % N2O colorbar levels

caDO2 = [-230 110];
clevelDO2 = [caDO2(1):1:caDO2(2)];

caPDen = [15.7 18.7]; % PDen axis limits
clevelPDen = [caPDen(1):0.02:caPDen(2)]; % PDen colorbar levels

fig=figure(24);
clf;
sp=tight_subplot(nr,nc,[.025 .025],[.08 .04],[.08 .04]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 6]);
set(gcf,'renderer','painters');
set(gcf,'GraphicsSmoothing','on');

% rows: Aug / Oct / May
% columns: CH4 / N2O / O2 / Density
% SUBPLOT 1
subplot(sp(1))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridA,dl_grid,DCH4iA,clevelDCH4,'edgecolor','none');
    plot(LISA.dn_local(asA),LISA.Depth(asA),'+k', 'linewidth', 2, 'markersize',ms)

    caxis(caDCH4)
    %c = colorbar('location','eastoutside');
    %c.Label.String = 'CH_4 (nmol/kg)';
    %c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xlA);
    set(gca,'ytick',yt);
    set(gca,'yticklabel',yt);
    set(gca,'xtick',xtA);
    %set(gca,'xticklabel',xtlA);
    %xlabel('Transect distance [km]');
    ylabel('Depth [m]');
    title('August');

    %add in plot of bathymetry
    dmA = [yl(2) DmaxA(1) DmaxA' DmaxA(end) yl(2)];
    kcA = [x_gridA(1,1) x_gridA(1,:) x_gridA(1,end)];
    fill(datenum(kcA),dmA,[0.5 0.5 0.5],'edgecolor','none');
    axis ij;

% SUBPLOT 2
subplot(sp(2))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridO,dl_grid,DCH4iO,clevelDCH4,'edgecolor','none');
    plot(LISO.dn_local(asO),LISO.Depth(asO),'+k', 'linewidth', 2, 'markersize',ms)

    %xlabel('Transect distance [km]');
    %ylabel('Depth [m]');
    caxis(caDCH4)
    %c = colorbar('location','eastoutside');
    %c.Label.String = 'CH_4 (nmol/kg)';
    %c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xlO);
    set(gca,'ytick',yt);
    %set(gca,'yticklabel',yt);
    set(gca,'xtick',xtO);
    %set(gca,'xticklabel',xtlM);
    title('October');

    %add in plot of bathymetry
    dmO = [yl(2) DmaxO(1) DmaxO' DmaxO(end) yl(2)];
    kcO = [x_gridO(1,1) x_gridO(1,:) x_gridO(1,end)];
    fill(kcO,dmO,[0.5 0.5 0.5],'edgecolor','none');
    axis ij;


% SUBPLOT 3
subplot(sp(3))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    clevel = [25:1:460];
    ca = [25 460]; %colorbar limits
    C = contourf(x_gridM,dl_grid,DCH4iM,clevelDCH4,'edgecolor','none');
    plot(LISM.dn_local(asM),LISM.Depth(asM),'+k', 'linewidth', 2, 'markersize',ms)

    %xlabel('Transect distance [km]');
    %ylabel('Depth [m]');
    caxis(caDCH4)
    c = colorbar('location','eastoutside');
    c.Label.String = 'CH_4 (nmol/kg)';
    c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xlM);
    set(gca,'ytick',yt);
    %set(gca,'yticklabel',yt);
    set(gca,'xtick',xtM);
    %set(gca,'xticklabel',xt);
    title(['May'])

    %add in plot of bathymetry
    dmM = [yl(2) DmaxM(1) DmaxM' DmaxM(end) yl(2)];
    kcM = [x_gridM(1,1) x_gridM(1,:) x_gridM(1,end)];
    fill(kcM,dmM,[0.5 0.5 0.5],'edgecolor','none');
    axis ij;

%SUBPLOT 4
subplot(sp(4))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridA,dl_grid,DN2OiA,clevelDN2O,'edgecolor','none');
    plot(LISA.dn_local(asA),LISA.Depth(asA),'+k', 'linewidth', 2, 'markersize',ms)

    caxis(caDN2O)
    %c = colorbar('location','eastoutside');
    %c.Label.String = 'N_2O (nmol/kg)';
    %c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    %xlim(xl);
    set(gca,'ytick',yt);
    set(gca,'yticklabel',yt);
    set(gca,'xtick',xtA);
    %set(gca,'xticklabel',xt);
    %xlabel('Transect distance [km]');
    ylabel('Depth [m]');
    %title('August');

    %add in plot of bathymetry
    dmA = [yl(2) DmaxA(1) DmaxA' DmaxA(end) yl(2)];
    kcA = [x_gridA(1,1) x_gridA(1,:) x_gridA(1,end)];
    fill(kcA,dmA,[0.5 0.5 0.5],'edgecolor','none');
    axis ij;


% SUBPLOT 5
subplot(sp(5))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridO,dl_grid,DN2OiO,clevelDN2O,'edgecolor','none');
    plot(LISO.dn_local(asO),LISO.Depth(asO),'+k', 'linewidth', 2, 'markersize',ms)

    %xlabel('Transect distance [km]');
    %ylabel('Depth [m]');
    caxis(caDN2O)
    %c = colorbar('location','eastoutside');
    %c.Label.String = 'CH_4 (nmol/kg)';
    %c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xlO);
    set(gca,'ytick',yt);
    %set(gca,'yticklabel',yt);
    set(gca,'xtick',xtO);
    %set(gca,'xticklabel',xt);

    %add in plot of bathymetry
    dmM = [yl(2) DmaxO(1) DmaxO' DmaxO(end) yl(2)];
    kcO = [x_gridO(1,1) x_gridO(1,:) x_gridO(1,end)];
    fill(kcO,dmO,[0.5 0.5 0.5],'edgecolor','none');
    axis ij;

% SUBPLOT 6
subplot(sp(6))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridM,dl_grid,DN2OiM,clevelDN2O,'edgecolor','none');
    plot(LISM.dn_local(asM),LISM.Depth(asM),'+k', 'linewidth', 2, 'markersize',ms)

    %xlabel('Transect distance [km]');
    %ylabel('Depth [m]');
    caxis(caDN2O)
    c = colorbar('location','eastoutside');
    c.Label.String = 'N_2O (nmol/kg)';
    c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xlM);
    set(gca,'ytick',yt);
    %set(gca,'yticklabel',yt);
    set(gca,'xtick',xtM);
    %set(gca,'xticklabel',xt);

    %add in plot of bathymetry
    dmM = [yl(2) DmaxM(1) DmaxM' DmaxM(end) yl(2)];
    kcM = [x_gridM(1,1) x_gridM(1,:) x_gridM(1,end)];
    fill(kcM,dmM,[0.5 0.5 0.5],'edgecolor','none');
    axis ij;

%SUBPLOT 7
subplot(sp(7))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridA,dl_grid,DO2iA,clevelDO2,'edgecolor','none');
    plot(LISA.dn_local(asA),LISA.Depth(asA),'+k', 'linewidth', 2, 'markersize',ms)

    caxis(caDO2)
    %c = colorbar('location','eastoutside');
    %c.Label.String = 'N_2O (nmol/kg)';
    %c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xlA);
    set(gca,'ytick',yt);
    set(gca,'yticklabel',yt);
    set(gca,'xtick',xtA);
    %set(gca,'xticklabel',xt);
    %xlabel('Transect distance [km]');
    ylabel('Depth [m]');
    %title('August');

    %add in plot of bathymetry
    dmA = [yl(2) DmaxA(1) DmaxA' DmaxA(end) yl(2)];
    kcA = [x_gridA(1,1) x_gridA(1,:) x_gridA(1,end)];
    fill(kcA,dmA,[0.5 0.5 0.5],'edgecolor','none');
    axis ij;


% SUBPLOT 8
subplot(sp(8))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridO,dl_grid,DO2iO,clevelDO2,'edgecolor','none');
    plot(LISO.dn_local(asO),LISO.Depth(asO),'+k', 'linewidth', 2, 'markersize',ms)

    %xlabel('Transect distance [km]');
    %ylabel('Depth [m]');
    caxis(caDO2)
    %c = colorbar('location','eastoutside');
    %c.Label.String = 'CH_4 (nmol/kg)';
    %c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xlO);
    set(gca,'ytick',yt);
    %set(gca,'yticklabel',yt);
    set(gca,'xtick',xtO);
    %set(gca,'xticklabel',xt);

    %add in plot of bathymetry
    dmO = [yl(2) DmaxO(1) DmaxO' DmaxO(end) yl(2)];
    kcO = [x_gridO(1,1) x_gridO(1,:) x_gridO(1,end)];
    fill(kcO,dmO,[0.5 0.5 0.5],'edgecolor','none');
    axis ij;

% SUBPLOT 9
subplot(sp(9))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridM,dl_grid,DO2iM,clevelDO2,'edgecolor','none');
    plot(LISM.dn_local(asM),LISM.Depth(asM),'+k', 'linewidth', 2, 'markersize',ms)

    %xlabel('Transect distance [km]');
    %ylabel('Depth [m]');
    caxis(caDO2)
    c = colorbar('location','eastoutside');
    c.Label.String = 'O_2 (\mumol/kg)';
    c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xlM);
    set(gca,'ytick',yt);
    %set(gca,'yticklabel',yt);
    set(gca,'xtick',xtM);
    %set(gca,'xticklabel',xt);

    %add in plot of bathymetry
    dmM = [yl(2) DmaxM(1) DmaxM' DmaxM(end) yl(2)];
    kcM = [x_gridM(1,1) x_gridM(1,:) x_gridM(1,end)];
    fill(kcM,dmM,[0.5 0.5 0.5],'edgecolor','none');
    axis ij;


%SUBPLOT 10
subplot(sp(10))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridA,dl_grid,PDeniA,clevelPDen,'edgecolor','none');
    plot(LISA.dn_local(asA),LISA.Depth(asA),'+k', 'linewidth', 2, 'markersize',ms)

    caxis(caPDen)
    %c = colorbar('location','eastoutside');
    %c.Label.String = 'N_2O (nmol/kg)';
    %c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xlA);
    set(gca,'ytick',yt);
    set(gca,'yticklabel',yt);
    set(gca,'xtick',xtA);
    set(gca,'xticklabel',xtlA);
    xlabel('local time [EST]');
    ylabel('Depth [m]');
    %title('August');

    %add in plot of bathymetry
    dmA = [yl(2) DmaxA(1) DmaxA' DmaxA(end) yl(2)];
    kcA = [x_gridA(1,1) x_gridA(1,:) x_gridA(1,end)];
    fill(kcA,dmA,[0.5 0.5 0.5],'edgecolor','none');
    axis ij;


% SUBPLOT 11
subplot(sp(11))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridO,dl_grid,PDeniO,clevelPDen,'edgecolor','none');
    plot(LISO.dn_local(asO),LISA.Depth(asO),'+k', 'linewidth', 2, 'markersize',ms)

    %xlabel('Transect distance [km]');
    %ylabel('Depth [m]');
    clim(caPDen)
    %c = colorbar('location','eastoutside');
    %c.Label.String = 'CH_4 (nmol/kg)';
    %c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xlO);
    set(gca,'ytick',yt);
    %set(gca,'yticklabel',yt);
    set(gca,'xtick',xtO);
    set(gca,'xticklabel',xtlO);
    xlabel('local time [EST]');

    %add in plot of bathymetry
    dmM = [yl(2) DmaxO(1) DmaxO' DmaxO(end) yl(2)];
    kcO = [x_gridO(1,1) x_gridO(1,:) x_gridO(1,end)];
    fill(kcO,dmO,[0.5 0.5 0.5],'edgecolor','none');
    axis ij;

% SUBPLOT 12
subplot(sp(12))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridM,dl_grid,PDeniM,clevelPDen,'edgecolor','none');
  %  C = surf(x_grid,dl_grid,PDeniM,clevelPDen,'facecolor','interp','edgecolor','interp');
  %  surf instead, with 'FaceColor','interp', 'EdgeColor','interp' and %view(0,90)’.  
    plot(LISM.dn_local(asM),LISM.Depth(asM),'+k', 'linewidth', 2, 'markersize',ms)

    %xlabel('Transect distance [km]');
    %ylabel('Depth [m]');
    clim(caPDen)
    c = colorbar('location','eastoutside');
    c.Label.String = '\sigma_{\theta} (kg/m^3)';
    c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xlM);
    set(gca,'ytick',yt);
    %set(gca,'yticklabel',yt);
    set(gca,'xtick',xtM);
    set(gca,'xticklabel',xtlM);
    xlabel('local time [EST]');

    %add in plot of bathymetry
    dmM = [yl(2) DmaxM(1) DmaxM' DmaxM(end) yl(2)];
    kcM = [x_gridM(1,1) x_gridM(1,:) x_gridM(1,end)];
    fill(kcM,dmM,[0.5 0.5 0.5],'edgecolor','none');
    axis ij;
    
    wysiwyg;

%print(gcf, '-dpng', '-r300', 'MID4_DCH4_DN2O_DO2_PDen.png');
%print(gcf,'-depsc','-vector','MID4_DCH4_DN2O_DO2_PDen.eps');
%epsclean('MID4_DCH4_DN2O_DO2_PDen.eps','MID4_DCH4_DN2O_DO2_PDen.eps');


%%

nr = 4; % number of rows
nc = 3; % number of columns
lw = 2; % default line width
fs = 10; % default font size
ms = 3; % default marker size
yt = [0 10 20 30]; %y-axis ticks;
xtA = [datenum(2023,08,02,08,0,0), datenum(2023,08,02,12,0,0) datenum(2023,08,02,16,0,0) datenum(2023,08,02,20,0,0), datenum(2023,08,03,0,0,0), datenum(2023,08,03,4,0,0), datenum(2023,08,03,8,0,0)];
xtlA = ['08:00';'12:00';'16:00';'20:00';'00:00';'04:00';'08:00'];

xtO = [datenum(2023,10,19,08,0,0), datenum(2023,10,19,12,0,0) datenum(2023,10,19,16,0,0) datenum(2023,10,19,20,0,0), datenum(2023,10,20,0,0,0), datenum(2023,10,20,4,0,0), datenum(2023,10,20,8,0,0)];
xtlO = ['08:00';'12:00';'16:00';'20:00';'00:00';'04:00';'08:00'];

xtM = [datenum(2024,05,22,08,0,0), datenum(2024,05,22,12,0,0) datenum(2024,05,22,16,0,0) datenum(2024,05,22,20,0,0), datenum(2024,05,23,0,0,0), datenum(2024,05,23,4,0,0), datenum(2024,05,23,8,0,0)];
xtlM = ['08:00';'12:00';'16:00';'20:00';'00:00';'04:00';'08:00'];


yl = [0 20]; % y-axis limit in m
xlA = [min(min(x_gridA)) max(max(x_gridA))]; %x-axis limit in time
xlO = [min(min(x_gridO)) max(max(x_gridO))]; %x-axis limit in time
xlM = [min(min(x_gridM)) max(max(x_gridM))]; %x-axis limit in time
caCH4 = [25 100]; % CH4 axis limits
clevelCH4 = [caCH4(1):1:caCH4(2)]; %CH4 colorbar levels

caN2O = [8.5 15]; % N2O axis limits
clevelN2O = [8.5:0.02:14]; % N2O colorbar levels

caO2 = [0 316];
clevelO2 = [caO2(1):1:caO2(2)];

caPDen = [15.7 18.8]; % PDen axis limits
clevelPDen = [caPDen(1):0.02:caPDen(2)]; % PDen colorbar levels

fig=figure(23);
clf;
sp=tight_subplot(nr,nc,[.025 .025],[.08 .04],[.08 .04]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 6]);
set(gcf,'renderer','painters');
set(gcf,'GraphicsSmoothing','on');

% rows: Aug / Oct / May
% columns: CH4 / N2O / O2 / Density
% SUBPLOT 1
subplot(sp(1))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridA,dl_grid,CH4iA,clevelCH4,'edgecolor','none');
    plot(LISA.dn_local(asA),LISA.Depth(asA),'+k', 'linewidth', 2, 'markersize',ms)

    caxis(caCH4)
    %c = colorbar('location','eastoutside');
    %c.Label.String = 'CH_4 (nmol/kg)';
    %c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xlA);
    set(gca,'ytick',yt);
    set(gca,'yticklabel',yt);
    set(gca,'xtick',xtA);
    %set(gca,'xticklabel',xtlA);
    %xlabel('Transect distance [km]');
    ylabel('Depth [m]');
    title('August');

    %add in plot of bathymetry
    dmA = [yl(2) DmaxA(1) DmaxA' DmaxA(end) yl(2)];
    kcA = [x_gridA(1,1) x_gridA(1,:) x_gridA(1,end)];
    fill(datenum(kcA),dmA,[0.5 0.5 0.5],'edgecolor','none');
    axis ij;

% SUBPLOT 2
subplot(sp(2))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridO,dl_grid,CH4iO,clevelCH4,'edgecolor','none');
    plot(LISO.dn_local(asO),LISO.Depth(asO),'+k', 'linewidth', 2, 'markersize',ms)

    %xlabel('Transect distance [km]');
    %ylabel('Depth [m]');
    caxis(caCH4)
    %c = colorbar('location','eastoutside');
    %c.Label.String = 'CH_4 (nmol/kg)';
    %c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xlO);
    set(gca,'ytick',yt);
    %set(gca,'yticklabel',yt);
    set(gca,'xtick',xtO);
    %set(gca,'xticklabel',xtlM);
    title('October');

    %add in plot of bathymetry
    dmO = [yl(2) DmaxO(1) DmaxO' DmaxO(end) yl(2)];
    kcO = [x_gridO(1,1) x_gridO(1,:) x_gridO(1,end)];
    fill(kcO,dmO,[0.5 0.5 0.5],'edgecolor','none');
    axis ij;


% SUBPLOT 3
subplot(sp(3))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    clevel = [25:1:460];
    ca = [25 460]; %colorbar limits
    C = contourf(x_gridM,dl_grid,CH4iM,clevelCH4,'edgecolor','none');
    plot(LISM.dn_local(asM),LISM.Depth(asM),'+k', 'linewidth', 2, 'markersize',ms)

    %xlabel('Transect distance [km]');
    %ylabel('Depth [m]');
    caxis(caCH4)
    c = colorbar('location','eastoutside');
    c.Label.String = 'CH_4 (nmol/kg)';
    c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xlM);
    set(gca,'ytick',yt);
    %set(gca,'yticklabel',yt);
    set(gca,'xtick',xtM);
    %set(gca,'xticklabel',xt);
    title(['May'])

    %add in plot of bathymetry
    dmM = [yl(2) DmaxM(1) DmaxM' DmaxM(end) yl(2)];
    kcM = [x_gridM(1,1) x_gridM(1,:) x_gridM(1,end)];
    fill(kcM,dmM,[0.5 0.5 0.5],'edgecolor','none');
    axis ij;

%SUBPLOT 4
subplot(sp(4))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridA,dl_grid,N2OiA,clevelN2O,'edgecolor','none');
    plot(LISA.dn_local(asA),LISA.Depth(asA),'+k', 'linewidth', 2, 'markersize',ms)

    caxis(caN2O)
    %c = colorbar('location','eastoutside');
    %c.Label.String = 'N_2O (nmol/kg)';
    %c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    %xlim(xl);
    set(gca,'ytick',yt);
    set(gca,'yticklabel',yt);
    set(gca,'xtick',xtA);
    %set(gca,'xticklabel',xt);
    %xlabel('Transect distance [km]');
    ylabel('Depth [m]');
    %title('August');

    %add in plot of bathymetry
    dmA = [yl(2) DmaxA(1) DmaxA' DmaxA(end) yl(2)];
    kcA = [x_gridA(1,1) x_gridA(1,:) x_gridA(1,end)];
    fill(kcA,dmA,[0.5 0.5 0.5],'edgecolor','none');
    axis ij;


% SUBPLOT 5
subplot(sp(5))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridO,dl_grid,N2OiO,clevelN2O,'edgecolor','none');
    plot(LISO.dn_local(asO),LISO.Depth(asO),'+k', 'linewidth', 2, 'markersize',ms)

    %xlabel('Transect distance [km]');
    %ylabel('Depth [m]');
    caxis(caN2O)
    %c = colorbar('location','eastoutside');
    %c.Label.String = 'CH_4 (nmol/kg)';
    %c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xlO);
    set(gca,'ytick',yt);
    %set(gca,'yticklabel',yt);
    set(gca,'xtick',xtO);
    %set(gca,'xticklabel',xt);

    %add in plot of bathymetry
    dmM = [yl(2) DmaxO(1) DmaxO' DmaxO(end) yl(2)];
    kcO = [x_gridO(1,1) x_gridO(1,:) x_gridO(1,end)];
    fill(kcO,dmO,[0.5 0.5 0.5],'edgecolor','none');
    axis ij;

% SUBPLOT 6
subplot(sp(6))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridM,dl_grid,N2OiM,clevelN2O,'edgecolor','none');
    plot(LISM.dn_local(asM),LISM.Depth(asM),'+k', 'linewidth', 2, 'markersize',ms)

    %xlabel('Transect distance [km]');
    %ylabel('Depth [m]');
    caxis(caN2O)
    c = colorbar('location','eastoutside');
    c.Label.String = 'N_2O (nmol/kg)';
    c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xlM);
    set(gca,'ytick',yt);
    %set(gca,'yticklabel',yt);
    set(gca,'xtick',xtM);
    %set(gca,'xticklabel',xt);

    %add in plot of bathymetry
    dmM = [yl(2) DmaxM(1) DmaxM' DmaxM(end) yl(2)];
    kcM = [x_gridM(1,1) x_gridM(1,:) x_gridM(1,end)];
    fill(kcM,dmM,[0.5 0.5 0.5],'edgecolor','none');
    axis ij;

%SUBPLOT 7
subplot(sp(7))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridA,dl_grid,O2iA,clevelO2,'edgecolor','none');
    plot(LISA.dn_local(asA),LISA.Depth(asA),'+k', 'linewidth', 2, 'markersize',ms)

    caxis(caO2)
    %c = colorbar('location','eastoutside');
    %c.Label.String = 'N_2O (nmol/kg)';
    %c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xlA);
    set(gca,'ytick',yt);
    set(gca,'yticklabel',yt);
    set(gca,'xtick',xtA);
    %set(gca,'xticklabel',xt);
    %xlabel('Transect distance [km]');
    ylabel('Depth [m]');
    %title('August');

    %add in plot of bathymetry
    dmA = [yl(2) DmaxA(1) DmaxA' DmaxA(end) yl(2)];
    kcA = [x_gridA(1,1) x_gridA(1,:) x_gridA(1,end)];
    fill(kcA,dmA,[0.5 0.5 0.5],'edgecolor','none');
    axis ij;


% SUBPLOT 8
subplot(sp(8))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridO,dl_grid,O2iO,clevelO2,'edgecolor','none');
    plot(LISO.dn_local(asO),LISO.Depth(asO),'+k', 'linewidth', 2, 'markersize',ms)

    %xlabel('Transect distance [km]');
    %ylabel('Depth [m]');
    caxis(caO2)
    %c = colorbar('location','eastoutside');
    %c.Label.String = 'CH_4 (nmol/kg)';
    %c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xlO);
    set(gca,'ytick',yt);
    %set(gca,'yticklabel',yt);
    set(gca,'xtick',xtO);
    %set(gca,'xticklabel',xt);

    %add in plot of bathymetry
    dmO = [yl(2) DmaxO(1) DmaxO' DmaxO(end) yl(2)];
    kcO = [x_gridO(1,1) x_gridO(1,:) x_gridO(1,end)];
    fill(kcO,dmO,[0.5 0.5 0.5],'edgecolor','none');
    axis ij;

% SUBPLOT 9
subplot(sp(9))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridM,dl_grid,O2iM,clevelO2,'edgecolor','none');
    plot(LISM.dn_local(asM),LISM.Depth(asM),'+k', 'linewidth', 2, 'markersize',ms)

    %xlabel('Transect distance [km]');
    %ylabel('Depth [m]');
    caxis(caO2)
    c = colorbar('location','eastoutside');
    c.Label.String = 'O_2 (\mumol/kg)';
    c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xlM);
    set(gca,'ytick',yt);
    %set(gca,'yticklabel',yt);
    set(gca,'xtick',xtM);
    %set(gca,'xticklabel',xt);

    %add in plot of bathymetry
    dmM = [yl(2) DmaxM(1) DmaxM' DmaxM(end) yl(2)];
    kcM = [x_gridM(1,1) x_gridM(1,:) x_gridM(1,end)];
    fill(kcM,dmM,[0.5 0.5 0.5],'edgecolor','none');
    axis ij;


%SUBPLOT 10
subplot(sp(10))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridA,dl_grid,PDeniA,clevelPDen,'edgecolor','none');
    plot(LISA.dn_local(asA),LISA.Depth(asA),'+k', 'linewidth', 2, 'markersize',ms)

    caxis(caPDen)
    %c = colorbar('location','eastoutside');
    %c.Label.String = 'N_2O (nmol/kg)';
    %c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xlA);
    set(gca,'ytick',yt);
    set(gca,'yticklabel',yt);
    set(gca,'xtick',xtA);
    set(gca,'xticklabel',xtlA);
    xlabel('local time [EST]');
    ylabel('Depth [m]');
    %title('August');

    %add in plot of bathymetry
    dmA = [yl(2) DmaxA(1) DmaxA' DmaxA(end) yl(2)];
    kcA = [x_gridA(1,1) x_gridA(1,:) x_gridA(1,end)];
    fill(kcA,dmA,[0.5 0.5 0.5],'edgecolor','none');
    axis ij;


% SUBPLOT 11
subplot(sp(11))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridO,dl_grid,PDeniO,clevelPDen,'edgecolor','none');
    plot(LISO.dn_local(asO),LISA.Depth(asO),'+k', 'linewidth', 2, 'markersize',ms)

    %xlabel('Transect distance [km]');
    %ylabel('Depth [m]');
    clim(caPDen)
    %c = colorbar('location','eastoutside');
    %c.Label.String = 'CH_4 (nmol/kg)';
    %c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xlO);
    set(gca,'ytick',yt);
    %set(gca,'yticklabel',yt);
    set(gca,'xtick',xtO);
    set(gca,'xticklabel',xtlO);
    xlabel('local time [EST]');

    %add in plot of bathymetry
    dmM = [yl(2) DmaxO(1) DmaxO' DmaxO(end) yl(2)];
    kcO = [x_gridO(1,1) x_gridO(1,:) x_gridO(1,end)];
    fill(kcO,dmO,[0.5 0.5 0.5],'edgecolor','none');
    axis ij;

% SUBPLOT 12
subplot(sp(12))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridM,dl_grid,PDeniM,clevelPDen,'edgecolor','none');
  %  C = surf(x_grid,dl_grid,PDeniM,clevelPDen,'facecolor','interp','edgecolor','interp');
  %  surf instead, with 'FaceColor','interp', 'EdgeColor','interp' and %view(0,90)’.  
    plot(LISM.dn_local(asM),LISM.Depth(asM),'+k', 'linewidth', 2, 'markersize',ms)

    %xlabel('Transect distance [km]');
    %ylabel('Depth [m]');
    clim(caPDen)
    c = colorbar('location','eastoutside');
    c.Label.String = '\sigma_{\theta} (kg/m^3)';
    c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xlM);
    set(gca,'ytick',yt);
    %set(gca,'yticklabel',yt);
    set(gca,'xtick',xtM);
    set(gca,'xticklabel',xtlM);
    xlabel('local time [EST]');

    %add in plot of bathymetry
    dmM = [yl(2) DmaxM(1) DmaxM' DmaxM(end) yl(2)];
    kcM = [x_gridM(1,1) x_gridM(1,:) x_gridM(1,end)];
    fill(kcM,dmM,[0.5 0.5 0.5],'edgecolor','none');
    axis ij;
    
    wysiwyg;

%print(gcf, '-dpng', '-r300', 'MID4_CH4_N2O_O2_PDen.png');
%print(gcf,'-depsc','-vector','MID4_CH4_N2O_O2_PDen.eps');
%epsclean('MID4_CH4_N2O_O2_PDen.eps','MID4_CH4_N2O_O2_PDen.eps');
%%

% for OCTOBER set to EDT
% for AUGUST set to EDT

ti_EDT = ti - 4/24;
%ti_EST = ti - 5/24;

t_grid = repmat(ti_EDT,length(dl),1);
t_grid = datenum(t_grid);

%t_grid = repmat(ti_EDT,length(dl),1);


ms=8;
figure(2)
clf; hold on;
box on;
set(gca,'tickdir','out');
%clevel = [20:1:70]; % Oct
clevel = [20:1:85]; %Aug

C = contourf(t_grid,dl_grid,CH4i,clevel,'edgecolor','none');
plot(datenum(LIS.datetime_local(as)),LIS.Depth(as),'.k','markersize',ms)
xlabel('time (EDT)');
ylabel('Depth (m)');
c = colorbar;
c.Label.String = 'CH_4 (nM)';
c.Label.FontSize = 16;
datetick;
axis ij;
%print -dpng -r300 plot_MID4_CH4_Aug.png;
% -dpng -r300 plot_MID4_CH4_Oct.png;


figure(3)
clf; hold on;
box on;
set(gca,'tickdir','out');
clevel = [8:0.1:12];
clevel = [8:0.1:14.5];
C = contourf(t_grid,dl_grid,N2Oi,clevel,'edgecolor','none');
plot(datenum(LIS.datetime_local(as)),LIS.Depth(as),'.k','markersize',ms)
xlabel('time (EDT)');
ylabel('Depth (m)');
c = colorbar;
c.Label.String = 'N_2O (nM)';
c.Label.FontSize = 16;
datetick;
axis ij;
%print -dpng -r300 plot_MID4_N2O_Aug.png;
%print -dpng -r300 plot_MID4_N2O_Oct.png;
%%
figure(4)
clf; hold on;
box on;
set(gca,'tickdir','out');
%clevel = [190:260]; %Oct
clevel = [15:230]; %Aug
C = contourf(t_grid,dl_grid,O2i,clevel,'edgecolor','none');
plot(datenum(LIS.datetime_local(as)),LIS.Depth(as),'.k','markersize',ms)
xlabel('time (EDT)');
ylabel('Depth (m)');
c = colorbar;
c.Label.String = 'O_2 (\mumol kg^{-1})';
c.Label.FontSize = 16;
datetick;
axis ij;
%print -dpng -r300 plot_MID4_O2_Aug.png;
%print -dpng -r300 plot_MID4_O2_Oct.png;

figure(5)
clf; hold on;
box on;
set(gca,'tickdir','out');
%clevel = [190:260]; %Oct
clevel = [26:0.01:27]; %Aug
C = contourf(t_grid,dl_grid,Si,clevel,'edgecolor','none');
plot(datenum(LIS.datetime_local(as)),LIS.Depth(as),'.k','markersize',ms)
xlabel('time (EDT)');
ylabel('Depth (m)');
c = colorbar;
c.Label.String = 'Salinity [pss]';
c.Label.FontSize = 16;
datetick;
axis ij;
%print -dpng -r300 plot_MID4_S_Aug.png;
%print -dpng -r300 plot_MID4_S_Oct.png;



%PDen
figure(6)
clf; hold on;
box on;
set(gca,'tickdir','out');
%clevel = [190:260]; %Oct
clevel = [15:0.01:20]; %Aug
C = contourf(t_grid,dl_grid,PDeni,clevel,'edgecolor','none');
plot(datenum(LIS.datetime_local(as)),LIS.Depth(as),'.k','markersize',ms)
xlabel('time (EDT)');
ylabel('Depth (m)');
c = colorbar;
c.Label.String = '\sigma_{\theta} [kg/m^3]';
c.Label.FontSize = 16;
datetick;
axis ij;
%print -dpng -r300 plot_MID4_PDen_Aug.png;
%print -dpng -r300 plot_MID4_PDen_Oct.png;


%%
% Now plot as vertical profiles
stnlist2 = unique(LIS.Station);

stnlist = stnlist2(1);
yl = [0 24];
figure(2)
clf; hold on;

%sp=tight_subplot(nr,nc,[0.12 .03],[.15 .06],[.1 .04]);
hold on; box on;
%set(gcf,'color','w');
set(gcf, 'PaperUnits', 'inches');
set(gcf,'renderer','painters');
set(gcf, 'PaperPosition', [0 0 12 3]);
set(gca,'linewidth',0.75);
    set(gcf,'GraphicsSmoothing','on')

s1=subplot(1,4,1);
hold on; box on;
set(gca,'fontsize',16);
set(gca,'tickdir','out');
set(gca,'linewidth',1)
s1.Position = [0.06 0.15 0.19 0.7];

for i = 1:length(stnlist)
    stn = stnlist(i);
    A = find(LIS.Station==stn);
    [~,B] = sort(LIS.Depth(A));
    errorbar(LIS.mean_CH4_nM(A(B)),LIS.Depth(A(B)),LIS.std_CH4_nM(A(B)),'horizontal','o-','linewidth',1.5);
end;
axis ij;
title([stnlist(1)]);
ylabel('Depth [m]');
xlabel('CH_4');

s2=subplot(1,4,2);
hold on; box on;
set(gca,'fontsize',16);
set(gca,'tickdir','out');
set(gca,'linewidth',1)
s2.Position = [0.28 0.15 0.19 0.7];
stn = 'EXRX-cast01';

for i = 1:length(stnlist)
    stn = stnlist(i);
A = find(LIS.Station==stn);
[~,B] = sort(LIS.Depth(A));
errorbar(LIS.mean_N2O_nM(A(B)),LIS.Depth(A(B)),LIS.std_N2O_nM(A(B)),'horizontal','o-','linewidth',1.5);
end;
axis ij;
title('N_2O');
set(gca,'YTickLabel',[]);

s3=subplot(1,4,3);
hold on; box on;
set(gca,'fontsize',16);
set(gca,'tickdir','out');
set(gca,'linewidth',1)
s3.Position = [0.5 0.15 0.19 0.7];

for i = 1:length(stnlist)
    stn = stnlist(i);
    A = find(LIS.Station==stn);
    [~,B] = sort(LIS.Depth(A));
    plot(LIS.O2_umolkg(A(B)),LIS.Depth(A(B)),'o-','linewidth',1.5);
end;
axis ij;
title('O_2');
set(gca,'YTickLabel',[]);

s4=subplot(1,4,4);
hold on; box on;
set(gca,'fontsize',16);
set(gca,'tickdir','out');
set(gca,'linewidth',1)
s4.Position = [0.72 0.15 0.19 0.7];

for i = 1:length(stnlist)
    stn = stnlist(i);
    A = find(LIS.Station==stn);
    [~,B] = sort(LIS.Depth(A));
    plot(LIS.S(A(B)),LIS.Depth(A(B)),'o-','linewidth',1.5);
end;
axis ij;
title('S');
set(gca,'YTickLabel',[]);

%print -dpng -r300 MID4_Aug23_hourly_profiles.png;
%print -dpng -r300 MID4_Oct23_hourly_profiles.png;


%%

% 3 mg/L * 1 mol/32000 mg * 1e6 umol/mol

% plot O2, S, Density
stnlist2 = unique(LIS.Station);

stnlist = stnlist2(12);
yl = [0 24];
figure(2)
clf; hold on;

%sp=tight_subplot(nr,nc,[0.12 .03],[.15 .06],[.1 .04]);
hold on; box on;
%set(gcf,'color','w');
set(gcf, 'PaperUnits', 'inches');
set(gcf,'renderer','painters');
set(gcf, 'PaperPosition', [0 0 10 6]);
set(gca,'linewidth',0.75);
    set(gcf,'GraphicsSmoothing','on')


s1=subplot(1,3,1);
hold on; box on;
set(gca,'fontsize',16);
set(gca,'tickdir','out');
set(gca,'linewidth',1)
s1.Position = [0.1 0.25 0.25 0.65];

for i = 1:length(stnlist)
    stn = stnlist(i);
    A = find(LIS.Station==stn);
    [~,B] = sort(LIS.Depth(A));
    B(numel(B)) = []; % remove the deepest point
    plot(LIS.O2_umolkg(A(B)),LIS.Depth(A(B)),'o-','linewidth',1.5);
end;

%plot([90 90], [min(LIS.Depth(A(B))),max(LIS.Depth(A(B)))],'--k');
axis ij;
title([stnlist(1)]);
xlabel('O_2 [\mumol/kg]')
%set(gca,'YTickLabel',[]);
ylabel('Depth [m]')
set(gca,'xlim',[0 240]);
set(gca,'ylim',yl);

s2=subplot(1,3,2);
hold on; box on;
set(gca,'fontsize',16);
set(gca,'tickdir','out');
set(gca,'linewidth',1)
s2.Position = [0.4 0.25 0.25 0.65];

for i = 1:length(stnlist)
    stn = stnlist(i);
    A = find(LIS.Station==stn);
    [~,B] = sort(LIS.Depth(A));
    B(numel(B)) = []; % remove the deepest
    plot(LIS.S(A(B)),LIS.Depth(A(B)),'o-','linewidth',1.5);
end;
axis ij;
xlabel('Salinity [pss]');
set(gca,'YTickLabel',[]);
set(gca,'xlim',[26 27.2]);
set(gca,'ylim',yl);

% s3=subplot(1,3,3);
% hold on; box on;
% set(gca,'fontsize',16);
% set(gca,'tickdir','out');
% set(gca,'linewidth',1)
% s3.Position = [0.7 0.25 0.25 0.65];
% 
% for i = 1:length(stnlist)
%     stn = stnlist(i);
%     A = find(LIS.Station==stn);
%     [~,B] = sort(LIS.Depth(A));
%     plot(LIS.PDen(A(B)),LIS.Depth(A(B)),'o-','linewidth',1.5);
% end;
% axis ij;
% xlabel('\sigma_{\theta} [kg/m^3]');
% set(gca,'YTickLabel',[]);
% %set(gca,'xlim',[26 27]);

s3=subplot(1,3,3);
hold on; box on;
set(gca,'fontsize',16);
set(gca,'tickdir','out');
set(gca,'linewidth',1)
s3.Position = [0.7 0.25 0.25 0.65];

for i = 1:length(stnlist)
    stn = stnlist(i);
    A = find(LIS.Station==stn);
    [~,B] = sort(LIS.Depth(A));
    B(numel(B)) = []; % remove the deepest
    plot(LIS.T(A(B)),LIS.Depth(A(B)),'o-','linewidth',1.5);
end;
axis ij;
xlabel('Temperature [^oC]');
set(gca,'YTickLabel',[]);
set(gca,'xlim',[19.5 24]);
set(gca,'ylim',yl);

wysiwyg;

print -dpng -r300 profile_Aug_MID4-cast14.png;
%print -dpng -r300 MID4_Oct23_hourly_profiles.png;


%%

stn = stnlist(2);
A = find(LIS.Station==stn);
[~,B] = sort(LIS.Depth(A));
plot(LIS.T(A(B)),LIS.Depth(A(B)),'s-r','linewidth',1.5, 'markerfacecolor','r');

stn = stnlist(2);
A = find(LIS.Station==stn);
[~,B] = sort(LIS.Depth(A));
plot(LIS.T(A(B)),LIS.Depth(A(B)),'s-r','linewidth',1.5, 'markerfacecolor','r');
