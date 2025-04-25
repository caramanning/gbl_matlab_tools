% load in August 2023 dataOct23
load LISAug23_CH4N2O_CTD.mat
LISA = LISAug23_CH4N2O_CTD; % August

load LISAug2023CastData.mat;
LISCDA.Cast = LISAug2023CastData.Cast;
LISCDA.Lat = LISAug2023CastData.Lat;
LISCDA.Lon = LISAug2023CastData.Lon;
LISCDA.Dmax = sw_dpth(LISAug2023CastData.Pmax_dbar,LISCDA.Lat);

LISA.StationDepth = NaN.*LISA.Lat;

for i = 1:numel(LISA.StationDepth)    
    A = find(LISCDA.Cast==LISA.CastNum(i));
    LISA.StationDepth(i) = LISCDA.Dmax(A);
end;

ncast = 16:22;
km_between = m_lldist(LISCDA.Lon(ncast),LISCDA.Lat(ncast)); % distance between consecutive stations
km_between = [0; km_between]; % add on zero for first station

km_cumulativeA = cumsum(km_between); % consecutive distance
DmaxA = LISCDA.Dmax(ncast); % max depth for station in transect based on deep cast

% rename stations to be consistent
% stn = 'EXCR-cast01';
% A = find(LISA.Station==stn);
% LISA.Station(A) = 'EXRX-cast01';
% 
% stn = 'EXCR-cast02';
% A = find(LISA.Station==stn);
% LISA.Station(A) = 'EXRX-cast02';
% 
% LISAug23_CH4N2O_CTD = LISA;
% save LISAug23_CH4N2O_CTD.mat LISAug23_CH4N2O_CTD;

stnlist = ["EXR1-cast01"
    "EXRX-cast02"
    "MID3-cast01"
    "MID4-cast14"
    "MID5-cast01"
    "WLIS-cast02"
    "WLI6-cast01"]; % this matches cast 16-22 and x as listed above


    asA = [];
    for i = 1:length(stnlist)
        A = find(LISA.Station==stnlist(i));
        asA = [asA; A];
    end;

    CH4A = LISA.CH4_mean_nmolkg(asA);
    N2OA = LISA.N2O_mean_nmolkg(asA); 
    O2A = LISA.O2_umolkg(asA);
    DCH4A = LISA.DCH4(asA);
    DN2OA = LISA.DN2O(asA);
    DO2A = LISA.DO2(asA);
        
    TA = LISA.T(asA);

    [median(CH4A) prctile(CH4A,25) prctile(CH4A,75)]
    [mean(CH4A) std(CH4A)]

    [median(DCH4A) prctile(DCH4A,25) prctile(DCH4A,75)]
    [mean(DCH4A) std(DCH4A)]

    [median(N2OA) prctile(N2OA,25) prctile(N2OA,75)]
    [mean(N2OA) std(N2OA)]

    [median(DN2OA) prctile(DN2OA,25) prctile(DN2OA,75)]
    [mean(DN2OA) std(DN2OA)]

    [median(O2A) prctile(O2A,25) prctile(O2A,75)]
    [mean(O2A) std(O2A)]

    [median(DO2A) prctile(DO2A,25) prctile(DO2A,75)]
    [mean(DO2A) std(DO2A)]

%    [min(O2A) max(O2A) mean(O2A) median(O2A)]
%    [min(DCH4A) max(DCH4A) mean(DCH4A) median(DCH4A)]

%%
% load in October 2023 data
load LISOct23_CH4N2O_CTD.mat
LISO = LISOct23_CH4N2O_CTD; % October

load LISOct2023CastData.mat;
LISCDO.Cast = LISOct2023CastData.Cast;
LISCDO.Lat = LISOct2023CastData.Lat;
LISCDO.Lon = LISOct2023CastData.Lon;
LISCDO.Dmax = sw_dpth(LISOct2023CastData.Pmax_dbar,LISCDO.Lat);


ncast = 16:22;
Cast = LISO.Cast(ncast);
Lat = LISO.Lat(ncast);
Lon = LISO.Lon(ncast);

ncast = 16:22; % casts for transect;
km_between = m_lldist(LISCDO.Lon(ncast),LISCDO.Lat(ncast)); % distance between consecutive stations
km_between = [0; km_between]; % add on zero for first station

km_cumulativeO = cumsum(km_between); % consecutive distance
DmaxO = LISCDO.Dmax(ncast); % max depth for station in transect based on deep cast


    asA = [];
    for i = 1:length(stnlist)
        A = find(LISO.Station==stnlist(i));
        asA = [asA; A];
    end;

    CH4O = LISO.CH4_mean_nmolkg(asA);
    N2OO = LISO.N2O_mean_nmolkg(asA); 
    O2O = LISO.O2_umolkg(asA);
    DCH4O = LISO.DCH4(asA);
    DN2OO = LISO.DN2O(asA);
    DO2O = LISO.DO2(asA);

    [median(CH4O) prctile(CH4O,25) prctile(CH4O,75)]
    [mean(CH4O) std(CH4O)]

    [median(DCH4O) prctile(DCH4O,25) prctile(DCH4O,75)]
    [mean(DCH4O) std(DCH4O)]

    [median(N2OO) prctile(N2OO,25) prctile(N2OO,75)]
    [mean(N2OO) std(N2OO)]

    [median(DN2OO) prctile(DN2OO,25) prctile(DN2OO,75)]
    [mean(DN2OO) std(DN2OO)]

    [median(O2O) prctile(O2O,25) prctile(O2O,75)]
    [mean(O2O) std(O2O)]

    [median(DO2O) prctile(DO2O,25) prctile(DO2O,75)]
    [mean(DO2O) std(DO2O)]



% rename stations to be consistent
% stn = 'EXCR1-cast01';
% A = find(LISO.Station==stn);
% LISO.Station(A) = 'EXR1-cast01';
% 
% LISOct23_CH4N2O_CTD = LISO;
% save LISOct23_CH4N2O_CTD.mat LISOct23_CH4N2O_CTD;

%%

% load in May 2024 data
load LISMay24_CH4N2O_CTD.mat
LISM = LISMay24_CH4N2O_CTD; % May

load LISMay2024CastData.mat;
LISCDM.Cast = LISMay2024CastData.Cast;
LISCDM.Lat = LISMay2024CastData.Lat;
LISCDM.Lon = LISMay2024CastData.Lon;
LISCDM.Dmax = sw_dpth(LISMay2024CastData.Pmax_dbar,LISCDM.Lat);


ncast = 16:22;
Cast = LISO.Cast(ncast);
Lat = LISO.Lat(ncast);
Lon = LISO.Lon(ncast);

ncast = 16:22; % casts for transect;
km_between = m_lldist(LISCDO.Lon(ncast),LISCDO.Lat(ncast)); % distance between consecutive stations
km_between = [0; km_between]; % add on zero for first station

km_cumulativeM = cumsum(km_between); % consecutive distance
DmaxM = LISCDM.Dmax(ncast); % max depth for station in transect based on deep cast

    asA = [];
    for i = 1:length(stnlist)
        A = find(LISM.Station==stnlist(i));
        asA = [asA; A];
    end;

    CH4M = LISM.CH4_mean_nmolkg(asA);
    N2OM = LISM.N2O_mean_nmolkg(asA); 
    O2M = LISM.O2_umolkg(asA);
    DCH4M = LISM.DCH4(asA);
    DN2OM = LISM.DN2O(asA);
    TM = LISM.T(asA);
    DO2M = LISM.DO2(asA);

    [median(CH4M) prctile(CH4M,25) prctile(CH4M,75)]
    [mean(CH4M) std(CH4M)]

    [median(DCH4M) prctile(DCH4M,25) prctile(DCH4M,75)]
    [mean(DCH4M) std(DCH4M)]

    [median(N2OM) prctile(N2OM,25) prctile(N2OM,75)]
    [mean(N2OM) std(N2OM)]

    [median(DN2OM) prctile(DN2OM,25) prctile(DN2OM,75)]
    [mean(DN2OM) std(DN2OM)]
%%
    [median(O2M) prctile(O2M,25) prctile(O2M,75)]
    [mean(O2M) std(O2M)]

    [median(DO2M) prctile(DO2M,25) prctile(DO2M,75)]
    [mean(DO2M) std(DO2M)]


%%
% we need to make a grid that is evenly spaced so that all the casts are
% interpolated onto the same spacing

dl = [0:0.1:34]'; % depth spacing for interpolation
%x = km_cumulativeA'; % temporary x values for interpolation
x = [-0.5 km_cumulativeA' 19]; % now set values to have a bit of extra start and end
n_stn = numel(x);

dl_grid = repmat(dl,1,n_stn);

% modify the station list to add a bit of distance extending past the first
% station
stnlist = ["EXR1-cast01"
    "EXRX-cast02"
    "MID3-cast01"
    "MID4-cast14"
    "MID5-cast01"
    "WLIS-cast02"
    "WLI6-cast01"]; % this matches cast 16-22 and x as listed above

% AUGUST 2023
km_gridA = repmat(x,length(dl),1);
CH4iA = nan(numel(dl),n_stn); % make a blank grid for storing CH4
N2OiA = nan(numel(dl),n_stn); % make a blank grid for storing CH4
O2iA = nan(numel(dl),n_stn); % make a blank grid for storing CH4

CH4iA = nan(numel(dl),n_stn); % make a blank grid for storing CH4
N2OiA = nan(numel(dl),n_stn); % make a blank grid for storing N2O
O2iA = nan(numel(dl),n_stn); % make a blank grid for storing O2
SiA = nan(numel(dl),n_stn); % make a blank grid for storing S
PDeniA = nan(numel(dl),n_stn); % make a blank grid for storing PDen
tiA = repmat(datetime(0,0,0), 1, n_stn);
DCH4iA = nan(numel(dl),n_stn); % make a blank grid for storing CH4
DN2OiA = nan(numel(dl),n_stn); % make a blank grid for storing N2O
DO2iA = nan(numel(dl),n_stn); % make a blank grid for storing O2
DCH4_nmolkgiA = nan(numel(dl),n_stn); % make a blank grid for storing CH4
DN2O_nmolkgiA = nan(numel(dl),n_stn); % make a blank grid for storing N2O
DO2_umolkgiA = nan(numel(dl),n_stn); % make a blank grid for storing O2

asA = []; % indices containing samples we want to use

LISA.km_cumulative = nan.*LISA.T;

dl_gridA = repmat(dl,1,n_stn);
x_gridA = repmat(x,length(dl),1);

CH4iA = nan(numel(dl),n_stn); % make a blank grid for storing CH4
N2OiA = nan(numel(dl),n_stn); % make a blank grid for storing N2O
O2iA = nan(numel(dl),n_stn); % make a blank grid for storing O2
SiA = nan(numel(dl),n_stn); % make a blank grid for storing S
PDeniA = nan(numel(dl),n_stn); % make a blank grid for storing PDen
tiA = repmat(datetime(0,0,0), 1, n_stn);

CH4iavgA = nan(1,n_stn);
N2OiavgA = nan(1,n_stn);
O2iavgA = nan(1,n_stn);

CH4avgA = nan(1,n_stn);
N2OavgA = nan(1,n_stn);
O2avgA = nan(1,n_stn);


%for i = 1:numel(stnlist)
for i = 1:numel(stnlist)    
    q = find(LISA.Station==stnlist(i));
    asA = [asA; q];
    LISA.km_cumulative(q) = x(i+1); % make sure we only do this for the real stations and not fillers

    tiA(i+1) = LISA.datetime(q(1));
    CH4iA(:,i+1) = interp1(LISA.Depth(q),LISA.mean_CH4_nM(q),dl);
    N2OiA(:,i+1) = interp1(LISA.Depth(q),LISA.mean_N2O_nM(q),dl);
    SiA(:,i+1) = interp1(LISA.Depth(q),LISA.S(q),dl);
    O2iA(:,i+1) = interp1(LISA.Depth(q),LISA.O2_umolkg(q),dl);
    PDeniA(:,i+1) = interp1(LISA.Depth(q),LISA.PDen(q),dl); 
    DCH4iA(:,i+1) = interp1(LISA.Depth(q),LISA.DCH4(q),dl);
    DN2OiA(:,i+1) = interp1(LISA.Depth(q),LISA.DN2O(q),dl);
    DO2iA(:,i+1) = interp1(LISA.Depth(q),LISA.DO2(q),dl);    
    DCH4_nmolkgiA(:,i+1) = interp1(LISA.Depth(q),LISA.DCH4_nmolkg(q),dl);
    DN2O_nmolkgiA(:,i+1) = interp1(LISA.Depth(q),LISA.DN2O_nmolkg(q),dl);
    DO2_umolkgiA(:,i+1) = interp1(LISA.Depth(q),LISA.DO2_umolkg(q),dl); 

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

    nnC = find(~isnan(DCH4_nmolkgiA(:,i+1)));
    DCH4_nmolkgiA(1:min(nnC),i+1) = DCH4_nmolkgiA(min(nnC),i+1);
    DCH4_nmolkgiA(max(nnC):end,i+1) = DCH4_nmolkgiA(max(nnC),i+1);

    nnC = find(~isnan(DN2O_nmolkgiA(:,i+1)));
    DN2O_nmolkgiA(1:min(nnC),i+1) = DN2O_nmolkgiA(min(nnC),i+1);
    DN2O_nmolkgiA(max(nnC):end,i+1) = DN2O_nmolkgiA(max(nnC),i+1);

    nnC = find(~isnan(DO2_umolkgiA(:,i+1)));
    DO2_umolkgiA(1:min(nnC),i+1) = DO2_umolkgiA(min(nnC),i+1);
    DO2_umolkgiA(max(nnC):end,i+1) = DO2_umolkgiA(max(nnC),i+1);      
end;

dmin = 1.6; % minimum depth to plot 
CH4iA(dl_gridA<dmin) = NaN;
N2OiA(dl_gridA<dmin) = NaN;
SiA(dl_gridA<dmin) = NaN;
O2iA(dl_gridA<dmin) = NaN;
PDeniA(dl_gridA<dmin) = NaN;
DCH4iA(dl_gridA<dmin) = NaN;
DN2OiA(dl_gridA<dmin) = NaN;
DO2iA(dl_gridA<dmin) = NaN;
DCH4_nmolkgiA(dl_gridA<dmin) = NaN;
DN2O_nmolkgiA(dl_gridA<dmin) = NaN;
DO2_umolkgiA(dl_gridA<dmin) = NaN;


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

DCH4_nmolkgiA(:,1) = DCH4_nmolkgiA(:,2);
DCH4_nmolkgiA(:,end) = DCH4_nmolkgiA(:,end-1);

DN2O_nmolkgiA(:,1) = DN2O_nmolkgiA(:,2);
DN2O_nmolkgiA(:,end) = DN2O_nmolkgiA(:,end-1);

DO2_umolkgiA(:,1) = DO2_umolkgiA(:,2);
DO2_umolkgiA(:,end) = DO2_umolkgiA(:,end-1);


%DmaxA(1) = DmaxA(2);

% now apply extrapolation at surface and bottom
%CH4iA(~isnan(CH4iA));


%%

% OCTOBER 2023
km_gridO = repmat(x,length(dl),1);
CH4iO = nan(numel(dl),n_stn); % make a blank grid for storing CH4
N2OiO = nan(numel(dl),n_stn); % make a blank grid for storing CH4
O2iO = nan(numel(dl),n_stn); % make a blank grid for storing CH4

CH4iO = nan(numel(dl),n_stn); % make a blank grid for storing CH4
N2OiO = nan(numel(dl),n_stn); % make a blank grid for storing N2O
O2iO = nan(numel(dl),n_stn); % make a blank grid for storing O2
SiO = nan(numel(dl),n_stn); % make a blank grid for storing S
PDeniO = nan(numel(dl),n_stn); % make a blank grid for storing PDen
tiO = repmat(datetime(0,0,0), 1, n_stn);
DCH4iO = nan(numel(dl),n_stn); % make a blank grid for storing CH4
DN2OiO = nan(numel(dl),n_stn); % make a blank grid for storing N2O
DO2iO = nan(numel(dl),n_stn); % make a blank grid for storing O2
DCH4_nmolkgiO = nan(numel(dl),n_stn); % make a blank grid for storing CH4
DN2O_nmolkgiO = nan(numel(dl),n_stn); % make a blank grid for storing N2O
DO2_umolkgiO = nan(numel(dl),n_stn); % make a blank grid for storing O2

asO = []; % indices containing samples we want to use

LISO.km_cumulative = nan.*LISO.T;

dl_gridO = repmat(dl,1,n_stn);
x_gridO = repmat(x,length(dl),1);

CH4iO = nan(numel(dl),n_stn); % make a blank grid for storing CH4
N2OiO = nan(numel(dl),n_stn); % make a blank grid for storing N2O
O2iO = nan(numel(dl),n_stn); % make a blank grid for storing O2
SiO = nan(numel(dl),n_stn); % make a blank grid for storing S
PDeniO = nan(numel(dl),n_stn); % make a blank grid for storing PDen
tiO = repmat(datetime(0,0,0), 1, n_stn);

for i = 1:numel(stnlist)
    q = find(LISO.Station==stnlist(i));
    asO = [asO; q];
    LISO.km_cumulative(q) = x(i+1); % make sure we only do this for the real stations and not fillers

    tiO(i+1) = LISO.datetime(q(1));
    CH4iO(:,i+1) = interp1(LISO.Depth(q),LISO.mean_CH4_nM(q),dl);
    N2OiO(:,i+1) = interp1(LISO.Depth(q),LISO.mean_N2O_nM(q),dl);
    SiO(:,i+1) = interp1(LISO.Depth(q),LISO.S(q),dl);
    O2iO(:,i+1) = interp1(LISO.Depth(q),LISO.O2_umolkg(q),dl);
    PDeniO(:,i+1) = interp1(LISO.Depth(q),LISO.PDen(q),dl); 
    DCH4iO(:,i+1) = interp1(LISO.Depth(q),LISO.DCH4(q),dl);
    DN2OiO(:,i+1) = interp1(LISO.Depth(q),LISO.DN2O(q),dl);
    DO2iO(:,i+1) = interp1(LISO.Depth(q),LISO.DO2(q),dl);        
    DCH4_nmolkgiO(:,i+1) = interp1(LISO.Depth(q),LISO.DCH4_nmolkg(q),dl);
    DN2O_nmolkgiO(:,i+1) = interp1(LISO.Depth(q),LISO.DN2O_nmolkg(q),dl);
    DO2_umolkgiO(:,i+1) = interp1(LISO.Depth(q),LISO.DO2_umolkg(q),dl); 


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

    nnC = find(~isnan(DCH4_nmolkgiO(:,i+1)));
    DCH4_nmolkgiO(1:min(nnC),i+1) = DCH4_nmolkgiO(min(nnC),i+1);
    DCH4_nmolkgiO(max(nnC):end,i+1) = DCH4_nmolkgiO(max(nnC),i+1);

    nnC = find(~isnan(DN2O_nmolkgiO(:,i+1)));
    DN2O_nmolkgiO(1:min(nnC),i+1) = DN2O_nmolkgiO(min(nnC),i+1);
    DN2O_nmolkgiO(max(nnC):end,i+1) = DN2O_nmolkgiO(max(nnC),i+1);

    nnC = find(~isnan(DO2_umolkgiO(:,i+1)));
    DO2_umolkgiO(1:min(nnC),i+1) = DO2_umolkgiO(min(nnC),i+1);
    DO2_umolkgiO(max(nnC):end,i+1) = DO2_umolkgiO(max(nnC),i+1);      
end;

dmin = 1.3; % minimum depth to plot 
CH4iO(dl_gridA<dmin) = NaN;
N2OiO(dl_gridA<dmin) = NaN;
SiO(dl_gridA<dmin) = NaN;
O2iO(dl_gridA<dmin) = NaN;
PDeniO(dl_gridA<dmin) = NaN;
DCH4iO(dl_gridO<dmin) = NaN;
DN2OiO(dl_gridO<dmin) = NaN;
DO2iO(dl_gridO<dmin) = NaN;
DCH4_nmolkgiO(dl_gridO<dmin) = NaN;
DN2O_nmolkgiO(dl_gridO<dmin) = NaN;
DO2_umolkgiO(dl_gridO<dmin) = NaN;

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

DCH4_nmolkgiO(:,1) = DCH4_nmolkgiO(:,2);
DCH4_nmolkgiO(:,end) = DCH4_nmolkgiO(:,end-1);

DN2O_nmolkgiO(:,1) = DN2O_nmolkgiO(:,2);
DN2O_nmolkgiO(:,end) = DN2O_nmolkgiO(:,end-1);

DO2_umolkgiO(:,1) = DO2_umolkgiO(:,2);
DO2_umolkgiO(:,end) = DO2_umolkgiO(:,end-1);

%% MAY 2024 

% MAY 2024
km_gridM = repmat(x,length(dl),1);
CH4iM = nan(numel(dl),n_stn); % make a blank grid for storing CH4
N2OiM = nan(numel(dl),n_stn); % make a blank grid for storing CH4
O2iM = nan(numel(dl),n_stn); % make a blank grid for storing CH4

CH4iM = nan(numel(dl),n_stn); % make a blank grid for storing CH4
N2OiM = nan(numel(dl),n_stn); % make a blank grid for storing N2O
O2iM = nan(numel(dl),n_stn); % make a blank grid for storing O2
SiM = nan(numel(dl),n_stn); % make a blank grid for storing S
PDeniM = nan(numel(dl),n_stn); % make a blank grid for storing PDen
tiM = repmat(datetime(0,0,0), 1, n_stn);

DCH4iM = nan(numel(dl),n_stn); % make a blank grid for storing CH4
DN2OiM = nan(numel(dl),n_stn); % make a blank grid for storing N2O
DO2iM = nan(numel(dl),n_stn); % make a blank grid for storing O
DCH4_nmolkgiM = nan(numel(dl),n_stn); % make a blank grid for storing CH4
DN2O_nmolkgiM = nan(numel(dl),n_stn); % make a blank grid for storing N2O
DO2_umolkgiM = nan(numel(dl),n_stn); % make a blank grid for storing O2

asM = []; % indices containing samples we want to use

LISM.km_cumulative = nan.*LISM.T;

dl_gridM = repmat(dl,1,n_stn);
x_gridM = repmat(x,length(dl),1);

CH4iM = nan(numel(dl),n_stn); % make a blank grid for storing CH4
N2OiM = nan(numel(dl),n_stn); % make a blank grid for storing CH4
O2iM = nan(numel(dl),n_stn); % make a blank grid for storing CH4

CH4iM = nan(numel(dl),n_stn); % make a blank grid for storing CH4
N2Oia = nan(numel(dl),n_stn); % make a blank grid for storing N2O
O2iM = nan(numel(dl),n_stn); % make a blank grid for storing O2
SiM = nan(numel(dl),n_stn); % make a blank grid for storing S
PDeniM = nan(numel(dl),n_stn); % make a blank grid for storing PDen
tiM = repmat(datetime(0,0,0), 1, n_stn);

for i = 1:numel(stnlist)
    q = find(LISM.Station==stnlist(i));
    asM = [asM; q];
    LISM.km_cumulative(q) = x(i+1); % make sure we only do this for the real stations and not fillers

    tiM(i+1) = LISM.datetime(q(1));
    CH4iM(:,i+1) = interp1(LISM.Depth(q),LISM.mean_CH4_nM(q),dl);
    N2OiM(:,i+1) = interp1(LISM.Depth(q),LISM.mean_N2O_nM(q),dl);
    SiM(:,i+1) = interp1(LISM.Depth(q),LISM.S(q),dl);
    O2iM(:,i+1) = interp1(LISM.Depth(q),LISM.O2_umolkg(q),dl);
    PDeniM(:,i+1) = interp1(LISM.Depth(q),LISM.PDen(q),dl); 
    DCH4iM(:,i+1) = interp1(LISM.Depth(q),LISM.DCH4(q),dl);
    DN2OiM(:,i+1) = interp1(LISM.Depth(q),LISM.DN2O(q),dl);
    DO2iM(:,i+1) = interp1(LISM.Depth(q),LISM.DO2(q),dl);        
    DCH4_nmolkgiM(:,i+1) = interp1(LISM.Depth(q),LISM.DCH4_nmolkg(q),dl);
    DN2O_nmolkgiM(:,i+1) = interp1(LISM.Depth(q),LISM.DN2O_nmolkg(q),dl);
    DO2_umolkgiM(:,i+1) = interp1(LISM.Depth(q),LISM.DO2_umolkg(q),dl);

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

    nnC = find(~isnan(DCH4_nmolkgiM(:,i+1)));
    DCH4_nmolkgiM(1:min(nnC),i+1) = DCH4_nmolkgiM(min(nnC),i+1);
    DCH4_nmolkgiM(max(nnC):end,i+1) = DCH4_nmolkgiM(max(nnC),i+1);

    nnC = find(~isnan(DN2O_nmolkgiM(:,i+1)));
    DN2O_nmolkgiM(1:min(nnC),i+1) = DN2O_nmolkgiM(min(nnC),i+1);
    DN2O_nmolkgiM(max(nnC):end,i+1) = DN2O_nmolkgiM(max(nnC),i+1);

    nnC = find(~isnan(DO2_umolkgiM(:,i+1)));
    DO2_umolkgiM(1:min(nnC),i+1) = DO2_umolkgiM(min(nnC),i+1);
    DO2_umolkgiM(max(nnC):end,i+1) = DO2_umolkgiM(max(nnC),i+1);       
end;

dmin = 1.3; % minimum depth to plot 
CH4iM(dl_gridA<dmin) = NaN;
N2OiM(dl_gridA<dmin) = NaN;
SiM(dl_gridA<dmin) = NaN;
O2iM(dl_gridA<dmin) = NaN;
PDeniM(dl_gridA<dmin) = NaN;
DCH4iM(dl_gridM<dmin) = NaN;
DN2OiM(dl_gridM<dmin) = NaN;
DO2iM(dl_gridM<dmin) = NaN;
DCH4_nmolkgiM(dl_gridM<dmin) = NaN;
DN2O_nmolkgiM(dl_gridM<dmin) = NaN;
DO2_umolkgiM(dl_gridM<dmin) = NaN;

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

DCH4_nmolkgiM(:,1) = DCH4_nmolkgiM(:,2);
DCH4_nmolkgiM(:,end) = DCH4_nmolkgiM(:,end-1);

DN2O_nmolkgiM(:,1) = DN2O_nmolkgiM(:,2);
DN2O_nmolkgiM(:,end) = DN2O_nmolkgiM(:,end-1);

DO2_umolkgiM(:,1) = DO2_umolkgiM(:,2);
DO2_umolkgiM(:,end) = DO2_umolkgiM(:,end-1);

%% PLOT THE DATA
% help for getting ranges
% Density
[min(LISO.PDen) max(LISO.PDen)]

[min(LISM.mean_N2O_nM) max(LISM.mean_N2O_nM)]
%%

[min(min(CH4iA)) max(max(CH4iA))
min(min(N2OiA)) max(max(N2OiA))
min(min(O2iA)) max(max(O2iA))]

% August
% 44.5834  445.6192
% 8.5353   13.7525
% 39.5337  211.9365


% October
% 25.7143  187.9974
% 9.4866   14.1598
% 197.4300  246.3629

% May
% 26.3480   87.7741
% 10.7882   16.5057
% 199.6387  323.3287

%%
nr = 4; % number of rows
nc = 3; % number of columns
lw = 2; % default line width
fs = 10; % default font size
ms = 3; % default marker size
yt = [0 10 20 30]; %y-axis ticks;
xt = [0 5 10 15];
yl = [0 32]; % y-axis limit in m
xl = [-0.5 19]; %x-axis limit in km

caCH4 = [25 350]; % CH4 axis limits
clevelCH4 = [caCH4(1):1:caCH4(2)]; %CH4 colorbar levels

caN2O = [8.5 16.6]; % N2O axis limits
clevelN2O = [caN2O(1):0.1:caN2O(2)]; % N2O colorbar levels

caO2 = [35 325];
clevelO2 = [caO2(1):1:caO2(2)];

caPDen = [16 18.8]; % PDen axis limits
clevelPDen = [caPDen(1):0.02:caPDen(2)]; % PDen colorbar levels

fig=figure(23);
clf;
sp=tight_subplot(nr,nc,[.025 .025],[.08 .04],[.08 .04]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 12 8]);
set(gcf,'renderer','painters');
    set(gcf,'GraphicsSmoothing','on')

% rows: Aug / Oct / May
% columns: CH4 / N2O / O2 / Density
% SUBPLOT 1
subplot(sp(1))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridA,dl_grid,CH4iA,clevelCH4,'edgecolor','none');
    plot(LISA.km_cumulative(asA),LISA.Depth(asA),'+k', 'linewidth', 2, 'markersize',ms)

    caxis(caCH4)
    %c = colorbar('location','eastoutside');
    %c.Label.String = 'CH_4 (nmol/kg)';
    %c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xl);
    set(gca,'ytick',yt);
    set(gca,'yticklabel',yt);
    set(gca,'xtick',xt);
    %set(gca,'xticklabel',xt);
    %xlabel('Transect distance [km]');
    ylabel('Depth [m]');
    title('August');

    %add in plot of bathymetry
    dmA = [yl(2) DmaxA(1) DmaxA' DmaxA(end) yl(2)];
    kcA = [x_gridA(1,1) x_gridA(1,:) x_gridA(1,end)];
    fill(kcA,dmA,[0.5 0.5 0.5],'edgecolor','none');
    axis ij;

% SUBPLOT 2
subplot(sp(2))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridO,dl_grid,CH4iO,clevelCH4,'edgecolor','none');
    plot(LISO.km_cumulative(asO),LISO.Depth(asO),'+k', 'linewidth', 2, 'markersize',ms)

    %xlabel('Transect distance [km]');
    %ylabel('Depth [m]');
    caxis(caCH4)
    %c = colorbar('location','eastoutside');
    %c.Label.String = 'CH_4 (nmol/kg)';
    %c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xl);
    set(gca,'ytick',yt);
    %set(gca,'yticklabel',yt);
    set(gca,'xtick',xt);
    %set(gca,'xticklabel',xt);
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
    plot(LISM.km_cumulative(asM),LISM.Depth(asM),'+k', 'linewidth', 2, 'markersize',ms)

    %xlabel('Transect distance [km]');
    %ylabel('Depth [m]');
    caxis(caCH4)
    c = colorbar('location','eastoutside');
    c.Label.String = 'CH_4 (nmol/kg)';
    c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xl);
    set(gca,'ytick',yt);
    %set(gca,'yticklabel',yt);
    set(gca,'xtick',xt);
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
    plot(LISA.km_cumulative(asA),LISA.Depth(asA),'+k', 'linewidth', 2, 'markersize',ms)

    caxis(caN2O)
    %c = colorbar('location','eastoutside');
    %c.Label.String = 'N_2O (nmol/kg)';
    %c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xl);
    set(gca,'ytick',yt);
    set(gca,'yticklabel',yt);
    set(gca,'xtick',xt);
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
    plot(LISO.km_cumulative(asO),LISO.Depth(asO),'+k', 'linewidth', 2, 'markersize',ms)

    %xlabel('Transect distance [km]');
    %ylabel('Depth [m]');
    caxis(caN2O)
    %c = colorbar('location','eastoutside');
    %c.Label.String = 'CH_4 (nmol/kg)';
    %c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xl);
    set(gca,'ytick',yt);
    %set(gca,'yticklabel',yt);
    set(gca,'xtick',xt);
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
    plot(LISM.km_cumulative(asM),LISM.Depth(asM),'+k', 'linewidth', 2, 'markersize',ms)

    %xlabel('Transect distance [km]');
    %ylabel('Depth [m]');
    caxis(caN2O)
    c = colorbar('location','eastoutside');
    c.Label.String = 'N_2O (nmol/kg)';
    c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xl);
    set(gca,'ytick',yt);
    %set(gca,'yticklabel',yt);
    set(gca,'xtick',xt);
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
    plot(LISA.km_cumulative(asA),LISA.Depth(asA),'+k', 'linewidth', 2, 'markersize',ms)

    caxis(caO2)
    %c = colorbar('location','eastoutside');
    %c.Label.String = 'N_2O (nmol/kg)';
    %c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xl);
    set(gca,'ytick',yt);
    set(gca,'yticklabel',yt);
    set(gca,'xtick',xt);
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
    plot(LISO.km_cumulative(asO),LISO.Depth(asO),'+k', 'linewidth', 2, 'markersize',ms)

    %xlabel('Transect distance [km]');
    %ylabel('Depth [m]');
    caxis(caO2)
    %c = colorbar('location','eastoutside');
    %c.Label.String = 'CH_4 (nmol/kg)';
    %c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xl);
    set(gca,'ytick',yt);
    %set(gca,'yticklabel',yt);
    set(gca,'xtick',xt);
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
    plot(LISM.km_cumulative(asM),LISM.Depth(asM),'+k', 'linewidth', 2, 'markersize',ms)

    %xlabel('Transect distance [km]');
    %ylabel('Depth [m]');
    caxis(caO2)
    c = colorbar('location','eastoutside');
    c.Label.String = 'O_2 (\mumol/kg)';
    c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xl);
    set(gca,'ytick',yt);
    %set(gca,'yticklabel',yt);
    set(gca,'xtick',xt);
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
    plot(LISA.km_cumulative(asA),LISA.Depth(asA),'+k', 'linewidth', 2, 'markersize',ms)

    caxis(caPDen)
    %c = colorbar('location','eastoutside');
    %c.Label.String = 'N_2O (nmol/kg)';
    %c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xl);
    set(gca,'ytick',yt);
    set(gca,'yticklabel',yt);
    set(gca,'xtick',xt);
    set(gca,'xticklabel',xt);
    xlabel('Transect distance [km]');
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
    plot(LISO.km_cumulative(asO),LISO.Depth(asO),'+k', 'linewidth', 2, 'markersize',ms)

    %xlabel('Transect distance [km]');
    %ylabel('Depth [m]');
    clim(caPDen)
    %c = colorbar('location','eastoutside');
    %c.Label.String = 'CH_4 (nmol/kg)';
    %c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xl);
    set(gca,'ytick',yt);
    %set(gca,'yticklabel',yt);
    set(gca,'xtick',xt);
    set(gca,'xticklabel',xt);
    xlabel('Transect distance [km]');

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
    plot(LISM.km_cumulative(asM),LISM.Depth(asM),'+k', 'linewidth', 2, 'markersize',ms)

    %xlabel('Transect distance [km]');
    %ylabel('Depth [m]');
    clim(caPDen)
    c = colorbar('location','eastoutside');
    c.Label.String = '\sigma_{\theta} (kg/m^3)';
    c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xl);
    set(gca,'ytick',yt);
    %set(gca,'yticklabel',yt);
    set(gca,'xtick',xt);
    set(gca,'xticklabel',xt);
    xlabel('Transect distance [km]');

    %add in plot of bathymetry
    dmM = [yl(2) DmaxM(1) DmaxM' DmaxM(end) yl(2)];
    kcM = [x_gridM(1,1) x_gridM(1,:) x_gridM(1,end)];
    fill(kcM,dmM,[0.5 0.5 0.5],'edgecolor','none');
    axis ij;
    
    wysiwyg;

%print(gcf, '-dpng', '-r300', 'LIS_Transect_CH4_N2O_O2_PDen.png');
%print(gcf,'-depsc','-vector','LIS_Transect_CH4_N2O_O2_PDen.eps');
%epsclean('LIS_Transect_CH4_N2O_O2_PDen.eps','LIS_Transect_CH4_N2O_O2_PDen.eps');

%%

[min(min(DCH4_nmolkgiA)) max(max(DCH4_nmolkgiA))
min(min(DN2O_nmolkgiA)) max(max(DN2O_nmolkgiA))
min(min(DO2_umolkgiA)) max(max(DO2_umolkgiA))]

[min(min(DCH4_nmolkgiO)) max(max(DCH4_nmolkgiO))
min(min(DN2O_nmolkgiO)) max(max(DN2O_nmolkgiO))
min(min(DO2_umolkgiO)) max(max(DO2_umolkgiO))]

[min(min(DCH4_nmolkgiM)) max(max(DCH4_nmolkgiM))
min(min(DN2O_nmolkgiM)) max(max(DN2O_nmolkgiM))
min(min(DO2_umolkgiM)) max(max(DO2_umolkgiM))]

[min(min(PDeniA)) max(max(PDeniA))
min(min(PDeniO)) max(max(PDeniO))
min(min(PDeniM)) max(max(PDeniM))]




%%
nr = 4; % number of rows
nc = 3; % number of columns
lw = 2; % default line width
fs = 10; % default font size
ms = 3; % default marker size
yt = [0 10 20 30]; %y-axis ticks;
xt = [0 5 10 15];
yl = [0 32]; % y-axis limit in m
xl = [-0.5 19]; %x-axis limit in km

fc = [0.75 0.75 0.75]; %fill color

fsl = 12; % font size for axis label

caCH4 = [20 435]; % CH4 axis limits
clevelCH4 = [caCH4(1):1:caCH4(2)]; %CH4 colorbar levels

caN2O = [0 7]; % N2O axis limits
clevelN2O = [caN2O(1):0.1:caN2O(2)]; % N2O colorbar levels

caO2 = [-200 70];
clevelO2 = [caO2(1):1:caO2(2)];

caPDen = [16 18.8]; % PDen axis limits
clevelPDen = [caPDen(1):0.02:caPDen(2)]; % PDen colorbar levels

fig=figure(23);
clf;
sp=tight_subplot(nr,nc,[.025 .025],[.08 .04],[.08 .04]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 12 8]);
set(gcf,'renderer','painters');
    set(gcf,'GraphicsSmoothing','on')

% rows: Aug / Oct / May
% columns: CH4 / N2O / O2 / Density
% SUBPLOT 1
subplot(sp(1))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridA,dl_grid,DCH4_nmolkgiA,clevelCH4,'edgecolor','none');
    plot(LISA.km_cumulative(asA),LISA.Depth(asA),'+k', 'linewidth', 2, 'markersize',ms)

    caxis(caCH4)
    %c = colorbar('location','eastoutside');
    %c.Label.String = 'CH_4 (nmol/kg)';
    %c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xl);
    set(gca,'ytick',yt);
    set(gca,'yticklabel',yt);
    set(gca,'xtick',xt);
    %set(gca,'xticklabel',xt);
    %xlabel('Transect distance [km]');
    ylabel('Depth [m]');
    title('August');

    %add in plot of bathymetry
    dmA = [yl(2) DmaxA(1) DmaxA' DmaxA(end) yl(2)];
    kcA = [x_gridA(1,1) x_gridA(1,:) x_gridA(1,end)];
    fill(kcA,dmA,fc,'edgecolor','none');
    text(17,27,'a','color',[0 0 0],'fontsize',fsl);
    axis ij;


% SUBPLOT 2
subplot(sp(2))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridO,dl_grid,DCH4_nmolkgiO,clevelCH4,'edgecolor','none');
    plot(LISO.km_cumulative(asO),LISO.Depth(asO),'+k', 'linewidth', 2, 'markersize',ms)

    %xlabel('Transect distance [km]');
    %ylabel('Depth [m]');
    caxis(caCH4)
    %c = colorbar('location','eastoutside');
    %c.Label.String = 'CH_4 (nmol/kg)';
    %c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xl);
    set(gca,'ytick',yt);
    %set(gca,'yticklabel',yt);
    set(gca,'xtick',xt);
    %set(gca,'xticklabel',xt);
    title('October');

    %add in plot of bathymetry
    dmO = [yl(2) DmaxO(1) DmaxO' DmaxO(end) yl(2)];
    kcO = [x_gridO(1,1) x_gridO(1,:) x_gridO(1,end)];
    fill(kcA,dmA,fc,'edgecolor','none');
    text(17,27,'b','color',[0 0 0],'fontsize',fsl);
    axis ij;


% SUBPLOT 3
subplot(sp(3))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    clevel = [25:1:460];
    ca = [25 460]; %colorbar limits
    C = contourf(x_gridM,dl_grid,DCH4_nmolkgiM,clevelCH4,'edgecolor','none');
    plot(LISM.km_cumulative(asM),LISM.Depth(asM),'+k', 'linewidth', 2, 'markersize',ms)

    %xlabel('Transect distance [km]');
    %ylabel('Depth [m]');
    caxis(caCH4)
    c = colorbar('location','eastoutside');
    c.Label.String = '\DeltaCH_4 (nmol/kg)';
    c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xl);
    set(gca,'ytick',yt);
    %set(gca,'yticklabel',yt);
    set(gca,'xtick',xt);
    %set(gca,'xticklabel',xt);
    title(['May'])

    %add in plot of bathymetry
    dmM = [yl(2) DmaxM(1) DmaxM' DmaxM(end) yl(2)];
    kcM = [x_gridM(1,1) x_gridM(1,:) x_gridM(1,end)];
    fill(kcA,dmA,fc,'edgecolor','none');
    text(17,27,'c','color',[0 0 0],'fontsize',fsl);
    axis ij;

%SUBPLOT 4
subplot(sp(4))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridA,dl_grid,DN2O_nmolkgiA,clevelN2O,'edgecolor','none');
    plot(LISA.km_cumulative(asA),LISA.Depth(asA),'+k', 'linewidth', 2, 'markersize',ms)

    caxis(caN2O)
    %c = colorbar('location','eastoutside');
    %c.Label.String = 'N_2O (nmol/kg)';
    %c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xl);
    set(gca,'ytick',yt);
    set(gca,'yticklabel',yt);
    set(gca,'xtick',xt);
    %set(gca,'xticklabel',xt);
    %xlabel('Transect distance [km]');
    ylabel('Depth [m]');
    %title('August');

    %add in plot of bathymetry
    dmA = [yl(2) DmaxA(1) DmaxA' DmaxA(end) yl(2)];
    kcA = [x_gridA(1,1) x_gridA(1,:) x_gridA(1,end)];
    fill(kcA,dmA,fc,'edgecolor','none');
    text(17,27,'d','color',[0 0 0],'fontsize',fsl);
    axis ij;


% SUBPLOT 5
subplot(sp(5))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridO,dl_grid,DN2O_nmolkgiO,clevelN2O,'edgecolor','none');
    plot(LISO.km_cumulative(asO),LISO.Depth(asO),'+k', 'linewidth', 2, 'markersize',ms)

    %xlabel('Transect distance [km]');
    %ylabel('Depth [m]');
    caxis(caN2O)
    %c = colorbar('location','eastoutside');
    %c.Label.String = 'CH_4 (nmol/kg)';
    %c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xl);
    set(gca,'ytick',yt);
    %set(gca,'yticklabel',yt);
    set(gca,'xtick',xt);
    %set(gca,'xticklabel',xt);

    %add in plot of bathymetry
    dmM = [yl(2) DmaxO(1) DmaxO' DmaxO(end) yl(2)];
    kcO = [x_gridO(1,1) x_gridO(1,:) x_gridO(1,end)];
     fill(kcA,dmA,fc,'edgecolor','none');
    text(17,27,'e','color',[0 0 0],'fontsize',fsl);
    axis ij;

% SUBPLOT 6
subplot(sp(6))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridM,dl_grid,DN2O_nmolkgiM,clevelN2O,'edgecolor','none');
    plot(LISM.km_cumulative(asM),LISM.Depth(asM),'+k', 'linewidth', 2, 'markersize',ms)

    %xlabel('Transect distance [km]');
    %ylabel('Depth [m]');
    caxis(caN2O)
    c = colorbar('location','eastoutside');
    c.Label.String = '\DeltaN_2O (nmol/kg)';
    c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xl);
    set(gca,'ytick',yt);
    %set(gca,'yticklabel',yt);
    set(gca,'xtick',xt);
    %set(gca,'xticklabel',xt);

    %add in plot of bathymetry
    dmM = [yl(2) DmaxM(1) DmaxM' DmaxM(end) yl(2)];
    kcM = [x_gridM(1,1) x_gridM(1,:) x_gridM(1,end)];
    fill(kcA,dmA,fc,'edgecolor','none');
    text(17,27,'f','color',[0 0 0],'fontsize',fsl);
    axis ij;

%SUBPLOT 7
subplot(sp(7))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridA,dl_grid,DO2_umolkgiA,clevelO2,'edgecolor','none');
    contour(x_gridA,dl_grid,DO2_umolkgiA,0,'edgecolor','k');
    plot(LISA.km_cumulative(asA),LISA.Depth(asA),'+k', 'linewidth', 2, 'markersize',ms)

    caxis(caO2)
    %c = colorbar('location','eastoutside');
    %c.Label.String = 'N_2O (nmol/kg)';
    %c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xl);
    set(gca,'ytick',yt);
    set(gca,'yticklabel',yt);
    set(gca,'xtick',xt);
    %set(gca,'xticklabel',xt);
    %xlabel('Transect distance [km]');
    ylabel('Depth [m]');
    %title('August');

    %add in plot of bathymetry
    dmA = [yl(2) DmaxA(1) DmaxA' DmaxA(end) yl(2)];
    kcA = [x_gridA(1,1) x_gridA(1,:) x_gridA(1,end)];
    fill(kcA,dmA,fc,'edgecolor','none');
    text(17,27,'g','color',[0 0 0],'fontsize',fsl);
    axis ij;


% SUBPLOT 8
subplot(sp(8))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridO,dl_grid,DO2_umolkgiO,clevelO2,'edgecolor','none');
    contour(x_gridO,dl_grid,DO2_umolkgiO,0,'edgecolor','k');
    plot(LISO.km_cumulative(asO),LISO.Depth(asO),'+k', 'linewidth', 2, 'markersize',ms)

    %xlabel('Transect distance [km]');
    %ylabel('Depth [m]');
    caxis(caO2)
    %c = colorbar('location','eastoutside');
    %c.Label.String = 'CH_4 (nmol/kg)';
    %c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xl);
    set(gca,'ytick',yt);
    %set(gca,'yticklabel',yt);
    set(gca,'xtick',xt);
    %set(gca,'xticklabel',xt);

    %add in plot of bathymetry
    dmO = [yl(2) DmaxO(1) DmaxO' DmaxO(end) yl(2)];
    kcO = [x_gridO(1,1) x_gridO(1,:) x_gridO(1,end)];
    fill(kcA,dmA,fc,'edgecolor','none');
    text(17,27,'h','color',[0 0 0],'fontsize',fsl);
    axis ij;

% SUBPLOT 9
subplot(sp(9))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridM,dl_grid,DO2_umolkgiM,clevelO2,'edgecolor','none');
        contour(x_gridM,dl_grid,DO2_umolkgiM,[0 0],'-k');
    plot(LISM.km_cumulative(asM),LISM.Depth(asM),'+k', 'linewidth', 2, 'markersize',ms)

    %xlabel('Transect distance [km]');
    %ylabel('Depth [m]');
    caxis(caO2)
    c = colorbar('location','eastoutside');
    c.Label.String = '\DeltaO_2 (\mumol/kg)';
    c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xl);
    set(gca,'ytick',yt);
    %set(gca,'yticklabel',yt);
    set(gca,'xtick',xt);
    %set(gca,'xticklabel',xt);

    %add in plot of bathymetry
    dmM = [yl(2) DmaxM(1) DmaxM' DmaxM(end) yl(2)];
    kcM = [x_gridM(1,1) x_gridM(1,:) x_gridM(1,end)];
    fill(kcA,dmA,fc,'edgecolor','none');
    text(17,27,'i','color',[0 0 0],'fontsize',fsl);
    axis ij;


%SUBPLOT 10
subplot(sp(10))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridA,dl_grid,PDeniA,clevelPDen,'edgecolor','none');
    plot(LISA.km_cumulative(asA),LISA.Depth(asA),'+k', 'linewidth', 2, 'markersize',ms)

    caxis(caPDen)
    %c = colorbar('location','eastoutside');
    %c.Label.String = 'N_2O (nmol/kg)';
    %c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xl);
    set(gca,'ytick',yt);
    set(gca,'yticklabel',yt);
    set(gca,'xtick',xt);
    set(gca,'xticklabel',xt);
    xlabel('Transect distance [km]');
    ylabel('Depth [m]');
    %title('August');

    %add in plot of bathymetry
    dmA = [yl(2) DmaxA(1) DmaxA' DmaxA(end) yl(2)];
    kcA = [x_gridA(1,1) x_gridA(1,:) x_gridA(1,end)];
    fill(kcA,dmA,fc,'edgecolor','none');
    text(17,27,'j','color',[0 0 0],'fontsize',fsl);
    axis ij;


% SUBPLOT 11
subplot(sp(11))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridO,dl_grid,PDeniO,clevelPDen,'edgecolor','none');
    plot(LISO.km_cumulative(asO),LISO.Depth(asO),'+k', 'linewidth', 2, 'markersize',ms)

    %xlabel('Transect distance [km]');
    %ylabel('Depth [m]');
    clim(caPDen)
    %c = colorbar('location','eastoutside');
    %c.Label.String = 'CH_4 (nmol/kg)';
    %c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xl);
    set(gca,'ytick',yt);
    %set(gca,'yticklabel',yt);
    set(gca,'xtick',xt);
    set(gca,'xticklabel',xt);
    xlabel('Transect distance [km]');

    %add in plot of bathymetry
    dmM = [yl(2) DmaxO(1) DmaxO' DmaxO(end) yl(2)];
    kcO = [x_gridO(1,1) x_gridO(1,:) x_gridO(1,end)];
    fill(kcA,dmA,fc,'edgecolor','none');
    text(17,27,'k','color',[0 0 0],'fontsize',fsl);
    axis ij;

% SUBPLOT 12
subplot(sp(12))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridM,dl_grid,PDeniM,clevelPDen,'edgecolor','none');
  %  C = surf(x_grid,dl_grid,PDeniM,clevelPDen,'facecolor','interp','edgecolor','interp');
  %  surf instead, with 'FaceColor','interp', 'EdgeColor','interp' and %view(0,90)’.  
    plot(LISM.km_cumulative(asM),LISM.Depth(asM),'+k', 'linewidth', 2, 'markersize',ms)

    %xlabel('Transect distance [km]');
    %ylabel('Depth [m]');
    clim(caPDen)
    c = colorbar('location','eastoutside');
    c.Label.String = '\sigma_{\theta} (kg/m^3)';
    c.FontSize = fs;
    set(gca,'layer','top');
    ylim(yl);
    xlim(xl);
    set(gca,'ytick',yt);
    %set(gca,'yticklabel',yt);
    set(gca,'xtick',xt);
    set(gca,'xticklabel',xt);
    xlabel('Transect distance [km]');

    %add in plot of bathymetry
    dmM = [yl(2) DmaxM(1) DmaxM' DmaxM(end) yl(2)];
    kcM = [x_gridM(1,1) x_gridM(1,:) x_gridM(1,end)];
    fill(kcA,dmA,fc,'edgecolor','none');
    text(17,27,'l','color',[0 0 0],'fontsize',fsl);
    axis ij;
    
    wysiwyg;

print(gcf, '-dpng', '-r300', 'LIS_Transect_DCH4_DN2O_DO2_PDen.png');
print(gcf,'-depsc','-vector','LIS_Transect_DCH4_DN2O_DO2_PDen.eps');
epsclean('LIS_Transect_DCH4_DN2O_DO2_PDen.eps','LIS_Transect_DCH4_DN2O_DO2_PDen.eps');

%%

figure(1)
clf; hold on;

subplot(1,3,1)
hold on; box on;
stn = 'EXCR-cast01';
A = find(LIS.Station==stn);
errorbar(LIS.mean_CH4_nM(A),LIS.Depth(A),LIS.std_CH4_nM(A),'horizontal','o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXCR-cast02';
A = find(LIS.Station==stn);
errorbar(LIS.mean_CH4_nM(A),LIS.Depth(A),LIS.std_CH4_nM(A),'horizontal','s-r','linewidth',1.5, 'markerfacecolor','r');
legend('EXCR-cast01','EXCR-cast02','location','northeast')

xlabel('CH_4 (nM)');
ylabel('Depth (m)');
axis ij;

subplot(1,3,2)
hold on; box on;
stn = 'EXCR-cast01';
A = find(LIS.Station==stn);
plot(LIS.T(A),LIS.Depth(A),'o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXCR-cast02';
A = find(LIS.Station==stn);
plot(LIS.T(A),LIS.Depth(A),'s-r','linewidth',1.5, 'markerfacecolor','r');
legend('EXCR-cast01','EXCR-cast02','location','northeast')

xlabel('Temp (^oC)');
ylabel('Depth (m)');
axis ij;

subplot(1,3,3)
hold on; box on;
stn = 'EXCR-cast01';
A = find(LIS.Station==stn);
plot(LIS.O2_umolkg(A),LIS.Depth(A),'o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXCR-cast02';
A = find(LIS.Station==stn);
plot(LIS.O2_umolkg(A),LIS.Depth(A),'s-r','linewidth',1.5, 'markerfacecolor','r');
legend('EXCR-cast01','EXCR-cast02','location','northeast')

xlabel('O_2 (\mumol/kg)');
ylabel('Depth (m)');
axis ij;


%%
% we need to make a grid that is evenly spaced so that all the casts are
% interpolated onto the same spacing

dl = [0:0.1:18]'; % depth spacing for interpolation
x = km_cumulative'; % temporary x values for interpolation
n_stn = numel(x);

dl_grid = repmat(dl,1,n_stn);

x_grid = repmat(x,length(dl),1);

CH4i = nan(numel(dl),n_stn); % make a blank grid for storing CH4
N2Oi = nan(numel(dl),n_stn); % make a blank grid for storing CH4
O2i = nan(numel(dl),n_stn); % make a blank grid for storing CH4

CH4i = nan(numel(dl),n_stn); % make a blank grid for storing CH4
N2Oi = nan(numel(dl),n_stn); % make a blank grid for storing N2O
O2i = nan(numel(dl),n_stn); % make a blank grid for storing O2
Si = nan(numel(dl),n_stn); % make a blank grid for storing S
PDeni = nan(numel(dl),n_stn); % make a blank grid for storing PDen
ti = repmat(datetime(0,0,0), 1, n_stn);

stnlist = ["EXR1-cast01"
    "EXCR-cast02"
    "MID3-cast01"
    "MID4-cast14"
    "MID5-cast01"
    "WLIS-cast02"
    "WLI6-cast01"]; % this matches cast 16-22 and x as listed above

as = [];

LIS.km_cumulative = nan.*LIS.T;

for i = 1:numel(stnlist)
    A = find(LIS.Station==stnlist(i));
    LIS.km_cumulative(A) = x(i);
    as = [as; A];
    ti(i) = LIS.datetime(A(1));
    CH4i(:,i) = interp1(LIS.Depth(A),LIS.mean_CH4_nM(A),dl);
    N2Oi(:,i) = interp1(LIS.Depth(A),LIS.mean_N2O_nM(A),dl);
    Si(:,i) = interp1(LIS.Depth(A),LIS.S(A),dl);
    O2i(:,i) = interp1(LIS.Depth(A),LIS.O2_umolkg(A),dl);
    PDeni(:,i) = interp1(LIS.Depth(A),LIS.PDen(A),dl);    
end;


fs=22;
ms=4;
figure(2)
clf; hold on;
box on;
set(gca,'tickdir','out');
set(gca,'fontsize',fs);
set(gcf, 'PaperUnits', 'inches');
set(gcf,'renderer','painters');
set(gcf, 'PaperPosition', [0 0 6 8]);
set(gca,'linewidth',0.75);
    set(gcf,'GraphicsSmoothing','on')

clevel = [25:1:460];
ca = [25 460]; %colorbar limits
C = contourf(x_grid,dl_grid,CH4i,clevel,'edgecolor','none');
plot(LIS.km_cumulative(as),LIS.Depth(as),'+k', 'linewidth', 2, 'markersize',ms)

xlabel('Transect distance [km]');
ylabel('Depth [m]');
caxis(ca)
c = colorbar('location','southoutside');
c.Label.String = 'CH_4 (nmol/kg)';
c.FontSize = fs;
set(gca,'layer','top');
set(gca,'ytick',[0 5 10 15]);
ylim([0 18]);
axis ij;

print -dpng -r300 LIS_CH4_transect_Oct.png;
%print -dpng -r300 LIS_CH4_transect_Aug.png;

wysiwyg;
%%
figure(3)
clf; hold on;
set(gca,'tickdir','out');
set(gca,'fontsize',fs);
set(gcf, 'PaperUnits', 'inches');
set(gcf,'renderer','painters');
set(gcf, 'PaperPosition', [0 0 6 8]);
set(gca,'linewidth',0.75);
set(gcf,'GraphicsSmoothing','on');
box on;

clevel = [8.5:0.02:14];
ca = [8.5 14];
C = contourf(x_grid,dl_grid,N2Oi,clevel,'edgecolor','none');
plot(LIS.km_cumulative(as),LIS.Depth(as),'+k', 'linewidth', 2, 'markersize',ms)

xlabel('Transect distance [km]');
ylabel('Depth [m]');
caxis(ca)
c = colorbar('location','southoutside');
c.Label.String = 'N_2O (nmol/kg)';
c.FontSize = fs;
set(gca,'layer','top');
set(gca,'ytick',[0 5 10 15]);
ylim([0 18]);
axis ij;
print -dpng -r300 LIS_N2O_transect_Oct.png;
%print -dpng -r300 LIS_N2O_transect_Aug.png;

wysiwyg;
%%

figure(4)
clf; hold on;
set(gca,'tickdir','out');
set(gca,'fontsize',fs);
set(gcf, 'PaperUnits', 'inches');
set(gcf,'renderer','painters');
set(gcf, 'PaperPosition', [0 0 6 8]);
set(gca,'linewidth',0.75);
set(gcf,'GraphicsSmoothing','on');
box on;

clevel = [40:1:250];
ca = [40 250];
C = contourf(x_grid,dl_grid,O2i,clevel,'edgecolor','none');

stn = 'EXR1-cast01';
A = find(LIS.Station==stn);
plot(x(1).*ones(numel(A)),LIS.Depth(A), '+k', 'LineWidth', 2, 'MarkerSize', 4);

stn = 'EXCR-cast02';
A = find(LIS.Station==stn);
plot(x(2).*ones(numel(A)),LIS.Depth(A), '+k', 'LineWidth', 2, 'MarkerSize', 4);

stn = 'MID3-cast01';
A = find(LIS.Station==stn);
plot(x(3).*ones(numel(A)),LIS.Depth(A), '+k', 'LineWidth', 2, 'MarkerSize', 4);

stn = 'MID4-cast14';
A = find(LIS.Station==stn);
plot(x(4).*ones(numel(A)),LIS.Depth(A), '+k', 'LineWidth', 2, 'MarkerSize', 4);

stn = 'MID5-cast01';
A = find(LIS.Station==stn);
plot(x(5).*ones(numel(A)),LIS.Depth(A), '+k', 'LineWidth', 2, 'MarkerSize', 4);

stn = 'WLIS-cast02';
A = find(LIS.Station==stn);
plot(x(6).*ones(numel(A)),LIS.Depth(A), '+k', 'LineWidth', 2, 'MarkerSize', 4);

stn = 'WLI6-cast01';
A = find(LIS.Station==stn);
plot(x(7).*ones(numel(A)),LIS.Depth(A), '+k', 'LineWidth', 2, 'MarkerSize', 4);


xlabel('Transect distance [km]');
ylabel('Depth [m]');
caxis([ca]);
c = colorbar('location','southoutside');
c.Label.String = 'O_2 (\mumol/kg)';
c.FontSize = fs;
set(gca,'layer','top');
ylim([0 18]);
axis ij;

%print -dpng -r300 LIS_O2_transect_Aug.png;
%print -dpng -r300 LIS_O2_transect_Oct.png;

wysiwyg;