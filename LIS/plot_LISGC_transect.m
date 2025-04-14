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
    asA = [];
    for i = 1:length(stnlist)
        A = find(LISA.Station==stnlist(i));
        asA = [asA; A];
    end;

    CH4A = LISA.CH4_mean_nmolkg(asA);
    O2A = LISA.O2_umolkg(asA);

    [min(CH4A) max(CH4A) mean(CH4A) median(CH4A)]
    [min(O2A) max(O2A) mean(O2A) median(O2A)]

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

    asM = [];
    for i = 1:length(stnlist)
        M = find(LISM.Station==stnlist(i));
        asM = [asM; M];
    end;

    CH4M = LISM.CH4_mean_nmolkg(asM);
    O2M = LISM.O2_umolkg(asM);

    [min(CH4M) max(CH4M) mean(CH4M) median(CH4M)]
    [min(O2M) max(O2M) mean(O2M) median(O2M)]


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

%for i = 1:numel(stnlist)
for i = 1:numel(stnlist)    
    A = find(LISA.Station==stnlist(i));
    asA = [asA; A];
    LISA.km_cumulative(A) = x(i+1); % make sure we only do this for the real stations and not fillers

    tiA(i+1) = LISA.datetime(A(1));
    CH4iA(:,i+1) = interp1(LISA.Depth(A),LISA.mean_CH4_nM(A),dl);
    N2OiA(:,i+1) = interp1(LISA.Depth(A),LISA.mean_N2O_nM(A),dl);
    SiA(:,i+1) = interp1(LISA.Depth(A),LISA.S(A),dl);
    O2iA(:,i+1) = interp1(LISA.Depth(A),LISA.O2_umolkg(A),dl);
    PDeniA(:,i+1) = interp1(LISA.Depth(A),LISA.PDen(A),dl); 

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
end;

dmin = 1.6; % minimum depth to plot 
CH4iA(dl_gridA<dmin) = NaN;
N2OiA(dl_gridA<dmin) = NaN;
SiA(dl_gridA<dmin) = NaN;
O2iA(dl_gridA<dmin) = NaN;
PDeniA(dl_gridA<dmin) = NaN;


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
    s = find(LISO.Station==stnlist(i));
    asO = [asO; s];
    LISO.km_cumulative(s) = x(i+1); % make sure we only do this for the real stations and not fillers

    tiO(i+1) = LISO.datetime(s(1));
    CH4iO(:,i+1) = interp1(LISO.Depth(s),LISO.mean_CH4_nM(s),dl);
    N2OiO(:,i+1) = interp1(LISO.Depth(s),LISO.mean_N2O_nM(s),dl);
    SiO(:,i+1) = interp1(LISO.Depth(s),LISO.S(s),dl);
    O2iO(:,i+1) = interp1(LISO.Depth(s),LISO.O2_umolkg(s),dl);
    PDeniO(:,i+1) = interp1(LISO.Depth(s),LISO.PDen(s),dl); 

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
end;

dmin = 1.3; % minimum depth to plot 
CH4iO(dl_gridA<dmin) = NaN;
N2OiO(dl_gridA<dmin) = NaN;
SiO(dl_gridA<dmin) = NaN;
O2iO(dl_gridA<dmin) = NaN;
PDeniO(dl_gridA<dmin) = NaN;

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
    s = find(LISM.Station==stnlist(i));
    asM = [asM; s];
    LISM.km_cumulative(s) = x(i+1); % make sure we only do this for the real stations and not fillers

    tiM(i+1) = LISM.datetime(s(1));
    CH4iM(:,i+1) = interp1(LISM.Depth(s),LISM.mean_CH4_nM(s),dl);
    N2OiM(:,i+1) = interp1(LISM.Depth(s),LISM.mean_N2O_nM(s),dl);
    SiM(:,i+1) = interp1(LISM.Depth(s),LISM.S(s),dl);
    O2iM(:,i+1) = interp1(LISM.Depth(s),LISM.O2_umolkg(s),dl);
    PDeniM(:,i+1) = interp1(LISM.Depth(s),LISM.PDen(s),dl); 

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
end;

dmin = 1.3; % minimum depth to plot 
CH4iM(dl_gridA<dmin) = NaN;
N2OiM(dl_gridA<dmin) = NaN;
SiM(dl_gridA<dmin) = NaN;
O2iM(dl_gridA<dmin) = NaN;
PDeniM(dl_gridA<dmin) = NaN;

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

caPDen = [15.7 18.8]; % PDen axis limits
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