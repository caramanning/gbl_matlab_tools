% load in August 2023 data
load LISAug23_CH4N2O_CTD.mat
LISA = LISAug23_CH4N2O_CTD; % August

load LISAug2023CastData.mat;
LISCDA.Cast = LISAug2023CastData.Cast;
LISCDA.Lat = LISAug2023CastData.Lat;
LISCDA.Lon = LISAug2023CastData.Lon;

ncast = 16:22;
km_between = m_lldist(LISCDA.Lon(ncast),LISCDA.Lat(ncast)); % distance between consecutive stations
km_between = [0; km_between]; % add on zero for first station

km_cumulativeA = cumsum(km_between); % consecutive distance



load LISOct23_CH4N2O_CTD.mat
LISO = LISOct23_CH4N2O_CTD; % October
LIS = LISOct23_CH4N2O_CTD;

load LISOct2023CastData.mat;
LISCDO.Cast = LISOct2023CastData.Cast;
LISCDO.Lat = LISOct2023CastData.Lat;
LISCDO.Lon = LISOct2023CastData.Lon;

ncast = 16:22;
Cast = LISO.Cast(ncast);
Lat = LISO.Lat(ncast);
Lon = LISO.Lon(ncast);

ncast = 16:22; % casts for transect;
km_between = m_lldist(LISCDO.Lon(ncast),LISCDO.Lat(ncast)); % distance between consecutive stations
km_between = [0; km_between]; % add on zero for first station

km_cumulativeO = cumsum(km_between); % consecutive distance


% for October the station names are different, need to change later
stn = 'EXRX-cast01';
A = find(LIS.Station==stn);
LISO.Station(A) = 'EXCR-cast01';
LIS.Station(A) = 'EXCR-cast01';

stn = 'EXRX-cast02';
A = find(LIS.Station==stn);
LISO.Station(A) = 'EXCR-cast02';
LIS.Station(A) = 'EXCR-cast02';

stn = 'EXCR1-cast01';
A = find(LIS.Station==stn);
LISO.Station(A) = 'EXR1-cast01';
LIS.Station(A) = 'EXR1-cast01';

%%
% we need to make a grid that is evenly spaced so that all the casts are
% interpolated onto the same spacing

dl = [0:0.1:18]'; % depth spacing for interpolation
%x = km_cumulativeA'; % temporary x values for interpolation
x = [-0.5 km_cumulativeA' 19]; % now set values to have a bit of extra start and end
n_stn = numel(x);

dl_grid = repmat(dl,1,n_stn);

% modify the station list to add a bit of distance extending past the first
% station
stnlist = ["EXR1-cast01"
    "EXR1-cast01"
    "EXCR-cast02"
    "MID3-cast01"
    "MID4-cast14"
    "MID5-cast01"
    "WLIS-cast02"
    "WLI6-cast01"
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
N2OiA = nan(numel(dl),n_stn); % make a blank grid for storing CH4
O2iA = nan(numel(dl),n_stn); % make a blank grid for storing CH4

CH4iA = nan(numel(dl),n_stn); % make a blank grid for storing CH4
N2Oia = nan(numel(dl),n_stn); % make a blank grid for storing N2O
O2iA = nan(numel(dl),n_stn); % make a blank grid for storing O2
SiA = nan(numel(dl),n_stn); % make a blank grid for storing S
PDeniA = nan(numel(dl),n_stn); % make a blank grid for storing PDen
tiA = repmat(datetime(0,0,0), 1, n_stn);

for i = 1:numel(stnlist)
    A = find(LISA.Station==stnlist(i));
    asA = [asA; A];
    if (i ~=1 && i~= numel(stnlist))
        LISA.km_cumulative(A) = x(i); % make sure we only do this for the real stations and not fillers
    end;

    tiA(i) = LISA.datetime(A(1));
    CH4iA(:,i) = interp1(LISA.Depth(A),LISA.mean_CH4_nM(A),dl);
    N2OiA(:,i) = interp1(LISA.Depth(A),LISA.mean_N2O_nM(A),dl);
    SiA(:,i) = interp1(LISA.Depth(A),LISA.S(A),dl);
    O2iA(:,i) = interp1(LISA.Depth(A),LISA.O2_umolkg(A),dl);
    PDeniA(:,i) = interp1(LISA.Depth(A),LISA.PDen(A),dl);    
end;

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
N2OiO = nan(numel(dl),n_stn); % make a blank grid for storing CH4
O2iO = nan(numel(dl),n_stn); % make a blank grid for storing CH4

CH4iO = nan(numel(dl),n_stn); % make a blank grid for storing CH4
N2Oia = nan(numel(dl),n_stn); % make a blank grid for storing N2O
O2iO = nan(numel(dl),n_stn); % make a blank grid for storing O2
SiO = nan(numel(dl),n_stn); % make a blank grid for storing S
PDeniO = nan(numel(dl),n_stn); % make a blank grid for storing PDen
tiO = repmat(datetime(0,0,0), 1, n_stn);

for i = 1:numel(stnlist)
    s = find(LISO.Station==stnlist(i));
    asO = [asO; s];
    if (i ~=1 && i~= numel(stnlist))
        LISO.km_cumulative(s) = x(i); % make sure we only do this for the real stations and not fillers
    end;

    tiO(i) = LISO.datetime(s(1));
    CH4iO(:,i) = interp1(LISO.Depth(s),LISO.mean_CH4_nM(s),dl);
    N2OiO(:,i) = interp1(LISO.Depth(s),LISO.mean_N2O_nM(s),dl);
    SiO(:,i) = interp1(LISO.Depth(s),LISO.S(s),dl);
    O2iO(:,i) = interp1(LISO.Depth(s),LISO.O2_umolkg(s),dl);
    PDeniO(:,i) = interp1(LISO.Depth(s),LISO.PDen(s),dl);    
end;

% MAY 2024 will be added in here once processed


%%
nr = 4; % number of rows
nc = 3; % number of columns
lw = 2; % default line width
fs = 10; % default font size
ms = 3; % default marker size
yt = [0 5 10 15 20]; %y-axis ticks;
xt = [0 5 10 15];
yl = [0 20]; % y-axis limit in m
xl = [-0.5 19]; %x-axis limit in km
caCH4 = [25 460]; % CH4 axis limits
clevelCH4 = [25:1:460]; %CH4 colorbar levels

caN2O = [8.5 14]; % N2O axis limits
clevelN2O = [8.5:0.02:14]; % N2O colorbar levels

caO2 = [40 250];
clevelO2 = [40:1:250];

caPDen = [17.35 18.8]; % N2O axis limits
clevelPDen = [17.35:0.02:18.8]; % N2O colorbar levels

fig=figure(23);
clf;
sp=tight_subplot(nr,nc,[.025 .025],[.08 .04],[.08 .04]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 6]);
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
    axis ij;
    title('October');

% SUBPLOT 3
subplot(sp(3))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    clevel = [25:1:460];
    ca = [25 460]; %colorbar limits
    C = contourf(x_gridO,dl_grid,CH4iO,clevelCH4,'edgecolor','none');
    plot(LISO.km_cumulative(asO),LISO.Depth(asO),'+k', 'linewidth', 2, 'markersize',ms)

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
    axis ij;
    title(['May - fake data'])

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
    axis ij;
    %title('October');

% SUBPLOT 6
subplot(sp(6))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridO,dl_grid,N2OiO,clevelN2O,'edgecolor','none');
    plot(LISO.km_cumulative(asO),LISO.Depth(asO),'+k', 'linewidth', 2, 'markersize',ms)

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
    axis ij;
    %title(['May - temp fake data'])

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
    axis ij;
    %title('October');

% SUBPLOT 9
subplot(sp(9))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridO,dl_grid,O2iO,clevelO2,'edgecolor','none');
    plot(LISO.km_cumulative(asO),LISO.Depth(asO),'+k', 'linewidth', 2, 'markersize',ms)

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
    axis ij;
    %title(['May - temp fake data'])


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
    axis ij;
    %title('October');

% SUBPLOT 12
subplot(sp(12))
    hold on; box on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',fs);

    C = contourf(x_gridO,dl_grid,PDeniO,clevelPDen,'edgecolor','none');
  %  C = surf(x_grid,dl_grid,PDeniO,clevelPDen,'facecolor','interp','edgecolor','interp');
  %  surf instead, with 'FaceColor','interp', 'EdgeColor','interp' and %view(0,90)’.  
    plot(LISO.km_cumulative(asO),LISO.Depth(asO),'+k', 'linewidth', 2, 'markersize',ms)

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
    axis ij;
    %title(['May - temp fake data'])

    
    wysiwyg;

print(gcf, '-dpng', '-r300', 'LIS_Transect_CH4_N2O_O2_PDen.png');
print(gcf,'-depsc','-vector','LIS_Transect_CH4_N2O_O2_PDen.eps');
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
print -dpng -r300 LIS_O2_transect_Oct.png;

wysiwyg;