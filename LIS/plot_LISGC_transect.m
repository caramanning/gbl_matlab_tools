% load LISAug23_CH4N2O_CTD.mat
% LIS = LISAug23_CH4N2O_CTD;
% 
% load LISAug2023CastData.mat;
% 
% ncast = 16:22;
% Cast = LISAug2023CastData.Cast(ncast);
% Lat = LISAug2023CastData.Lat(ncast);
% Lon = LISAug2023CastData.Lon(ncast);


load LISOct23_CH4N2O_CTD.mat
LIS = LISOct23_CH4N2O_CTD;

load LISOct2023CastData.mat;

ncast = 16:22;
Cast = LISOct2023CastData.Cast(ncast);
Lat = LISOct2023CastData.Lat(ncast);
Lon = LISOct2023CastData.Lon(ncast);


km_between = m_lldist(Lon,Lat); % distance between consecutive stations
km_between = [0; km_between]; % add on zero for first station

km_cumulative = cumsum(km_between); % consecutive distance

% for October the station names are different
stn = 'EXRX-cast01';
A = find(LIS.Station==stn);
LIS.Station(A) = 'EXCR-cast01';

stn = 'EXRX-cast02';
A = find(LIS.Station==stn);
LIS.Station(A) = 'EXCR-cast02';

stn = 'EXCR1-cast01';
A = find(LIS.Station==stn);
LIS.Station(A) = 'EXR1-cast01';

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


stn = 'EXR1-cast01';
A = find(LIS.Station==stn);
CH4i(:,1) = interp1(LIS.Depth(A),LIS.mean_CH4_nM(A),dl);
N2Oi(:,1) = interp1(LIS.Depth(A),LIS.mean_N2O_nM(A),dl);
O2i(:,1) = interp1(LIS.Depth(A),LIS.O2_umolkg(A),dl);


stn = 'EXCR-cast02';
A = find(LIS.Station==stn);
CH4i(:,2) = interp1(LIS.Depth(A),LIS.mean_CH4_nM(A),dl);
N2Oi(:,2) = interp1(LIS.Depth(A),LIS.mean_N2O_nM(A),dl);
O2i(:,2) = interp1(LIS.Depth(A),LIS.O2_umolkg(A),dl);


stn = 'MID3-cast01';
A = find(LIS.Station==stn);
CH4i(:,3) = interp1(LIS.Depth(A),LIS.mean_CH4_nM(A),dl);
N2Oi(:,3) = interp1(LIS.Depth(A),LIS.mean_N2O_nM(A),dl);
O2i(:,3) = interp1(LIS.Depth(A),LIS.O2_umolkg(A),dl);

stn = 'MID4-cast14';
A = find(LIS.Station==stn);
CH4i(:,4) = interp1(LIS.Depth(A),LIS.mean_CH4_nM(A),dl);
N2Oi(:,4) = interp1(LIS.Depth(A),LIS.mean_N2O_nM(A),dl);
O2i(:,4) = interp1(LIS.Depth(A),LIS.O2_umolkg(A),dl);

stn = 'MID5-cast01';
A = find(LIS.Station==stn);
CH4i(:,5) = interp1(LIS.Depth(A),LIS.mean_CH4_nM(A),dl);
N2Oi(:,5) = interp1(LIS.Depth(A),LIS.mean_N2O_nM(A),dl);
O2i(:,5) = interp1(LIS.Depth(A),LIS.O2_umolkg(A),dl);

stn = 'WLIS-cast02';
A = find(LIS.Station==stn);
CH4i(:,6) = interp1(LIS.Depth(A),LIS.mean_CH4_nM(A),dl);
N2Oi(:,6) = interp1(LIS.Depth(A),LIS.mean_N2O_nM(A),dl);
O2i(:,6) = interp1(LIS.Depth(A),LIS.O2_umolkg(A),dl);


stn = 'WLI6-cast01';
A = find(LIS.Station==stn);
[~,Ai] = sort(LIS.Depth(A));

CH4i(:,7) = interp1(LIS.Depth(A(Ai)),LIS.mean_CH4_nM(A(Ai)),dl);
N2Oi(:,7) = interp1(LIS.Depth(A(Ai)),LIS.mean_N2O_nM(A(Ai)),dl);
O2i(:,7) = interp1(LIS.Depth(A),LIS.O2_umolkg(A),dl);



fs=22;
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
c.Label.String = 'N_2O (nmol/kg)';
c.FontSize = fs;
ylim([0 18]);
set(gca,'layer','top');
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


