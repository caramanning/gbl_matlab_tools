load LISAug23_CH4N2O_CTD.mat;
LIS = LISAug23_CH4N2O_CTD;
LIS.Depth(isnan(LIS.Depth)) = 20;

load LISOct23_CH4N2O_CTD.mat;
LISO = LISOct23_CH4N2O_CTD;


UTC_to_local = -4/24;
LIS.datetime_local = LIS.datetime + UTC_to_local;


LIS = LISO;
%%
%
% we need to make a grid that is evenly spaced so that all the casts are
% interpolated onto the same spacing

dl = [0:0.1:20]'; % depth spacing for interpolation
x = [1 3 5 7 9 11 13 14]; % temporary x values for interpolation

n_stn = numel(x);

dl_grid = repmat(dl,1,n_stn);

x_grid = repmat(x,length(dl),1);

CH4i = nan(numel(dl),n_stn); % make a blank grid for storing CH4
N2Oi = nan(numel(dl),n_stn); % make a blank grid for storing N2O
O2i = nan(numel(dl),n_stn); % make a blank grid for storing O2
Si = nan(numel(dl),n_stn); % make a blank grid for storing S
PDeni = nan(numel(dl),n_stn); % make a blank grid for storing PDen
ti = repmat(datetime(0,0,0), 1, n_stn);

stnlist = ["MID4-cast01"
    "MID4-cast03"
    "MID4-cast05"
    "MID4-cast07"
    "MID4-cast09"
    "MID4-cast11"
    "MID4-cast13"
    "MID4-cast14"];

as = [];

for i = 1:numel(stnlist)
    A = find(LIS.Station==stnlist(i));
    as = [as; A];
    ti(i) = LIS.datetime(A(1));
    CH4i(:,i) = interp1(LIS.Depth(A),LIS.mean_CH4_nM(A),dl);
    N2Oi(:,i) = interp1(LIS.Depth(A),LIS.mean_N2O_nM(A),dl);
    Si(:,i) = interp1(LIS.Depth(A),LIS.S(A),dl);
    O2i(:,i) = interp1(LIS.Depth(A),LIS.O2_umolkg(A),dl);
    PDeni(:,i) = interp1(LIS.Depth(A),LIS.PDen(A),dl);    
end;


%%

% for OCTOBER set to EDT
% for AUGUST set to EDT

ti_EDT = ti - 4/24;
ti_EST = ti - 5/24;

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
print -dpng -r300 plot_MID4_CH4_Oct.png;


figure(3)
clf; hold on;
box on;
set(gca,'tickdir','out');
clevel = [8:0.1:12];
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
print -dpng -r300 plot_MID4_N2O_Oct.png;

figure(4)
clf; hold on;
box on;
set(gca,'tickdir','out');
%clevel = [190:260]; %Oct
clevel = [20:260]; %Aug
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
print -dpng -r300 plot_MID4_O2_Oct.png;

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
print -dpng -r300 plot_MID4_S_Oct.png;



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
print -dpng -r300 plot_MID4_PDen_Oct.png;


