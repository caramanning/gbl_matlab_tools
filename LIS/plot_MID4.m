%load LISAug23_CH4N2O_CTD.mat;
%LIS = LISAug23_CH4N2O_CTD;
%LIS.Depth(isnan(LIS.Depth)) = 20;

load LISOct23_CH4N2O_CTD.mat;
LIS = LISOct23_CH4N2O_CTD;

load LISMay24_CH4N2O_CTD.mat;
LIS = LISMay24_CH4N2O_CTD;


UTC_to_local = -4/24;
LIS.datetime_local = LIS.datetime + UTC_to_local;


%LIS = LISO;
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
% set values above the first sample to the value of the top sample
sdl = find(dl>1.5,1);
for i = 1:length(CH4i(1,:))
    % find first non-NaN value
    fnn = find(~isnan(CH4i(:,i)),1);
    CH4i(sdl:fnn-1,i) = CH4i(fnn,i);
    N2Oi(sdl:fnn-1,i) = N2Oi(fnn,i);
    Si(sdl:fnn-1,i) = Si(fnn,i);
    O2i(sdl:fnn-1,i) = O2i(fnn,i);
    PDeni(sdl:fnn-1,i) = PDeni(fnn,i);
end;
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
