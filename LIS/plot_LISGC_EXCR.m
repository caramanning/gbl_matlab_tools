load LISAug23_CH4N2O_CTD.mat
LIS = LISAug23_CH4N2O_CTD;



figure(1)
clf; hold on;

clf; 

%sp=tight_subplot(nr,nc,[0.12 .03],[.15 .06],[.1 .04]);
hold on; box on;
%set(gcf,'color','w');
set(gcf, 'PaperUnits', 'inches');
set(gcf,'renderer','painters');
set(gcf, 'PaperPosition', [0 0 12 6]);
set(gca,'linewidth',0.75);
    set(gcf,'GraphicsSmoothing','on')

subplot(1,4,1)
hold on; box on;
set(gca,'fontsize',16);
stn = 'EXCR-cast01';

A = find(LIS.Station==stn);

CH4eq = CH4sol(LIS.S(A),LIS.T(A),1920e-9).*1000;
plot(CH4eq,LIS.Depth(A),'--k','linewidth',1.5);

errorbar(LIS.mean_CH4_nM(A),LIS.Depth(A),LIS.std_CH4_nM(A),'horizontal','o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXCR-cast02';
A = find(LIS.Station==stn);
errorbar(LIS.mean_CH4_nM(A),LIS.Depth(A),LIS.std_CH4_nM(A),'horizontal','s-r','linewidth',1.5, 'markerfacecolor','r');
%legend('EXCR-cast01','EXCR-cast02','location','northeast')

xlabel('CH_4 (nmol/kg)');
ylabel('Depth (m)');
set(gca,'xlim',[0 500]);
axis ij;


subplot(1,4,2)
hold on; box on;
set(gca,'fontsize',16);
stn = 'EXCR-cast01';
A = find(LIS.Station==stn);

N2Oeq = N2Osol(LIS.S(A),LIS.T(A),330e-9).*1000;

plot(N2Oeq,LIS.Depth(A),'--k','linewidth',1.5);
errorbar(LIS.mean_N2O_nM(A),LIS.Depth(A),LIS.std_N2O_nM(A),'horizontal','o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXCR-cast02';
A = find(LIS.Station==stn);
errorbar(LIS.mean_N2O_nM(A),LIS.Depth(A),LIS.std_N2O_nM(A),'horizontal','s-r','linewidth',1.5, 'markerfacecolor','r');
%legend('EXCR-cast01','EXCR-cast02','location','northeast')

xlabel('N_2O (nmol/kg)');
xlim([7 12]);
%ylabel('Depth (m)');
axis ij;


subplot(1,4,3)
hold on; box on;
set(gca,'fontsize',16);
stn = 'EXCR-cast01';
A = find(LIS.Station==stn);
plot(O2sol(LIS.S(A),LIS.T(A)),LIS.Depth(A),'--k','linewidth',1.5);
plot(LIS.O2_umolkg(A),LIS.Depth(A),'o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXCR-cast02';
A = find(LIS.Station==stn);
plot(LIS.O2_umolkg(A),LIS.Depth(A),'s-r','linewidth',1.5, 'markerfacecolor','r');
%legend('EXCR-cast01','EXCR-cast02','location','northeast')

xlabel('O_2 (\mumol/kg)');
xlim([0 250]);
%ylabel('Depth (m)');
axis ij;

subplot(1,4,4)
hold on; box on;
set(gca,'fontsize',16);
stn = 'EXCR-cast01';
A = find(LIS.Station==stn);
plot(LIS.T(A),LIS.Depth(A),'o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXCR-cast02';
A = find(LIS.Station==stn);
plot(LIS.T(A),LIS.Depth(A),'s-r','linewidth',1.5, 'markerfacecolor','r');
%legend('EXCR-cast01','EXCR-cast02','location','northeast')

xlabel('Temp. (^oC)');
%ylabel('Depth (m)');
axis ij;

wysiwyg;

print -dpng -r300 LISAug23_EXCR.png;


%%
% we need to make a grid that is evenly spaced so that all the casts are
% interpolated onto the same spacing

dl = [0:0.5:20]'; % depth spacing for interpolation
x = [1 3 5 7 9 11 13 14]; % temporary x values for interpolation
n_stn = numel(x);

dl_grid = repmat(dl,1,n_stn);

x_grid = repmat(x,length(dl),1);

CH4i = nan(numel(dl),n_stn); % make a blank grid for storing CH4
N2Oi = nan(numel(dl),n_stn); % make a blank grid for storing CH4


stn = 'MID4-cast01';
A = find(LIS.Station==stn);
CH4i(:,1) = interp1(LIS.Depth(A),LIS.mean_CH4_nM(A),dl);
N2Oi(:,1) = interp1(LIS.Depth(A),LIS.mean_N2O_nM(A),dl);


stn = 'MID4-cast03';
A = find(LIS.Station==stn);
CH4i(:,2) = interp1(LIS.Depth(A),LIS.mean_CH4_nM(A),dl);
N2Oi(:,2) = interp1(LIS.Depth(A),LIS.mean_N2O_nM(A),dl);


stn = 'MID4-cast05';
A = find(LIS.Station==stn);
CH4i(:,3) = interp1(LIS.Depth(A),LIS.mean_CH4_nM(A),dl);
N2Oi(:,3) = interp1(LIS.Depth(A),LIS.mean_N2O_nM(A),dl);

stn = 'MID4-cast07';
A = find(LIS.Station==stn);
CH4i(:,4) = interp1(LIS.Depth(A),LIS.mean_CH4_nM(A),dl);
N2Oi(:,4) = interp1(LIS.Depth(A),LIS.mean_N2O_nM(A),dl);


stn = 'MID4-cast09';
A = find(LIS.Station==stn);
CH4i(:,5) = interp1(LIS.Depth(A),LIS.mean_CH4_nM(A),dl);
N2Oi(:,5) = interp1(LIS.Depth(A),LIS.mean_N2O_nM(A),dl);

stn = 'MID4-cast11';
A = find(LIS.Station==stn);
CH4i(:,6) = interp1(LIS.Depth(A),LIS.mean_CH4_nM(A),dl);
N2Oi(:,6) = interp1(LIS.Depth(A),LIS.mean_N2O_nM(A),dl);

stn = 'MID4-cast13';
A = find(LIS.Station==stn);
CH4i(:,7) = interp1(LIS.Depth(A),LIS.mean_CH4_nM(A),dl);
N2Oi(:,7) = interp1(LIS.Depth(A),LIS.mean_N2O_nM(A),dl);

stn = 'MID4-cast14';
A = find(LIS.Station==stn);
CH4i(:,8) = interp1(LIS.Depth(A),LIS.mean_CH4_nM(A),dl);
N2Oi(:,8) = interp1(LIS.Depth(A),LIS.mean_N2O_nM(A),dl);


figure(2)
clf; hold on;
box on;
set(gca,'tickdir','out');
clevel = [40:1:80];
C = contourf(x_grid,dl_grid,CH4i,clevel,'edgecolor','none');
xlabel('cast number');
ylabel('Depth (m)');
c = colorbar;
c.Label.String = 'CH_4 (nM)';
c.Label.FontSize = 16;

axis ij;

figure(3)
clf; hold on;
box on;
set(gca,'tickdir','out');
clevel = [8:0.1:12];
C = contourf(x_grid,dl_grid,N2Oi,clevel,'edgecolor','none');
xlabel('cast number');
ylabel('Depth (m)');
c = colorbar;
c.Label.String = 'N_2O (nM)';
c.Label.FontSize = 16;

axis ij;


