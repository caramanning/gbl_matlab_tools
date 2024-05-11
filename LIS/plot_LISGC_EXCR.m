load LISAug23_CH4N2O_CTD.mat;
LIS = LISAug23_CH4N2O_CTD;

load LISOct23_CH4N2O_CTD.mat;
LISO = LISOct23_CH4N2O_CTD;


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
xlim([7 12.5]);
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

A = find(LISO.Station==stn);

CH4eq = CH4sol(LISO.S(A),LISO.T(A),1920e-9).*1000;
plot(CH4eq,LISO.Depth(A),'--k','linewidth',1.5);

errorbar(LISO.mean_CH4_nM(A),LISO.Depth(A),LISO.std_CH4_nM(A),'horizontal','o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXCR-cast02';
A = find(LISO.Station==stn);
errorbar(LISO.mean_CH4_nM(A),LISO.Depth(A),LISO.std_CH4_nM(A),'horizontal','s-r','linewidth',1.5, 'markerfacecolor','r');
%legend('EXCR-cast01','EXCR-cast02','location','northeast')

xlabel('CH_4 (nmol/kg)');
ylabel('Depth (m)');
set(gca,'xlim',[0 500]);
axis ij;


subplot(1,4,2)
hold on; box on;
set(gca,'fontsize',16);
stn = 'EXCR-cast01';
A = find(LISO.Station==stn);

N2Oeq = N2Osol(LISO.S(A),LISO.T(A),330e-9).*1000;

plot(N2Oeq,LISO.Depth(A),'--k','linewidth',1.5);
errorbar(LISO.mean_N2O_nM(A),LISO.Depth(A),LISO.std_N2O_nM(A),'horizontal','o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXCR-cast02';
A = find(LISO.Station==stn);
errorbar(LISO.mean_N2O_nM(A),LISO.Depth(A),LISO.std_N2O_nM(A),'horizontal','s-r','linewidth',1.5, 'markerfacecolor','r');
%legend('EXCR-cast01','EXCR-cast02','location','northeast')

xlabel('N_2O (nmol/kg)');
xlim([7 12.5]);
%ylabel('Depth (m)');
axis ij;


subplot(1,4,3)
hold on; box on;
set(gca,'fontsize',16);
stn = 'EXCR-cast01';
A = find(LISO.Station==stn);
plot(O2sol(LISO.S(A),LISO.T(A)),LISO.Depth(A),'--k','linewidth',1.5);
plot(LISO.O2_umolkg(A),LISO.Depth(A),'o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXCR-cast02';
A = find(LISO.Station==stn);
plot(LISO.O2_umolkg(A),LISO.Depth(A),'s-r','linewidth',1.5, 'markerfacecolor','r');
%legend('EXCR-cast01','EXCR-cast02','location','northeast')

xlabel('O_2 (\mumol/kg)');
xlim([0 250]);
%ylabel('Depth (m)');
axis ij;

subplot(1,4,4)
hold on; box on;
set(gca,'fontsize',16);
stn = 'EXCR-cast01';
A = find(LISO.Station==stn);
plot(LISO.T(A),LISO.Depth(A),'o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXCR-cast02';
A = find(LISO.Station==stn);
plot(LISO.T(A),LISO.Depth(A),'s-r','linewidth',1.5, 'markerfacecolor','r');
%legend('EXCR-cast01','EXCR-cast02','location','northeast')

xlabel('Temp. (^oC)');
%ylabel('Depth (m)');
axis ij;

wysiwyg;

print -dpng -r300 LISOct23_EXCR.png;


