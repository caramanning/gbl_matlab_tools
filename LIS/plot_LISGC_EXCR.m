load LISAug23_CH4N2O_CTD.mat;
LIS = LISAug23_CH4N2O_CTD;

load LISOct23_CH4N2O_CTD.mat;
LISO = LISOct23_CH4N2O_CTD;

LIS.CH4_eq_nM = CH4sol(LIS.S,LIS.T,1920e-9)'.*LIS.Dens./1000.*1000;
LIS.DeltaCH4 = (LIS.mean_CH4_nM - LIS.CH4_eq_nM)./LIS.CH4_eq_nM.*100;

LIS.N2O_eq_nM = N2Osol(LIS.S,LIS.T,332e-9).*LIS.Dens./1000.*1000;
LIS.DeltaN2O = (LIS.mean_N2O_nM - LIS.N2O_eq_nM)./LIS.N2O_eq_nM.*100;

LIS.O2_eq_umolkg = O2sol(LIS.S,LIS.T);
LIS.DeltaO2 = (LIS.O2_umolkg - LIS.O2_eq_umolkg)./LIS.O2_eq_umolkg.*100;


% x-axis limits
xlc = [0 500]; %ch4
xln = [7 14]; %n2o
xlo = [0 260]; %o2
xlt = [16.5 22.5]; %temp

%x-axis ticks
xtc = [0 250 500]; %

yl = [0 22.2];

figure(1)
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
s1.Position = [0.06 0.15 0.19 0.8];
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
set(gca,'xlim',xlc);
set(gca,'ylim',yl);
xticks(xtc);
%xticklabels([0 250 500]);
axis ij;


s2=subplot(1,4,2);
hold on; box on;
set(gca,'fontsize',16);
set(gca,'tickdir','out');
set(gca,'linewidth',1)
s2.Position = [0.28 0.15 0.19 0.8];
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
xlim(xln);
set(gca,'yticklabel',{[]});
set(gca,'ylim',yl);
%ylabel('Depth (m)');
axis ij;

s3=subplot(1,4,3);
hold on; box on;
set(gca,'fontsize',16);
set(gca,'tickdir','out');
set(gca,'linewidth',1)
s3.Position = [0.5 0.15 0.19 0.8];
stn = 'EXCR-cast01';
A = find(LIS.Station==stn);
plot(O2sol(LIS.S(A),LIS.T(A)),LIS.Depth(A),'--k','linewidth',1.5);
plot(LIS.O2_umolkg(A),LIS.Depth(A),'o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXCR-cast02';
A = find(LIS.Station==stn);
plot(LIS.O2_umolkg(A),LIS.Depth(A),'s-r','linewidth',1.5, 'markerfacecolor','r');
%legend('EXCR-cast01','EXCR-cast02','location','northeast')

xlabel('O_2 (\mumol/kg)');
xlim(xlo);
%ylabel('Depth (m)');
set(gca,'yticklabel',{[]});
set(gca,'ylim',yl);
axis ij;

s4=subplot(1,4,4);
hold on; box on;
set(gca,'fontsize',16);
set(gca,'tickdir','out');
set(gca,'linewidth',1)
s4.Position = [0.72 0.15 0.19 0.8];
stn = 'EXCR-cast01';
A = find(LIS.Station==stn);
plot(LIS.T(A),LIS.Depth(A),'o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXCR-cast02';
A = find(LIS.Station==stn);
plot(LIS.T(A),LIS.Depth(A),'s-r','linewidth',1.5, 'markerfacecolor','r');
%legend('EXCR-cast01','EXCR-cast02','location','northeast')

xlabel('Temp. (^oC)');
%ylabel('Depth (m)');
xlim(xlt);
set(gca,'yticklabel',{[]});
set(gca,'ylim',yl);
axis ij;

wysiwyg;

print -dpng -r300 LISAug23_EXCR_short.png;

%%

yl = [0 22.2];
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
s1.Position = [0.06 0.15 0.19 0.8];
stn = 'EXRX-cast01';

A = find(LISO.Station==stn);

errorbar(LISO.mean_CH4_nM(A),LISO.Depth(A),LISO.std_CH4_nM(A),'horizontal','o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXRX-cast02';
A = find(LISO.Station==stn);
CH4eq = CH4sol(LISO.S(A),LISO.T(A),1920e-9).*1000;
plot(CH4eq,LISO.Depth(A),'--k','linewidth',1.5);

errorbar(LISO.mean_CH4_nM(A),LISO.Depth(A),LISO.std_CH4_nM(A),'horizontal','s-r','linewidth',1.5, 'markerfacecolor','r');
%legend('EXCR-cast01','EXCR-cast02','location','northeast')

xlabel('CH_4 (nmol/kg)');
ylabel('Depth (m)');
set(gca,'xlim',xlc);
set(gca,'ylim',yl);
xticks(xtc);
%xticklabels([0 250 500]);
axis ij;


s2=subplot(1,4,2);
hold on; box on;
set(gca,'fontsize',16);
set(gca,'tickdir','out');
set(gca,'linewidth',1)
s2.Position = [0.28 0.15 0.19 0.8];
stn = 'EXRX-cast01';
A = find(LISO.Station==stn);

errorbar(LISO.mean_N2O_nM(A),LISO.Depth(A),LISO.std_N2O_nM(A),'horizontal','o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXRX-cast02';
A = find(LISO.Station==stn);

N2Oeq = N2Osol(LISO.S(A),LISO.T(A),330e-9).*1000;
plot(N2Oeq,LISO.Depth(A),'--k','linewidth',1.5);


errorbar(LISO.mean_N2O_nM(A),LISO.Depth(A),LISO.std_N2O_nM(A),'horizontal','s-r','linewidth',1.5, 'markerfacecolor','r');
%legend('EXCR-cast01','EXCR-cast02','location','northeast')

xlabel('N_2O (nmol/kg)');
xlim(xln);
set(gca,'yticklabel',{[]});
set(gca,'ylim',yl);
%ylabel('Depth (m)');
axis ij;

s3=subplot(1,4,3);
hold on; box on;
set(gca,'fontsize',16);
set(gca,'tickdir','out');
set(gca,'linewidth',1)
s3.Position = [0.5 0.15 0.19 0.8];
stn = 'EXRX-cast01';
A = find(LISO.Station==stn);
plot(LISO.O2_umolkg(A),LISO.Depth(A),'o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXRX-cast02';
A = find(LISO.Station==stn);
plot(O2sol(LISO.S(A),LISO.T(A)),LISO.Depth(A),'--k','linewidth',1.5);
plot(LISO.O2_umolkg(A),LISO.Depth(A),'s-r','linewidth',1.5, 'markerfacecolor','r');
%legend('EXCR-cast01','EXCR-cast02','location','northeast')

xlabel('O_2 (\mumol/kg)');
xlim(xlo);
set(gca,'ylim',yl);
%ylabel('Depth (m)');
set(gca,'yticklabel',{[]})
axis ij;

s4=subplot(1,4,4);
hold on; box on;
set(gca,'fontsize',16);
set(gca,'tickdir','out');
set(gca,'linewidth',1)
s4.Position = [0.72 0.15 0.19 0.8];
stn = 'EXRX-cast01';
A = find(LISO.Station==stn);
plot(LISO.T(A),LISO.Depth(A),'o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXRX-cast02';
A = find(LISO.Station==stn);
plot(LISO.T(A),LISO.Depth(A),'s-r','linewidth',1.5, 'markerfacecolor','r');
%legend('EXCR-cast01','EXCR-cast02','location','northeast')

xlabel('Temp. (^oC)');
%ylabel('Depth (m)');
xlim(xlt);
set(gca,'yticklabel',{[]});
set(gca,'ylim',yl);
axis ij;

wysiwyg;

print -dpng -r300 LISOct23_EXCR_short.png;

%%
yl = [0 22.5]

figure(2)
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
stn = 'EXRX-cast01';

A = find(LISO.Station==stn);


errorbar(LISO.mean_CH4_nM(A),LISO.Depth(A),LISO.std_CH4_nM(A),'horizontal','o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXRX-cast02';
A = find(LISO.Station==stn);

CH4eq = CH4sol(LISO.S(A),LISO.T(A),1920e-9).*1000;
plot(CH4eq,LISO.Depth(A),'--k','linewidth',1.5);

errorbar(LISO.mean_CH4_nM(A),LISO.Depth(A),LISO.std_CH4_nM(A),'horizontal','s-r','linewidth',1.5, 'markerfacecolor','r');
%legend('EXCR-cast01','EXCR-cast02','location','northeast')

xlabel('CH_4 (nmol/kg)');
ylabel('Depth (m)');
set(gca,'xlim',xlc);
set(gca,'ylim',yl);
xticks(xtc);
axis ij;

subplot(1,4,2)
hold on; box on;
set(gca,'fontsize',16);
stn = 'EXRX-cast01';
A = find(LISO.Station==stn);


errorbar(LISO.mean_N2O_nM(A),LISO.Depth(A),LISO.std_N2O_nM(A),'horizontal','o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXRX-cast02';
A = find(LISO.Station==stn);

N2Oeq = N2Osol(LISO.S(A),LISO.T(A),330e-9).*1000;
plot(N2Oeq,LISO.Depth(A),'--k','linewidth',1.5);

errorbar(LISO.mean_N2O_nM(A),LISO.Depth(A),LISO.std_N2O_nM(A),'horizontal','s-r','linewidth',1.5, 'markerfacecolor','r');
%legend('EXCR-cast01','EXCR-cast02','location','northeast')

xlabel('N_2O (nmol/kg)');
xlim(xln);
set(gca,'ylim',yl);
%ylabel('Depth (m)');
axis ij;


subplot(1,4,3)
hold on; box on;
set(gca,'fontsize',16);
stn = 'EXRX-cast01';

A = find(LISO.Station==stn);
plot(LISO.O2_umolkg(A),LISO.Depth(A),'o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXRX-cast02';
A = find(LISO.Station==stn);
plot(O2sol(LISO.S(A),LISO.T(A)),LISO.Depth(A),'--k','linewidth',1.5);
plot(LISO.O2_umolkg(A),LISO.Depth(A),'s-r','linewidth',1.5, 'markerfacecolor','r');
%legend('EXCR-cast01','EXCR-cast02','location','northeast')

xlabel('O_2 (\mumol/kg)');
xlim(xlo);
set(gca,'ylim',yl);
%ylabel('Depth (m)');
axis ij;

subplot(1,4,4)
hold on; box on;
set(gca,'fontsize',16);
stn = 'EXRX-cast01';
A = find(LISO.Station==stn);
plot(LISO.T(A),LISO.Depth(A),'o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXRX-cast02';
A = find(LISO.Station==stn);
plot(LISO.T(A),LISO.Depth(A),'s-r','linewidth',1.5, 'markerfacecolor','r');
%legend('EXCR-cast01','EXCR-cast02','location','northeast')

xlabel('Temp. (^oC)');
%ylabel('Depth (m)');
xlim(xlt);
set(gca,'ylim',yl);
axis ij;

wysiwyg;

print -dpng -r300 LISOct23_EXCR.png;


