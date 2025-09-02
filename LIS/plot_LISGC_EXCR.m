load LISAug23_CH4N2O_CTD.mat;
LISA = LISAug23_CH4N2O_CTD;

load LISOct23_CH4N2O_CTD.mat;
LISO = LISOct23_CH4N2O_CTD;

load LISMay24_CH4N2O_CTD.mat;
LISM = LISMay24_CH4N2O_CTD;

LISA.CH4_eq_nM = CH4sol(LISA.S,LISA.T,1920e-9)'.*LISA.Dens./1000.*1000;
LISA.DeltaCH4 = (LISA.mean_CH4_nM - LISA.CH4_eq_nM)./LISA.CH4_eq_nM.*100;

LISA.N2O_eq_nM = N2Osol(LISA.S,LISA.T,332e-9).*LISA.Dens./1000.*1000;
LISA.DeltaN2O = (LISA.mean_N2O_nM - LISA.N2O_eq_nM)./LISA.N2O_eq_nM.*100;

LISA.O2_eq_umolkg = O2sol(LISA.S,LISA.T);
LISA.DeltaO2 = (LISA.O2_umolkg - LISA.O2_eq_umolkg)./LISA.O2_eq_umolkg.*100;

LISO.CH4_eq_nM = CH4sol(LISO.S,LISO.T,1920e-9)'.*LISO.Dens./1000.*1000;
LISO.DeltaCH4 = (LISO.mean_CH4_nM - LISO.CH4_eq_nM)./LISO.CH4_eq_nM.*100;

LISO.N2O_eq_nM = N2Osol(LISO.S,LISO.T,332e-9).*LISO.Dens./1000.*1000;
LISO.DeltaN2O = (LISO.mean_N2O_nM - LISO.N2O_eq_nM)./LISO.N2O_eq_nM.*100;

LISO.O2_eq_umolkg = O2sol(LISO.S,LISO.T);
LISO.DeltaO2 = (LISO.O2_umolkg - LISO.O2_eq_umolkg)./LISO.O2_eq_umolkg.*100;

LISM.CH4_eq_nM = CH4sol(LISM.S,LISM.T,1920e-9)'.*LISM.Dens./1000.*1000;
LISM.DeltaCH4 = (LISM.mean_CH4_nM - LISM.CH4_eq_nM)./LISM.CH4_eq_nM.*100;

LISM.N2O_eq_nM = N2Osol(LISM.S,LISM.T,332e-9).*LISM.Dens./1000.*1000;
LISM.DeltaN2O = (LISM.mean_N2O_nM - LISM.N2O_eq_nM)./LISM.N2O_eq_nM.*100;

LISM.O2_eq_umolkg = O2sol(LISM.S,LISM.T);
LISM.DeltaO2 = (LISM.O2_umolkg - LISM.O2_eq_umolkg)./LISM.O2_eq_umolkg.*100;

%%
% x-axis limits
xlc = [0 500]; %ch4
xln = [9 15.5]; %n2o
xlo = [0 350]; %o2
xlt = [11 22.5]; %temp
xld = [15.5 18.5];
fs = 12; % font size

nr=4; nc=3;

%x-axis ticks
xtc = [0 250 500]; %
xtn = [10 12 14];
xto = [0 100 200 300];
xtd = [16 17 18];

yl = [0 24];
yt = [0 10 20];

figure(6)
clf;

sp=tight_subplot(nr,nc,[.08 .03],[.08 .04],[.08 .04]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 8]);
set(gcf,'renderer','painters');
    set(gcf,'GraphicsSmoothing','on')

% ------- AUG CH4    
subplot(sp(1));
hold on; box on;
set(gca,'fontsize',fs);
set(gca,'tickdir','out');
set(gca,'linewidth',1);
title('August 2023');
%s1.Position = [0.06 0.15 0.19 0.8];
stn = 'EXRX-cast01';

A = find(LISA.Station==stn);
[~,B] = sort(LISA.Depth(A));
CH4eq = CH4sol(LISA.S(A(B)),LISA.T(A(B)),1920e-9).*1000;
%plot(CH4eq,LISA.Depth(A(B)),'--k','linewidth',1.5);

errorbar(LISA.mean_CH4_nM(A(B)),LISA.Depth(A(B)),LISA.std_CH4_nM(A(B)),'horizontal','o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXRX-cast02';
A = find(LISA.Station==stn);
[~,B] = sort(LISA.Depth(A));
errorbar(LISA.mean_CH4_nM(A(B)),LISA.Depth(A(B)),LISA.std_CH4_nM(A(B)),'horizontal','s-r','linewidth',1.5, 'markerfacecolor','r');
%legend('EXRX-cast01','EXRX-cast02','location','northeast')

xlabel('CH_4 (nmol/kg)');
ylabel('Depth (m)');
set(gca,'xlim',xlc);
set(gca,'ylim',yl);
xticks(xtc);
xticklabels(xtc);
yticks(yt);
yticklabels(yt);

axis ij;

% ------- OCT CH4
subplot(sp(2));
hold on; box on;
set(gca,'fontsize',fs);
set(gca,'tickdir','out');
set(gca,'linewidth',1);
title('October 2023');
%s1.Position = [0.06 0.15 0.19 0.8];
stn = 'EXRX-cast01';

A = find(LISO.Station==stn);
[~,B] = sort(LISO.Depth(A));
CH4eq = CH4sol(LISO.S(A(B)),LISO.T(A(B)),1920e-9).*1000;
%plot(CH4eq,LISO.Depth(A(B)),'--k','linewidth',1.5);

errorbar(LISO.mean_CH4_nM(A(B)),LISO.Depth(A(B)),LISO.std_CH4_nM(A(B)),'horizontal','o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXRX-cast02';
A = find(LISO.Station==stn);
[~,B] = sort(LISO.Depth(A));
errorbar(LISO.mean_CH4_nM(A(B)),LISO.Depth(A(B)),LISO.std_CH4_nM(A(B)),'horizontal','s-r','linewidth',1.5, 'markerfacecolor','r');
%legend('EXRX-cast01','EXRX-cast02','location','northeast')

xlabel('CH_4 (nmol/kg)');
%ylabel('Depth (m)');
set(gca,'xlim',xlc);
set(gca,'ylim',yl);
xticks(xtc);
xticklabels(xtc);
yticks(yt);
yticklabels([]);
axis ij;

% ------- MAY CH4
subplot(sp(3));
hold on; box on;
set(gca,'fontsize',fs);
set(gca,'tickdir','out');
set(gca,'linewidth',1)
%s1.Position = [0.06 0.15 0.19 0.8];
title('May 2024');
stn = 'EXRX-cast01';

A = find(LISM.Station==stn);
[~,B] = sort(LISM.Depth(A));
CH4eq = CH4sol(LISM.S(A(B)),LISM.T(A(B)),1920e-9).*1000;
%plot(CH4eq,LISM.Depth(A(B)),'--k','linewidth',1.5);

errorbar(LISM.mean_CH4_nM(A(B)),LISM.Depth(A(B)),LISM.std_CH4_nM(A(B)),'horizontal','o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXRX-cast02';
A = find(LISM.Station==stn);
[~,B] = sort(LISM.Depth(A));
errorbar(LISM.mean_CH4_nM(A(B)),LISM.Depth(A(B)),LISM.std_CH4_nM(A(B)),'horizontal','s-r','linewidth',1.5, 'markerfacecolor','r');
%legend('EXRX-cast01','EXRX-cast02','location','northeast')

xlabel('CH_4 (nmol/kg)');
%ylabel('Depth (m)');
set(gca,'xlim',xlc);
set(gca,'ylim',yl);
xticks(xtc);
xticklabels(xtc);
yticks(yt);
yticklabels([]);
axis ij;
subplot(sp(1));
hold on; box on;
set(gca,'fontsize',fs);
set(gca,'tickdir','out');
set(gca,'linewidth',1)
%s1.Position = [0.06 0.15 0.19 0.8];
stn = 'EXRX-cast01';

A = find(LISA.Station==stn);
[~,B] = sort(LISA.Depth(A));
CH4eq = CH4sol(LISA.S(A(B)),LISA.T(A(B)),1920e-9).*1000;
%plot(CH4eq,LISA.Depth(A(B)),'--k','linewidth',1.5);

errorbar(LISA.mean_CH4_nM(A(B)),LISA.Depth(A(B)),LISA.std_CH4_nM(A(B)),'horizontal','o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXRX-cast02';
A = find(LISA.Station==stn);
[~,B] = sort(LISA.Depth(A));
errorbar(LISA.mean_CH4_nM(A(B)),LISA.Depth(A(B)),LISA.std_CH4_nM(A(B)),'horizontal','s-r','linewidth',1.5, 'markerfacecolor','r');
%legend('EXRX-cast01','EXRX-cast02','location','northeast')

xlabel('CH_4 (nmol/kg)');
ylabel('Depth (m)');
set(gca,'xlim',xlc);
set(gca,'ylim',yl);
xticks(xtc);
xticklabels(xtc);
axis ij;

% ------- AUGUST N2O
subplot(sp(4));
hold on; box on;
set(gca,'fontsize',fs);
set(gca,'tickdir','out');
set(gca,'linewidth',1)
%s2.Position = [0.28 0.15 0.19 0.8];
stn = 'EXRX-cast01';
A = find(LISA.Station==stn);
[~,B] = sort(LISA.Depth(A));

N2Oeq = N2Osol(LISA.S(A(B)),LISA.T(A(B)),330e-9).*1000;

%plot(N2Oeq,LISA.Depth(A(B)),'--k','linewidth',1.5);
errorbar(LISA.mean_N2O_nM(A(B)),LISA.Depth(A(B)),LISA.std_N2O_nM(A(B)),'horizontal','o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXRX-cast02';
A = find(LISA.Station==stn);
[~,B] = sort(LISA.Depth(A));
errorbar(LISA.mean_N2O_nM(A(B)),LISA.Depth(A(B)),LISA.std_N2O_nM(A(B)),'horizontal','s-r','linewidth',1.5, 'markerfacecolor','r');
%legend('EXRX-cast01','EXRX-cast02','location','northeast')

xlabel('N_2O (nmol/kg)');
xlim(xln);
set(gca,'ylim',yl);
xticks(xtn);
xticklabels(xtn);
yticks(yt);
yticklabels(yt);

%set(gca,'yticklabel',{[]});
%set(gca,'ylim',yl);
%ylabel('Depth (m)');
axis ij;

% ------- OCTOBER N2O
subplot(sp(5));
hold on; box on;
set(gca,'fontsize',fs);
set(gca,'tickdir','out');
set(gca,'linewidth',1)
%s2.Position = [0.28 0.15 0.19 0.8];
stn = 'EXRX-cast01';
A = find(LISO.Station==stn);
[~,B] = sort(LISO.Depth(A));

N2Oeq = N2Osol(LISO.S(A(B)),LISO.T(A(B)),330e-9).*1000;

%plot(N2Oeq,LISO.Depth(A(B)),'--k','linewidth',1.5);
errorbar(LISO.mean_N2O_nM(A(B)),LISO.Depth(A(B)),LISO.std_N2O_nM(A(B)),'horizontal','o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXRX-cast02';
A = find(LISO.Station==stn);
[~,B] = sort(LISO.Depth(A));
errorbar(LISO.mean_N2O_nM(A(B)),LISO.Depth(A(B)),LISO.std_N2O_nM(A(B)),'horizontal','s-r','linewidth',1.5, 'markerfacecolor','r');
%legend('EXRX-cast01','EXRX-cast02','location','northeast')

xlabel('N_2O (nmol/kg)');
xlim(xln);
set(gca,'ylim',yl);
xticks(xtn);
xticklabels(xtn);
yticks(yt);
yticklabels([]);

%set(gca,'yticklabel',{[]});
%set(gca,'ylim',yl);
%ylabel('Depth (m)');
axis ij;

% ------- MAY N2O
subplot(sp(6));
hold on; box on;
set(gca,'fontsize',fs);
set(gca,'tickdir','out');
set(gca,'linewidth',1)
%s2.Position = [0.28 0.15 0.19 0.8];
stn = 'EXRX-cast01';
A = find(LISM.Station==stn);
[~,B] = sort(LISM.Depth(A));

N2Oeq = N2Osol(LISM.S(A(B)),LISM.T(A(B)),330e-9).*1000;

%plot(N2Oeq,LISM.Depth(A(B)),'--k','linewidth',1.5);
errorbar(LISM.mean_N2O_nM(A(B)),LISM.Depth(A(B)),LISM.std_N2O_nM(A(B)),'horizontal','o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXRX-cast02';
A = find(LISM.Station==stn);
[~,B] = sort(LISM.Depth(A));
errorbar(LISM.mean_N2O_nM(A(B)),LISM.Depth(A(B)),LISM.std_N2O_nM(A(B)),'horizontal','s-r','linewidth',1.5, 'markerfacecolor','r');
%legend('EXRX-cast01','EXRX-cast02','location','northeast')

xlabel('N_2O (nmol/kg)');
xlim(xln);
set(gca,'ylim',yl);
xticks(xtn);
xticklabels(xtn);
yticks(yt);
yticklabels([]);

%set(gca,'yticklabel',{[]});
%set(gca,'ylim',yl);
%ylabel('Depth (m)');
axis ij;

% ------- AUGUST OXYGEN
subplot(sp(7));
hold on; box on;
set(gca,'fontsize',fs);
set(gca,'tickdir','out');
set(gca,'linewidth',1)
%s3.Position = [0.5 0.15 0.19 0.8];
stn = 'EXRX-cast01';
A = find(LISA.Station==stn);
[~,B] = sort(LISA.Depth(A));
%plot(O2sol(LISA.S(A(B)),LISA.T(A)),LISA.Depth(A(B)),'--k','linewidth',1.5);
plot(LISA.O2_umolkg(A(B)),LISA.Depth(A(B)),'o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXRX-cast02';
A = find(LISA.Station==stn);
[~,B] = sort(LISA.Depth(A));
plot(LISA.O2_umolkg(A(B)),LISA.Depth(A(B)),'s-r','linewidth',1.5, 'markerfacecolor','r');
%legend('EXRX-cast01','EXRX-cast02','location','northeast')

xlabel('O_2 (\mumol/kg)');
xlim(xlo);
xticks(xto);
xticklabels(xto);
%ylabel('Depth (m)');
yticks(yt);
yticklabels(yt);
axis ij;

% ------- OCTOBER OXYGEN
subplot(sp(8));
hold on; box on;
set(gca,'fontsize',fs);
set(gca,'tickdir','out');
set(gca,'linewidth',1)
%s3.Position = [0.5 0.15 0.19 0.8];
stn = 'EXRX-cast01';
A = find(LISO.Station==stn);
[~,B] = sort(LISO.Depth(A));
%plot(O2sol(LISO.S(A(B)),LISO.T(A)),LISO.Depth(A(B)),'--k','linewidth',1.5);
plot(LISO.O2_umolkg(A(B)),LISO.Depth(A(B)),'o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXRX-cast02';
A = find(LISO.Station==stn);
[~,B] = sort(LISO.Depth(A));
plot(LISO.O2_umolkg(A(B)),LISO.Depth(A(B)),'s-r','linewidth',1.5, 'markerfacecolor','r');
%legend('EXRX-cast01','EXRX-cast02','location','northeast')

xlabel('O_2 (\mumol/kg)');
xlim(xlo);
xticks(xto);
xticklabels(xto);
%ylabel('Depth (m)');
set(gca,'yticklabel',{[]});
set(gca,'ylim',yl);
axis ij;

% ------- MAY OXYGEN
subplot(sp(9));
hold on; box on;
set(gca,'fontsize',fs);
set(gca,'tickdir','out');
set(gca,'linewidth',1)
%s3.Position = [0.5 0.15 0.19 0.8];
stn = 'EXRX-cast01';
A = find(LISM.Station==stn);
[~,B] = sort(LISM.Depth(A));
%plot(O2sol(LISM.S(A(B)),LISM.T(A)),LISM.Depth(A(B)),'--k','linewidth',1.5);
plot(LISM.O2_umolkg(A(B)),LISM.Depth(A(B)),'o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXRX-cast02';
A = find(LISM.Station==stn);
[~,B] = sort(LISM.Depth(A));
plot(LISM.O2_umolkg(A(B)),LISM.Depth(A(B)),'s-r','linewidth',1.5, 'markerfacecolor','r');
%legend('EXRX-cast01','EXRX-cast02','location','northeast')

xlabel('O_2 (\mumol/kg)');
xlim(xlo);
xticks(xto);
xticklabels(xto);
%ylabel('Depth (m)');
set(gca,'yticklabel',{[]});
set(gca,'ylim',yl);
axis ij;

% ------- AUGUST DENSITY
subplot(sp(10));
hold on; box on;
set(gca,'fontsize',fs);
set(gca,'tickdir','out');
set(gca,'linewidth',1)
%s4.Position = [0.72 0.15 0.19 0.8];
stn = 'EXRX-cast01';
A = find(LISA.Station==stn);
[~,B] = sort(LISA.Depth(A));
plot(LISA.PDen(A(B)),LISA.Depth(A(B)),'o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXRX-cast02';
A = find(LISA.Station==stn);
[~,B] = sort(LISA.Depth(A));
plot(LISA.PDen(A(B)),LISA.Depth(A(B)),'s-r','linewidth',1.5, 'markerfacecolor','r');
%legend('EXRX-cast01','EXRX-cast02','location','northeast')

xlabel('\sigma_{\theta} (kg/m^3)');
ylabel('Depth (m)');
xlim(xld);
xticks(xtd);
yticks(yt);
yticklabels(yt);
axis ij;

% ------- OCTOBER DENSITY
subplot(sp(11));
hold on; box on;
set(gca,'fontsize',fs);
set(gca,'tickdir','out');
set(gca,'linewidth',1)
%s4.Position = [0.72 0.15 0.19 0.8];
stn = 'EXRX-cast01';
A = find(LISO.Station==stn);
[~,B] = sort(LISO.Depth(A));
plot(LISO.PDen(A(B)),LISO.Depth(A(B)),'o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXRX-cast02';
A = find(LISO.Station==stn);
[~,B] = sort(LISO.Depth(A));
plot(LISO.PDen(A(B)),LISO.Depth(A(B)),'s-r','linewidth',1.5, 'markerfacecolor','r');
%legend('EXRX-cast01','EXRX-cast02','location','northeast')

xlabel('\sigma_{\theta} (kg/m^3)');
%ylabel('Depth (m)');
xlim(xld);
xticks(xtd);
yticks(yt);
yticklabels([]);
axis ij;

% ------- MAY DENSITY
subplot(sp(12));
hold on; box on;
set(gca,'fontsize',fs);
set(gca,'tickdir','out');
set(gca,'linewidth',1)
%s4.Position = [0.72 0.15 0.19 0.8];
stn = 'EXRX-cast01';
A = find(LISM.Station==stn);
[~,B] = sort(LISM.Depth(A));
plot(LISM.PDen(A(B)),LISM.Depth(A(B)),'o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXRX-cast02';
A = find(LISM.Station==stn);
[~,B] = sort(LISM.Depth(A));
plot(LISM.PDen(A(B)),LISM.Depth(A(B)),'s-r','linewidth',1.5, 'markerfacecolor','r');
%legend('EXRX-cast01','EXRX-cast02','location','northeast')

xlabel('\sigma_{\theta} (kg/m^3)');
%ylabel('Depth (m)');
xlim(xld);
xticks(xtd);
yticks(yt);
yticklabels([]);
axis ij;


wysiwyg;

print -dpng -r300 LIS_EXRX_AOM.png;

%%

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
s1.Position = [0.06 0.15 0.19 0.8];
stn = 'EXRX-cast01';

A = find(LISO.Station==stn);
[~,B] = sort(LISO.Depth(A));
errorbar(LISO.mean_CH4_nM(A(B)),LISO.Depth(A(B)),LISO.std_CH4_nM(A(B)),'horizontal','o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXRX-cast02';
A = find(LISO.Station==stn);
[~,B] = sort(LISO.Depth(A));
CH4eq = CH4sol(LISO.S(A(B)),LISO.T(A(B)),1920e-9).*1000;
%plot(CH4eq,LISO.Depth(A(B)),'--k','linewidth',1.5);

errorbar(LISO.mean_CH4_nM(A(B)),LISO.Depth(A(B)),LISO.std_CH4_nM(A(B)),'horizontal','s-r','linewidth',1.5, 'markerfacecolor','r');
%legend('EXRX-cast01','EXRX-cast02','location','northeast')

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
[~,B] = sort(LISO.Depth(A));
errorbar(LISO.mean_N2O_nM(A(B)),LISO.Depth(A(B)),LISO.std_N2O_nM(A(B)),'horizontal','o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXRX-cast02';
A = find(LISO.Station==stn);
[~,B] = sort(LISO.Depth(A));
N2Oeq = N2Osol(LISO.S(A(B)),LISO.T(A(B)),330e-9).*1000;
plot(N2Oeq,LISO.Depth(A(B)),'--k','linewidth',1.5);


errorbar(LISO.mean_N2O_nM(A(B)),LISO.Depth(A(B)),LISO.std_N2O_nM(A(B)),'horizontal','s-r','linewidth',1.5, 'markerfacecolor','r');
%legend('EXRX-cast01','EXRX-cast02','location','northeast')

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
[~,B] = sort(LISO.Depth(A));
plot(LISO.O2_umolkg(A(B)),LISO.Depth(A(B)),'o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXRX-cast02';
A = find(LISO.Station==stn);
[~,B] = sort(LISO.Depth(A));
plot(O2sol(LISO.S(A(B)),LISO.T(A(B))),LISO.Depth(A(B)),'--k','linewidth',1.5);
plot(LISO.O2_umolkg(A(B)),LISO.Depth(A(B)),'s-r','linewidth',1.5, 'markerfacecolor','r');
%legend('EXRX-cast01','EXRX-cast02','location','northeast')

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
[~,B] = sort(LISO.Depth(A));
plot(LISO.T(A(B)),LISO.Depth(A(B)),'o-b','linewidth',1.5, 'markerfacecolor','b');

stn = 'EXRX-cast02';
A = find(LISO.Station==stn);
[~,B] = sort(LISO.Depth(A));
plot(LISO.T(A(B)),LISO.Depth(A(B)),'s-r','linewidth',1.5, 'markerfacecolor','r');
%legend('EXRX-cast01','EXRX-cast02','location','northeast')

xlabel('Temp. (^oC)');
%ylabel('Depth (m)');
xlim(xlt);
set(gca,'yticklabel',{[]});
set(gca,'ylim',yl);
axis ij;

wysiwyg;

%print -dpng -r300 LISOct23_EXRX_short.png;

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
%legend('EXRX-cast01','EXRX-cast02','location','northeast')

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
%legend('EXRX-cast01','EXRX-cast02','location','northeast')

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
%legend('EXRX-cast01','EXRX-cast02','location','northeast')

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
%legend('EXRX-cast01','EXRX-cast02','location','northeast')

xlabel('Temp. (^oC)');
%ylabel('Depth (m)');
xlim(xlt);
set(gca,'ylim',yl);
axis ij;

wysiwyg;

%print -dpng -r300 LISOct23_EXRX.png;


