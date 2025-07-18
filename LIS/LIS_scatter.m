% figure out appropriate ranges for the plots
theta=[10 25];
s=[22 28];

smin=min(s)-0.01.*min(s);
smax=max(s)+0.01.*max(s);
thetamin=min(theta)-0.1*max(theta);
thetamax=max(theta)+0.1*max(theta);
xdim=round((smax-smin)./0.1+1);
ydim=round((thetamax-thetamin)+1);
dens=zeros(ydim,xdim);
thetai=((1:ydim)-1)*1+thetamin;
si=((1:xdim)-1)*0.1+smin;
disp(xdim);disp(ydim);
for j=1:ydim
    for i=1:xdim
        dens(j,i)=sw_dens(si(i),thetai(j),0);
    end
end

dens=dens-1000;

% figure out appropriate ranges for the plots

SAmin=min(s)-0.01.*min(s);
SAmax=max(s)+0.01.*max(s);
thetamin=min(theta)-0.1*max(theta);
thetamax=max(theta)+0.1*max(theta);
xdim=round((smax-smin)./0.1+1);
ydim=round((thetamax-thetamin)+1);
dens_SA_CT=zeros(ydim,xdim);
thetai=((1:ydim)-1)*1+thetamin;
si=((1:xdim)-1)*0.1+smin;
disp(xdim);disp(ydim);
for j=1:ydim
    for i=1:xdim
        dens_SA_CT(j,i)=gsw_rho(si(i),thetai(j),0);
    end
end

dens_SA_CT=dens_SA_CT-1000;


%%
load LISAug23_CH4N2O_CTD.mat
LISA = LISAug23_CH4N2O_CTD; % August

CH4atmdry = 2020e-9;
N2Oatmdry = 338e-9;

UTC_to_local = -4/24;
CH4atmdryA = CH4atmdry;
N2OatmdryA = N2Oatmdry;

load LISAug23_CH4N2O_CTD.mat
LISA = LISAug23_CH4N2O_CTD; % August
LISA.datetime_local = LISA.datetime + UTC_to_local;
LISA.dn_local = datenum(LISA.datetime_local);

LISA.CH4_mean_nmolkg = LISA.mean_CH4_nM./(1000+LISA.PDen).*1000;
LISA.N2O_mean_nmolkg = LISA.mean_N2O_nM./(1000+LISA.PDen).*1000;
LISA.CH4_std_nmolkg = LISA.std_CH4_nM./(1000+LISA.PDen).*1000;
LISA.N2O_std_nmolkg = LISA.std_N2O_nM./(1000+LISA.PDen).*1000;

LISA.N2Oatm_H2Osat = N2OatmdryA .* (1 - vpress(LISA.S,LISA.T));
LISA.N2O_eq_nmolkg = N2Osol(LISA.S,LISA.T,LISA.N2Oatm_H2Osat).*1000;

LISA.CH4atm_H2Osat = CH4atmdryA .* (1 - vpress(LISA.S,LISA.T));
LISA.CH4_eq_nmolkg = CH4sol(LISA.S,LISA.T,LISA.CH4atm_H2Osat)'.*1000;

LISA.DCH4_nmolkg = LISA.CH4_mean_nmolkg - LISA.CH4_eq_nmolkg;
LISA.DN2O_nmolkg = LISA.N2O_mean_nmolkg - LISA.N2O_eq_nmolkg;
LISA.DO2_umolkg = LISA.O2_umolkg - O2sol(LISA.S,LISA.T);

LISA.DCH4 = (LISA.CH4_mean_nmolkg./LISA.CH4_eq_nmolkg - 1).*100;
LISA.DN2O = (LISA.N2O_mean_nmolkg./LISA.N2O_eq_nmolkg - 1).* 100;
LISA.DO2 = (LISA.O2_umolkg./O2sol(LISA.S,LISA.T) - 1).*100;

LISA.SA = gsw_SA_from_SP(LISA.S,LISA.P,LISA.Lon,LISA.Lat);
LISA.CT =  gsw_CT_from_t(LISA.SA,LISA.T,LISA.P);
LISA.PDen_SA_CT = gsw_rho(LISA.SA,LISA.CT,0);

%%

UTC_to_local = -4/24;
CH4atmdryP = CH4atmdry;
N2OatmdryO = N2Oatmdry;

load LISOct23_CH4N2O_CTD.mat
LISO = LISOct23_CH4N2O_CTD; % Oct

LISO = LISOct23_CH4N2O_CTD; % October
LISO.datetime_local = LISO.datetime + UTC_to_local;
LISO.dn_local = datenum(LISO.datetime_local);

LISO.CH4_mean_nmolkg = LISO.mean_CH4_nM./(1000+LISO.PDen).*1000;
LISO.N2O_mean_nmolkg = LISO.mean_N2O_nM./(1000+LISO.PDen).*1000;
LISO.CH4_std_nmolkg = LISO.std_CH4_nM./(1000+LISO.PDen).*1000;
LISO.N2O_std_nmolkg = LISO.std_N2O_nM./(1000+LISO.PDen).*1000;

LISO.N2Oatm_H2Osat = N2OatmdryA .* (1 - vpress(LISO.S,LISO.T));
LISO.N2O_eq_nmolkg = N2Osol(LISO.S,LISO.T,LISO.N2Oatm_H2Osat).*1000;

LISO.CH4atm_H2Osat = CH4atmdryA .* (1 - vpress(LISO.S,LISO.T));
LISO.CH4_eq_nmolkg = CH4sol(LISO.S,LISO.T,LISO.CH4atm_H2Osat)'.*1000;

LISO.DCH4_nmolkg = LISO.CH4_mean_nmolkg - LISO.CH4_eq_nmolkg;
LISO.DN2O_nmolkg = LISO.N2O_mean_nmolkg - LISO.N2O_eq_nmolkg;
LISO.DO2_umolkg = LISO.O2_umolkg - O2sol(LISO.S,LISO.T);

LISO.DCH4 = (LISO.CH4_mean_nmolkg./LISO.CH4_eq_nmolkg - 1).*100;
LISO.DN2O = (LISO.N2O_mean_nmolkg./LISO.N2O_eq_nmolkg - 1).* 100;
LISO.DO2 = (LISO.O2_umolkg./O2sol(LISO.S,LISO.T) - 1).*100;

LISO.SA = gsw_SA_from_SP(LISO.S,LISO.P,LISO.Lon,LISO.Lat);
LISO.CT =  gsw_CT_from_t(LISO.SA,LISO.T,LISO.P);
LISO.PDen_SA_CT = gsw_rho(LISO.SA,LISO.CT,0);

%%

load LISMay24_CH4N2O_CTD.mat
LISM = LISMay24_CH4N2O_CTD; % August

% remove the casts that are in the Eastern Sound as these are not relevant
% to the scatter plot, these are CLIS and ARTG
easternstn = find(LISM.Station=='CLIS-cast01' | LISM.Station=='ARTG-cast01');
LISME = LISM(easternstn,:);
LISM([easternstn],:) = [];


LISM.datetime_local = LISM.datetime + UTC_to_local;
LISM.dn_local = datenum(LISM.datetime_local);

LISM.CH4_mean_nmolkg = LISM.mean_CH4_nM./(1000+LISM.PDen).*1000;
LISM.N2O_mean_nmolkg = LISM.mean_N2O_nM./(1000+LISM.PDen).*1000;
LISM.CH4_std_nmolkg = LISM.std_CH4_nM./(1000+LISM.PDen).*1000;
LISM.N2O_std_nmolkg = LISM.std_N2O_nM./(1000+LISM.PDen).*1000;

LISM.N2Oatm_H2Osat = N2OatmdryA .* (1 - vpress(LISM.S,LISM.T));
LISM.N2O_eq_nmolkg = N2Osol(LISM.S,LISM.T,LISM.N2Oatm_H2Osat).*1000;

LISM.CH4atm_H2Osat = CH4atmdryA .* (1 - vpress(LISM.S,LISM.T));
LISM.CH4_eq_nmolkg = CH4sol(LISM.S,LISM.T,LISM.CH4atm_H2Osat)'.*1000;

LISM.DCH4_nmolkg = LISM.CH4_mean_nmolkg - LISM.CH4_eq_nmolkg;
LISM.DN2O_nmolkg = LISM.N2O_mean_nmolkg - LISM.N2O_eq_nmolkg;
LISM.DO2_umolkg = LISM.O2_umolkg - O2sol(LISM.S,LISM.T);

LISM.DCH4 = (LISM.CH4_mean_nmolkg./LISM.CH4_eq_nmolkg - 1).*100;
LISM.DN2O = (LISM.N2O_mean_nmolkg./LISM.N2O_eq_nmolkg - 1).* 100;
LISM.DO2 = (LISM.O2_umolkg./O2sol(LISM.S,LISM.T) - 1).*100;

LISM.SA = gsw_SA_from_SP(LISM.S,LISM.P,LISM.Lon,LISM.Lat);
LISM.CT =  gsw_CT_from_t(LISM.SA,LISM.T,LISM.P);
LISM.PDen_SA_CT = gsw_rho(LISM.SA,LISM.CT,0);

%%

[mean(LISA.T) std(LISA.T) min(LISA.T) max(LISA.T)]
[mean(LISA.SA) std(LISA.SA) min(LISA.SA) max(LISA.SA)]
%%

[mean(LISO.T) std(LISO.T) min(LISO.T) max(LISO.T)]
[mean(LISO.SA) std(LISO.SA) min(LISO.SA) max(LISO.SA)]

%%
[mean(LISM.T) std(LISM.T) min(LISM.T) max(LISM.T)]
[mean(LISM.SA) std(LISM.SA) min(LISM.SA) max(LISM.SA)]
%%


figure(1)
clf; hold on;
scatter(LISA.O2_umolkg,LISA.mean_CH4_nM);
xlabel('O2'); ylabel('CH4')

figure(2)
clf; hold on;
scatter(LISA.O2_umolkg,LISA.mean_N2O_nM);
xlabel('O2'); ylabel('N2O');

figure(3)
clf; hold on;
scatter(LISA.Dens,LISA.mean_CH4_nM);
xlabel('Dens'); ylabel('CH4')

figure(4)
clf; hold on;
scatter(LISA.Dens,LISA.mean_N2O_nM);
xlabel('Dens'); ylabel('N2O');

%%

Cl = [0 500];
Ol = [-250 150];
Nl = [0 7];

figure(10)
clf; 
subplot(1,2,1)
hold on; box on;
[~,i] = sort(LISA.PDen); 
i = flipud(i);
scatter(LISA.DCH4_nmolkg(i),LISA.DO2_umolkg(i),20,LISA.PDen(i),'filled','o');

[~,i] = sort(LISO.PDen); 
i = flipud(i);
scatter(LISO.DCH4_nmolkg(i),LISO.DO2_umolkg(i),20,LISO.PDen(i),'filled','d');

[~,i] = sort(LISM.PDen); 
i = flipud(i);
scatter(LISM.DCH4_nmolkg(i),LISM.DO2_umolkg(i),20,LISM.PDen(i),'filled','^');
colorbar;
title('CH_4');
% set(gca,'xlim',xl);
% set(gca, 'ylim',yl);
xlim(xl);
ylim(Ol);
xlabel('\DeltaCH_4'); ylabel('\DeltaO_2');

subplot(1,2,2)
hold on; box on;
[~,i] = sort(LISA.PDen); 
i = flipud(i);
scatter(LISA.DN2O_nmolkg(i),LISA.DO2_umolkg(i),20,LISA.PDen(i),'filled','o');

[~,i] = sort(LISO.PDen); 
i = flipud(i);
scatter(LISO.DN2O_nmolkg(i),LISO.DO2_umolkg(i),20,LISO.PDen(i),'filled','d');

[~,i] = sort(LISM.PDen); 
i = flipud(i);
scatter(LISM.DN2O_nmolkg(i),LISM.DO2_umolkg(i),20,LISM.PDen(i),'filled','^');
colorbar;
title('N_2O');
%xlim(Sl);
%ylim(Tl);
xlabel('\DeltaN_2O'); ylabel('\DeltaO_2');

%%
Cl = [0 500];
Ol = [-250 150];
Nl = [0 7];

nr = 2; nc = 3;

figure(10)
clf;
sp=tight_subplot(nr,nc,[.05 .03],[.15 .1],[.08 .04]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 5.5]);

subplot(sp(1))
hold on; box on;
[~,i] = sort(LISA.PDen); 
i = flipud(i);
scatter(LISA.DO2_umolkg(i),LISA.DCH4_nmolkg(i),20,LISA.PDen(i),'filled','o');

colorbar;
title('Aug');
% set(gca,'xlim',xl);
% set(gca, 'ylim',yl);
xlim(Ol);
ylim(Cl);
ylabel('\DeltaCH_4'); xlabel('\DeltaO_2');

subplot(sp(2))
hold on; box on;
[~,i] = sort(LISO.PDen); 
i = flipud(i);
scatter(LISO.DO2_umolkg(i),LISO.DCH4_nmolkg(i),20,LISO.PDen(i),'filled','d');

colorbar;
title('Oct');
% set(gca,'xlim',xl);
% set(gca, 'ylim',yl);
xlim(Ol);
ylim(Cl);
ylabel('\DeltaCH_4');
xlabel('\DeltaO_2');

subplot(sp(3));
hold on; box on;
[~,i] = sort(LISM.PDen); 
i = flipud(i);
scatter(LISM.DO2_umolkg(i),LISM.DCH4_nmolkg(i),20,LISM.PDen(i),'filled','^');
colorbar;
title('CH_4');
% set(gca,'xlim',xl);
% set(gca, 'ylim',yl);
xlim(Ol);
ylim(Cl);
ylabel('\DeltaCH_4'); 
xlabel('\DeltaO_2');

subplot(sp(4))
hold on; box on;
[~,i] = sort(LISA.PDen); 
i = flipud(i);
scatter(LISA.DO2_umolkg(i),LISA.DN2O_nmolkg(i),20,LISA.PDen(i),'filled','o');

colorbar;
title('CH_4');
% set(gca,'xlim',xl);
% set(gca, 'ylim',yl);
xlim(Ol);
ylim(Nl);
ylabel('\DeltaN_2O'); xlabel('\DeltaO_2');

subplot(sp(5))
hold on; box on;
[~,i] = sort(LISO.PDen); 
i = flipud(i);
scatter(LISO.DO2_umolkg(i),LISO.DN2O_nmolkg(i),20,LISO.PDen(i),'filled','d');

colorbar;
title('CH_4');
% set(gca,'xlim',xl);
% set(gca, 'ylim',yl);
xlim(Ol);
ylim(Nl);
xlabel('\DeltaO_2');

subplot(sp(6))
hold on; box on;
[~,i] = sort(LISM.PDen); 
i = flipud(i);
scatter(LISM.DO2_umolkg(i),LISM.DN2O_nmolkg(i),20,LISM.PDen(i),'filled','^');
colorbar;
title('CH_4');
% set(gca,'xlim',xl);
% set(gca, 'ylim',yl);
xlim(Ol);
ylim(Nl);
xlabel('\DeltaO_2');

wysiwyg;

print -dpng -r300 20250712_LIS_scatter_O2_CH4_N2O.png;

print -depsc -r300 20250707_LIS_scatter_O2_CH4_N2O.eps;


%%

subplot(1,2,2)
hold on; box on;
[~,i] = sort(LISA.PDen); 
i = flipud(i);
scatter(LISA.DN2O_nmolkg(i),LISA.DO2_umolkg(i),20,LISA.PDen(i),'filled','o');

[~,i] = sort(LISO.PDen); 
i = flipud(i);
scatter(LISO.DN2O_nmolkg(i),LISO.DO2_umolkg(i),20,LISO.PDen(i),'filled','d');

[~,i] = sort(LISM.PDen); 
i = flipud(i);
scatter(LISM.DN2O_nmolkg(i),LISM.DO2_umolkg(i),20,LISM.PDen(i),'filled','^');
colorbar;
title('N_2O');
%xlim(Sl);
%ylim(Tl);
xlabel('\DeltaN_2O'); ylabel('\DeltaO_2');


%%
Sl = [22 27.4];
Tl = [10 25];

yt = [10 15 20 25]; % yaxis ticks
xt = [22 23 24 25 26 27]; %x axis ticks
yl = [10 25]; %yaxls limit
xl = [22 27.4]; %x axis limit

mec = [0.5 0.5 0.5]; % marker edge color
mec = 'none';
mfa = 0.5; % marker face alpha (transparency)

cch4 = [25 440];
cn2o = [8.3 16.4];

ms=60;

nr = 1; % number of rows
nc = 3;  % number of columns
lw=1;
fsl = 12;
fig=figure(11)
clf;
sp=tight_subplot(nr,nc,[.05 .03],[.15 .1],[.08 .04]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 5.5]);


subplot(sp(1))
hold on; box on;
set(gca,'linewidth',lw);
set(gca, 'FontName', 'Arial')
[c,h]=contour(si,thetai,dens,'--','color',[0.5 0.5 0.5]);
%[c,h]=contour(si,thetai,dens,'--','color','r');
%clabel(c,h,'LabelSpacing',1000);
set(gca,'fontsize',fsl);

colormap('turbo')
[~,i] = sort(LISA.CH4_mean_nmolkg);
s = scatter(LISA.S(i),LISA.T(i),ms,LISA.CH4_mean_nmolkg(i),'filled','o','markeredgecolor',mec);
s.MarkerFaceAlpha = mfa;
s.MarkerEdgeAlpha = mfa;
[~,i] = sort(LISO.CH4_mean_nmolkg);
s = scatter(LISO.S(i),LISO.T(i),ms,LISO.CH4_mean_nmolkg(i),'filled','d','markeredgecolor',mec);
s.MarkerFaceAlpha = mfa;
s.MarkerEdgeAlpha = mfa;
[~,i] = sort(LISM.CH4_mean_nmolkg);
s = scatter(LISM.S(i),LISM.T(i),ms,LISM.CH4_mean_nmolkg(i),'filled','^','markeredgecolor',mec);
s.MarkerFaceAlpha = mfa;
s.MarkerEdgeAlpha = mfa;
a = colorbar('location','southoutside');
a.TickDirection = 'out';
caxis(cch4);
a.Label.String = 'CH_4 [nmol kg^{-1}]';
a.Label.FontSize = fsl;
%title('CH_4');
% set(gca,'xlim',xl);
% set(gca, 'ylim',yl);
set(gca,'tickdir','out')
xlim(xl);
ylim(yl);
    set(gca,'ytick',yt);
    set(gca,'xtick',xt);
    set(gca,'yticklabel',yt,'fontsize',fsl);
    set(gca,'xticklabel',xt);
    xlabel('Salinity [PSS-78]');
    ylabel('Temperature [^oC]');

subplot(sp(2))
hold on; box on;
set(gca,'linewidth',lw);
set(gca, 'FontName', 'Arial')
[c,h]=contour(si,thetai,dens,'--','color',[0.5 0.5 0.5]);
%clabel(c,h,'LabelSpacing',1000);
set(gca,'fontsize',fsl);

[~,i] = sort(LISA.N2O_mean_nmolkg);
s = scatter(LISA.S(i),LISA.T(i),ms,LISA.N2O_mean_nmolkg(i),'filled','o','markeredgecolor',mec);
s.MarkerFaceAlpha = mfa;
s.MarkerEdgeAlpha = mfa;

[~,i] = sort(LISO.mean_N2O_nM);
s = scatter(LISO.S(i),LISO.T(i),ms,LISO.N2O_mean_nmolkg(i),'filled','d','markeredgecolor',mec);
s.MarkerFaceAlpha = mfa;
s.MarkerEdgeAlpha = mfa;
[~,i] = sort(LISM.mean_N2O_nM);
s = scatter(LISM.S(i),LISM.T(i),ms,LISM.N2O_mean_nmolkg(i),'filled','^','markeredgecolor',mec);
s.MarkerFaceAlpha = mfa;
s.MarkerEdgeAlpha = mfa;
a = colorbar('location','southoutside');
a.TickDirection = 'out';
caxis(cn2o);
a.Label.String = 'N_2O [nmol kg^{-1}]';
a.Label.FontSize = fsl;
%title('N_2O');
set(gca,'tickdir','out')
xlim(xl);
ylim(yl);
    set(gca,'ytick',yt);
    set(gca,'xtick',xt);
    set(gca,'xticklabel',xt);
    xlabel('Salinity [PSS-78]');

subplot(sp(3))
hold on; box on;
set(gca,'linewidth',lw);
set(gca, 'FontName', 'Arial')
[c,h]=contour(si,thetai,dens,'--','color',[0.5 0.5 0.5]);
%clabel(c,h,'LabelSpacing',1000);
set(gca,'fontsize',fsl);
[~,i] = sort(LISA.mean_N2O_nM);
s=scatter(LISA.S(i),LISA.T(i),ms,LISA.O2_umolkg(i),'filled','o','markeredgecolor',mec);
s.MarkerFaceAlpha = mfa;
s.MarkerEdgeAlpha = mfa;
[~,i] = sort(LISO.mean_N2O_nM);
s=scatter(LISO.S(i),LISO.T(i),ms,LISO.O2_umolkg(i),'filled','d','markeredgecolor',mec);
s.MarkerFaceAlpha = mfa;
s.MarkerEdgeAlpha = mfa;
[~,i] = sort(LISM.mean_N2O_nM);
s=scatter(LISM.S(i),LISM.T(i),ms,LISM.O2_umolkg(i),'filled','^','markeredgecolor',mec);
s.MarkerFaceAlpha = mfa;
s.MarkerEdgeAlpha = mfa;
colorbar('location','southoutside');
a = colorbar('location','southoutside');
a.TickDirection = 'out';
a.Label.String = 'O_2 [\mumol kg^{-1}]';
a.Label.FontSize = fsl;
%title('O_2')
% set(gca,'xlim',xl);
% set(gca, 'ylim',yl);
set(gca,'tickdir','out')
xlim(xl);
ylim(yl);
    set(gca,'ytick',yt);
    set(gca,'xtick',xt);
    set(gca,'xticklabel',xt);    
    xlabel('Salinity [PSS-78]');


wysiwyg;

print -dpng -r300 20250709_WLIS_TSplot_noclabel.png;
print(gcf,'-depsc','-vector','20250709_WLIS_TSplot_noclabel.eps');
print(gcf,'-dpdf','-vector','20250709_WLIS_TSplot_noclabel.pdf');
epsclean('20250709_WLIS_TSplot_noclabel.eps','20250709_WLIS_TSplot_noclabel_clean.eps');

%%
figure(2)
clf; hold on;
scatter(LISA.S,LISA.T,20,LISA.mean_N2O_nM,'filled','o');
colorbar;
title('N_2O');
xlabel('S'); ylabel('T');

figure(3)
clf; hold on; box on;
scatter(LISA.S,LISA.T,20,LISA.O2_umolkg,'filled','o');
colorbar;
title(O_2)
xlabel('S'); ylabel('T');

%%
figure(2)
clf; hold on;
scatter(LISA.O2_umolkg,LISA.mean_N2O_nM);
xlabel('O2'); ylabel('N2O');

figure(3)
clf; hold on;
scatter(LISA.Dens,LISA.mean_CH4_nM);
xlabel('Dens'); ylabel('CH4')

figure(4)
clf; hold on;
scatter(LISA.Dens,LISA.mean_N2O_nM);
xlabel('Dens'); ylabel('N2O');

