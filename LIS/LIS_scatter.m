load LISAug23_CH4N2O_CTD.mat
LISA = LISAug23_CH4N2O_CTD; % August

load LISOct23_CH4N2O_CTD.mat
LISO = LISOct23_CH4N2O_CTD; % Oct

load LISAug23_CH4N2O_CTD.mat
LISM = LISMay24_CH4N2O_CTD; % May


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

Sl = [22 27];
Tl = [8 25];
figure(1)
clf; 
subplot(1,3,1)
hold on; box on;
[~,i] = sort(LISA.mean_CH4_nM);
scatter(LISA.S(i),LISA.T(i),20,LISA.mean_CH4_nM(i),'filled','o');
colorbar;
title('CH_4');
% set(gca,'xlim',xl);
% set(gca, 'ylim',yl);
xlim(Sl);
ylim(Tl);
xlabel('S'); ylabel('T');

subplot(1,3,2)
hold on; box on;
scatter(LISO.S,LISO.T,20,LISO.mean_CH4_nM,'filled','o');
colorbar;
title('CH_4');
xlim(Sl);
ylim(Tl);
xlabel('S'); ylabel('T');

subplot(1,3,3)
hold on; box on;
scatter(LISM.S,LISM.T,20,LISM.mean_CH4_nM,'filled','o');
colorbar;
title('CH_4');
xlim(Sl);
ylim(Tl);
xlabel('S'); ylabel('T');

%%
Sl = [22 27];
Tl = [8 25];
ms=30;
figure(1)
clf; 
hold on; box on;
[~,i] = sort(LISA.mean_CH4_nM);
scatter(LISA.S(i),LISA.T(i),ms,LISA.mean_CH4_nM(i),'filled','o');
[~,i] = sort(LISO.mean_CH4_nM);
scatter(LISO.S(i),LISO.T(i),ms,LISO.mean_CH4_nM(i),'filled','o');
[~,i] = sort(LISM.mean_CH4_nM);
scatter(LISM.S(i),LISM.T(i),ms,LISM.mean_CH4_nM(i),'filled','o');
colorbar;
title('CH_4');
% set(gca,'xlim',xl);
% set(gca, 'ylim',yl);
xlim(Sl);
ylim(Tl);
xlabel('S'); ylabel('T');

Sl = [22 27];
Tl = [8 25];
figure(1)
clf; 
hold on; box on;
[~,i] = sort(LISA.mean_N2O_nM);
scatter(LISA.S(i),LISA.T(i),20,LISA.mean_N2O_nM(i),'filled','o');
[~,i] = sort(LISO.mean_N2O_nM);
scatter(LISO.S(i),LISO.T(i),20,LISO.mean_N2O_nM(i),'filled','o');
[~,i] = sort(LISM.mean_N2O_nM);
scatter(LISM.S(i),LISM.T(i),20,LISM.mean_N2O_nM(i),'filled','o');
colorbar;
title('N_2O');
% set(gca,'xlim',xl);
% set(gca, 'ylim',yl);
xlim(Sl);
ylim(Tl);
xlabel('S'); ylabel('T');

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

