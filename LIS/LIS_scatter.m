load LISAug23_CH4N2O_CTD.mat
LISA = LISAug23_CH4N2O_CTD; % August

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
xlabel('O2'); ylabel('CH4')

figure(4)
clf; hold on;
scatter(LISA.Dens,LISA.mean_N2O_nM);
xlabel('O2'); ylabel('N2O');
