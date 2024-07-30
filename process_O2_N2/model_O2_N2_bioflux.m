% calculate N2 flux
load O2_N2_Aug5_Aug12.mat 


S = S_i; T = T_i; 
u10 = 3.*ones(size(S));
dt = 60; % 60 seconds per time step

N2_eq = N2sol(S,T);
N2_meas = N2_i;
mol_per_umol = 1e-6;

N2_eq_molm3 = N2_eq.*sw_dens(S,T,0).* mol_per_umol;
N2_meas_molm3 = N2_meas.*sw_dens(S,T,0).* mol_per_umol;

[~,Sc_N2] = gasmoldiff(S,T,'N2'); % unitless
kgas_N2 = kgas(u10,Sc_N2,'W14'); % m/s
Fdiff_N2 = kgas_N2.*(N2_eq_molm3 - N2_meas_molm3); %mol/ m2 s

diffN2 = diff(N2_meas_molm3); % change in N2 between each time step
dN2dt = [diffN2./dt, 0]; % change in N2 in mol/m3 s

N2_meas_molm3_60 = movmean(N2_meas_molm3,60); %60 min running mean
diffN2_60 = diff(N2_meas_molm3_60); % change in N2 between each time step
dN2dt_60 = [diffN2_60./dt, 0]; % change in N2 in mol/m3 s


O2_eq = O2sol(S,T);
O2_meas = O2_i;

mol_per_umol = 1e-6;

O2_eq_molm3 = O2_eq.*sw_dens(S,T,0).* mol_per_umol;
O2_meas_molm3 = O2_meas.*sw_dens(S,T,0).* mol_per_umol;

[~,Sc_O2] = gasmoldiff(S,T,'O2'); % unitless
kgas_O2 = kgas(u10,Sc_O2,'W14'); % m/s
Fdiff_O2 = kgas_O2.*(O2_eq_molm3 - O2_meas_molm3); %mol/ m2 s

O2_meas_molm3_60 = movmean(O2_meas_molm3,60); %60 min running mean
diffO2 = diff(O2_meas_molm3); % change in O2 between each time step
dO2dt = [diffO2./dt, 0]; % change in O2 in mol/m3 s

diffO2_60 = diff(O2_meas_molm3_60);
dO2dt_60 = [diffO2_60./dt, 0]; % change in O2 in mol/m3 s


%%
Fbio_N2 = nan.*dN2dt;
Fbio_N2(1) = 0;
for i = 1:numel(N2_eq)
   Fbio_N2(i) = dN2dt(i) - Fdiff_N2(i);
end

Fbio_N2_60 = nan.*dN2dt_60;
Fbio_N2_60(1) = 0;
for i = 1:numel(N2_eq)
   Fbio_N2_60(i) = dN2dt_60(i) - Fdiff_N2(i);
end

Fbio_N2_mmol_d = Fbio_N2.*1000.*60*60*24;
Fdiff_N2_mmol_d = Fdiff_N2.*1000.*60*60*24;
Fbio_N2_mmol_d_movmean = movmean(Fbio_N2_mmol_d,60);

Fbio_N2_60_mmol_d = Fbio_N2_60.*1000.*60*60*24;

figure(1)
clf; hold on;
subplot(3,1,1)
hold on; box on;
plot(xti,N2_meas_molm3);
plot(xti,N2_eq_molm3);
plot(xti,N2_meas_molm3_60,'k');
ylabel('N_2 (\mumol/kg)');
legend('meas','eq','meas 60 min smooth');

subplot(3,1,2);
hold on; box on;
%plot(xti,Fbio_N2_mmol_d);
plot(xti,Fdiff_N2_mmol_d);
plot(xti,Fbio_N2_mmol_d_movmean);
plot(xti,Fbio_N2_60_mmol_d);
%legend('Fbio','Fdiff');
ylabel('N_2 flux (mmol m^{-2} d^{-1})');

subplot(3,1,3);
hold on; box on;
plot(xti,T);

%%
% calculate O2 bio flux using non-smoothed (once per minute) data
Fbio_O2 = nan.*dO2dt;
Fbio_O2(1) = 0;
for i = 1:numel(O2_eq)
   Fbio_O2(i) = dO2dt(i) - Fdiff_O2(i);
end

Fbio_O2_mmol_d = Fbio_O2.*1000.*60*60*24;
Fdiff_O2_mmol_d = Fdiff_O2.*1000.*60*60*24;

% calculate O2 bio flux from moving mean
Fbio_O2_mmol_d_movmean = movmean(Fbio_O2_mmol_d,60);

% calculate O2 bio flux using smoothed dO2dt data
for i = 1:numel(O2_eq)
   Fbio_O2_60(i) = dO2dt_60(i) - Fdiff_O2(i);
end

Fbio_O2_60_mmol_d = Fbio_O2_60.*1000.*60*60*24;
%Fdiff_O2_mmol_d = Fdiff_O2.*1000.*60*60*24;
%Fbio_O2_mmol_d_movmean = movmean(Fbio_O2_mmol_d,60);


figure(2)
clf; hold on;
subplot(3,1,1)
hold on; box on;
plot(xti,O2_meas_molm3);
plot(xti,O2_eq_molm3);
plot(xti,O2_meas_molm3_60,'k');
ylabel('O_2 (mol m^{-3})');
legend('meas','eq');

subplot(3,1,2);
hold on; box on;
plot(xti,Fbio_O2_mmol_d);
plot(xti,Fdiff_O2_mmol_d);
plot(xti,Fbio_O2_mmol_d_movmean);
plot(xti,Fbio_O2_60_mmol_d);
%plot(xti_1hr,Fbio_O2_mmol_d_interp);
legend('Fbio','Fdiff','Fbio movmean','Fbio 60 min smooth');
ylabel('O_2 flux (mmol m^{-2} d^{-1}');
set(gca,'ylim',[-5e3 5e3]);

subplot(3,1,3);
hold on; box on;
plot(xti,T);

%% calculating daily fluxes over a 24 hr period
a = find(xti>=datetime(2022,8,6,0,0,0)&xti<datetime(2022,8,7,0,0,0));
[mean(Fbio_O2_60_mmol_d(a)) median(Fbio_O2_60_mmol_d(a))]

b = find(xti>=datetime(2022,8,7,0,0,0)&xti<datetime(2022,8,8,0,0,0));
[mean(Fbio_O2_60_mmol_d(b)) median(Fbio_O2_60_mmol_d(b))]

c = find(xti>=datetime(2022,8,8,0,0,0)&xti<datetime(2022,8,9,0,0,0));
[mean(Fbio_O2_60_mmol_d(c)) median(Fbio_O2_60_mmol_d(c))]

d = find(xti>=datetime(2022,8,9,0,0,0)&xti<datetime(2022,8,10,0,0,0));
[mean(Fbio_O2_60_mmol_d(d)) median(Fbio_O2_60_mmol_d(d))]

e = find(xti>=datetime(2022,8,10,0,0,0)&xti<datetime(2022,8,11,0,0,0));
[mean(Fbio_O2_60_mmol_d(e)) median(Fbio_O2_60_mmol_d(e))]

f = find(xti>=datetime(2022,8,11,0,0,0)&xti<datetime(2022,8,12,0,0,0));
[mean(Fbio_O2_60_mmol_d(f)) median(Fbio_O2_60_mmol_d(f))]

g = find(xti>=datetime(2022,8,6,0,0,0)&xti<datetime(2022,8,12,0,0,0));
[mean(Fbio_O2_60_mmol_d(g)) median(Fbio_O2_60_mmol_d(g))]


%%
a = find(xti>=datetime(2022,8,6,0,0,0)&xti<datetime(2022,8,7,0,0,0));
[mean(Fbio_O2_mmol_d_movmean(a)) median(Fbio_O2_mmol_d_movmean(a))]

b = find(xti>=datetime(2022,8,7,0,0,0)&xti<datetime(2022,8,8,0,0,0));
[mean(Fbio_O2_mmol_d_movmean(b)) median(Fbio_O2_mmol_d_movmean(b))]

c = find(xti>=datetime(2022,8,8,0,0,0)&xti<datetime(2022,8,9,0,0,0));
[mean(Fbio_O2_mmol_d_movmean(c)) median(Fbio_O2_mmol_d_movmean(c))]

d = find(xti>=datetime(2022,8,9,0,0,0)&xti<datetime(2022,8,10,0,0,0));
[mean(Fbio_O2_mmol_d_movmean(d)) median(Fbio_O2_mmol_d_movmean(d))]

e = find(xti>=datetime(2022,8,10,0,0,0)&xti<datetime(2022,8,11,0,0,0));
[mean(Fbio_O2_mmol_d_movmean(e)) median(Fbio_O2_mmol_d_movmean(e))]

f = find(xti>=datetime(2022,8,11,0,0,0)&xti<datetime(2022,8,12,0,0,0));
[mean(Fbio_O2_mmol_d_movmean(f)) median(Fbio_O2_mmol_d_movmean(f))]

g = find(xti>=datetime(2022,8,6,0,0,0)&xti<datetime(2022,8,12,0,0,0));
[mean(Fbio_O2_mmol_d_movmean(g)) median(Fbio_O2_mmol_d_movmean(g))]


%%
h = find(xti>=datetime(2022,8,7,0,0,0)&xti<datetime(2022,8,12,0,0,0));
[mean(Fbio_O2_mmol_d_movmean(g)) median(Fbio_O2_mmol_d_movmean(g))]

[mean(Fbio_N2_mmol_d_movmean(g)) median(Fbio_N2_mmol_d_movmean(g))]

