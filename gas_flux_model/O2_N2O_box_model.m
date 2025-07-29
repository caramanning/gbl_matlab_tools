% model assumptions:
% start with a given N2O and O2 concentration based on an expected nitrification slope, calculate AOU vs DN2O, upwell to surface
% allow to partially equilibrate, calculate change in AOU to N2O ratio
% sequester again, allow more nitrification, calculate AOU to N2O

N2O_AOU_ratio = 0.02; % nmol N2O per umol O2
AOUi_umolkg = 50; % initial AOU in umol/kg
excess_DN2Oi_nmolkg = 0; % excess DeltaN2O beyond amount expected from AOU
DN2Oi_nmolkg = N2O_AOU_ratio .* AOUi_umolkg + excess_DN2Oi_nmolkg; % initial DeltaN2O in nmol/kg

fN2O = 335e-9; % N2O partial pressure in atm
dtd = 1/24; % time step in days;
dts = dtd .* 86400; % time step in seconds;

S=32;
T=14;
P = 0; % pressure
u10 = 8; % wind speed in m/s
mld = 20; % mixed layer depth in m
t_surf = 20; % time at surface in days
[~, ScN2O] = gasmoldiff(S,T,'N2O'); % Schmidt number for N2O
[~, ScO2] = gasmoldiff(S,T,'O2'); % Schmidt number for O2

kN2O = kgas(u10,ScN2O,'W14'); % N2O gas transfer velocity in m/s
kO2 = kgas(u10,ScO2,'W14');

%% 
mol_per_nmol = 1e-9;
mol_per_umol = 1e-6;

% convert all gas concentrations to mol/m3
dens = sw_dens(S,T,P);

N2Oeq = N2Osol(S,T,fN2O).* dens .* mol_per_umol;
O2eq = O2sol(S,T) .* dens .* mol_per_umol;
%%
AOUi = AOUi_umolkg .* dens .* mol_per_umol;
DN2Oi = DN2Oi_nmolkg .* dens .* mol_per_nmol;

N2Oi = N2Oeq + DN2Oi; % in nmol/kg
O2i = O2eq - AOUi; % in umol/kg

nts = t_surf./dtd; % number of time steps in model

N2Om = nan(nts,1);
O2m = nan(nts,1);

FN2O = nan(nts,1);
FO2 = nan(nts,1);

N2Om(1) = N2Oi;
O2m(1) = O2i;

tm = (0:nts-1).*dtd;

for i = 1:nts-1
    FN2O(i) = kN2O.*(N2Oeq - N2Om(i)) .* dts; % positive flux into ocean, negative flux out of ocean
    N2Om(i+1) = N2Om(i) + FN2O(i)./mld;

    FO2(i) = kO2.*(O2eq - O2m(i)) .* dts;
    O2m(i+1) = O2m(i) + FO2(i)./mld;
end;

DN2Om = N2Om - N2Oeq;
AOUm = O2eq - O2m;
N2O_AOU_ratiom = DN2Om./AOUm.*1000; % ratio in nmol/umol


figure(1)
clf; 
subplot(2,3,1)
hold on;
plot(tm,N2Oeq.*ones(numel(tm),1),'-.k')
plot(tm,N2Om,'-b');
legend('equilibrium','model','location','northeast');
ylabel('N_2O [mol/m^3]');
xlabel('time (d)')

subplot(2,3,2)
hold on;
plot(tm,O2eq.*ones(numel(tm),1),'-.k')
plot(tm,O2m,'-b');
ylabel('O_2 [mol/m^3]');

subplot(2,3,3)
hold on;
plot(tm,N2O_AOU_ratiom)
ylabel('\DeltaN_2O/AOU ratio [nmol/\mumol]');

subplot(2,3,4)
hold on;
plot(tm,DN2Om);
ylabel('\DeltaN2O mol/m^3');

subplot(2,3,5)
hold on;
plot(tm,AOUm);
ylabel('AOU mol/m^3');

subplot(2,3,6)
scatter(AOUm,DN2Om,5,tm, 'o');
c=colorbar;
%clabel(c,'time [d]');
xlabel('AOU mol/m3')
ylabel('\DeltaN_2O mol/m^3')

