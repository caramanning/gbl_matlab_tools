%foldername = "G:\Shared drives\Gas Biogeochemistry Lab\projects\CIRCA 2022 seed grant\test deployment August 1 to 3";
foldername = "G:\Shared drives\Gas Biogeochemistry Lab\projects\CIRCA 2022 seed grant\deployment Sep 12";

% define the time grid to use for interpolation
dt = 1/60/24; % time step is 1 minute
t_start = datetime(2022,9,12,18,0,0);
t_end = datetime(2022,11,21,12,0,0);
xti = t_start:dt:t_end;

fname_O2 =  'PME_cc_383325.mat';
filepath_O2 = fullfile(foldername,fname_O2);
load(filepath_O2);

fname_CTD =  'S12071_Sep12_Nov22.mat';
filepath_CTD = fullfile(foldername,fname_CTD);
load(filepath_CTD);

fname_TDGP = 'miniTDGP 20220912 to 20221122.mat';
filepath_TDGP = fullfile(foldername,fname_TDGP);
load(filepath_TDGP);
miniTDGP(102494:end,:) = []; % remove values at end from Nov 22

% O2_offset_drift
cal_times = [datetime(2022,9,11,10,16,0) datetime(2022,11,22,5,11,0)]; % times that cal samples were collected
cal_ratio = [1.00 1.00]; % ratio of titrator value/sensor value


% convert times from UTC to EST (4 hr offset)
PME_O2_cc.datetime_EST = PME_O2_cc.datetime - 4/24;

% calculate datetime from miniTDGP time
miniTDGP.datetime = datetime(miniTDGP.yyyy,miniTDGP.mm,miniTDGP.dd,miniTDGP.HH,miniTDGP.MM,miniTDGP.SS);

% calculate O2 in umol/L from mg/L
MW_O2 = 31.998;  % molar weight of O2 in g/mol
mg_per_g = 1000; % milligrams per gram
umol_per_mol = 1e6;
PME_O2_cc.DO_umolL = PME_O2_cc.DO ./ MW_O2 ./ mg_per_g .* umol_per_mol;

% calculate the interpolated S, T, O2, TDGP with time step matched
S_i_raw = interp1(S12071.sample_date,S12071.sal,xti);
T_i = interp1(PME_O2_cc.datetime_EST,PME_O2_cc.T,xti); % PME temperature
T_i_PO = interp1(miniTDGP.datetime,miniTDGP.T,xti);
T_i_SO = interp1(S12071.sample_date,S12071.temp,xti);
%T_i_SO = interp1(S12071.datetime,)
O2_umolL_i = interp1(PME_O2_cc.datetime_EST,PME_O2_cc.DO_umolL,xti); % PME temperature
TDGP_i = interp1(miniTDGP.datetime,miniTDGP.P,xti);
cal_i = interp1(cal_times,cal_ratio,xti); % interpolated offset

% manually correct the salinity values that are too low. 
dt_good = datetime(2022,08,06,13,27,0); % assume all values from here onward are good
mean_S_good = mean(S_i_raw(xti>=datetime(2022,08,06,13,27,0)));
S_i = S_i_raw; % create new variable for the corrected salinity data
S_i(xti<dt_good) = mean_S_good; % set all values at start of time series to mean from end of time series
S_i(S_i_raw<27.5) = mean_S_good; %set erroneous low values on Aug 9 to the mean value too

%create new salinity time series with low values converted to NaN
S_i_lowexcluded = S_i_raw;
S_i_lowexcluded(xti<dt_good) = nan;
S_i_lowexcluded(S_i_raw<27.5) = nan;


d_offset = 1; % approximate depth offset of 1 m
D_approx_i = interp1(S12071.sample_date,S12071.depth+d_offset,xti);

lat = 42.1;
P_approx_i = sw_pres(D_approx_i',lat);

O2_i_uncal = O2_umolL_i ./ sw_dens(S_i,T_i,P_approx_i') .*1000;

O2_i = O2_i_uncal.*cal_i; % convert from uncalibrated to calibrated

% calculate O2 saturation state
O2_eq = O2sol(S_i,T_i);
RO2_i = (O2_i./O2_eq); % ratio of measured to equilibrium concentration
O2sat_i = RO2_i*100;
DO2_i = (RO2_i -1) * 100;

%PME_O2_cc.DO_umolkg = PME_O2_cc.DO_umolL ./ sw_dens(S_i,T_i,D_approx_i);

% calculate N2 saturation state
slp = 1013.25.*ones(size(T_i)); 
air.slp = slp; % for now use 1013.25 mbar as air pressure, need to download real data
[dat,p,air] = n2_from_tp(TDGP_i,RO2_i,[],[],T_i,S_i,air);

DN2_i = (dat.n2sat - 1) * 100;
N2_i = dat.n2_molkg.*1e6; % N2 in umol/kg
N2_eq = N2sol(S_i,T_i);

% manually remove first 4 hours of N2 data as it is suspect
%N2_i(xti<datetime(2022,8,5,14,30,0))  = NaN;
%DN2_i(xti<datetime(2022,8,5,14,30,0))  = NaN;
%dat.n2sat(xti<datetime(2022,8,5,14,30,0))  = NaN;
% O2_gp = RO2_i.*0.21.*1014;

% O2_gp = RO2_i.*0.21.*1014;
% 
% N2_gp = TDGP_i - O2_gp - 0.01*1014;
% 
% RN2 = N2_gp./(1014.*0.78);
% DN2 = (RN2 -1) * 100;
% N2_eq = N2sol(S_i,T_i);
% N2_i = RN2.*N2_eq;

figure(1)
clf;
subplot(6,1,1)
hold on; box on;
plot(xti,O2_i);
plot(xti,O2_eq);
legend('O_2','O_{2,eq}','location','eastoutside');
axis tight;
ylabel('\mumol kg^{-1}');
%ylim([100 350]);
xlim([t_start t_end]);

subplot(6,1,2)
hold on; box on;
plot(xti,N2_i);
plot(xti,N2_eq);
legend('N_2','N_{2,eq}','location','eastoutside');
axis tight;
%set(gca,'xticklabel',{[]})
ylabel('\mumol kg^{-1}');
%ylim([350 460]);
xlim([t_start t_end]);

subplot(6,1,3)
hold on; box on;
plot(xti,DO2_i);
plot(xti,0.*T_i,':k');
legend('\DeltaO_2','location','eastoutside');
ylabel('\DeltaO_2 [%]');
axis tight;
%ylim([-60 60]);
xlim([t_start t_end]);

subplot(6,1,4)
hold on; box on;
plot(xti,DN2_i);
plot(xti,0.*T_i,':k');
legend('\DeltaN_{2}','location','eastoutside');
ylabel('\DeltaN_2 [%]');
axis tight;
%ylim([-8 15]);
xlim([t_start t_end]);


subplot(6,1,5)
hold on; box on;
plot(xti,T_i);
%plot(xti,T_i_PO);
plot(xti,T_i_SO,'k');
legend('PME', 'SO', 'location','eastoutside');
ylabel('temp [^oC]');
axis tight;
%ylim([20 29]);
xlim([t_start t_end]);

subplot(6,1,6)
hold on; box on;
plot(xti,S_i_lowexcluded);
legend('salin','location','eastoutside');
ylabel('sal [PSS]');
axis tight;
%ylim([28 30]);
xlim([t_start t_end]);

print -dpng -r300 O2_N2_sat_Sep12_Nov22.png;

figure(2)
clf; 
hold on; box on;
plot(xti,T_i);
%plot(xti,T_i_PO);
plot(xti,T_i_SO,'k');
legend('PME', StarOddi', 'location','north');
ylabel('temp [^oC]');
axis tight;
xlim([t_start t_end]);

print -dpng -r300 temp_Sep12_Nov22.png;

figure(3)
clf; 
subplot(3,1,1)
hold on; box on;
plot(xti,T_i);
plot(xti,T_i_PO);
plot(xti,T_i_SO,'k');
legend('PME', 'PO', 'SO', 'location','north');
ylabel('temp [^oC]');
axis tight;
xlim([t_start t_end]);

subplot(3,1,2)
hold on; box on;
plot(xti,TDGP_i);
%legend('TDGP', 'location','eastoutside');
ylabel('TDGP [mbar]');
axis tight;
xlim([t_start t_end]);

subplot(3,1,3)
hold on; box on;
plot(xti,O2_i);
%legend('O_2', 'location','eastoutside');
ylabel('O_2 [\mumol kg^{-1}]');
axis tight;
xlim([t_start t_end]);

print -dpng -r300 O2_N2_Sep12_Nov22.png;


