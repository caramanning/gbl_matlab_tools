%foldername = "G:\Shared drives\Gas Biogeochemistry Lab\projects\CIRCA 2022 seed grant\test deployment August 1 to 3";
%foldername = "H:\Shared drives\Gas Biogeochemistry Lab\projects\CIRCA 2022 seed grant\deployment Aug 5";
foldername = "G:\Shared drives\Gas Biogeochemistry Lab\projects\CIRCA 2022 seed grant\deployment Sep 12";

% deep sensors
fname_O2 =  'PME_cc_383325.mat';
filepath_O2 = fullfile(foldername,fname_O2);
load(filepath_O2);
O2_deep = PME_O2_cc;

fname_CTD =  'S12071_Sep12_Nov22.mat';
filepath_CTD = fullfile(foldername,fname_CTD);
load(filepath_CTD);
CTD_deep = S12071;
CTD_deep.sal = CTD_deep.sal + 0.69; % correct for offset observed between sensors


fname_TDGP = 'miniTDGP 20220912 to 20221122.mat';
filepath_TDGP = fullfile(foldername,fname_TDGP);
load(filepath_TDGP);
TDGP_deep = miniTDGP;

% shallow sensors
fname_O2 =  'PME_cc_364333.mat';
filepath_O2 = fullfile(foldername,fname_O2);
load(filepath_O2);
O2_shall = PME_O2_cc;

fname_CTD =  'S12073_Sep12_Nov22.mat';
filepath_CTD = fullfile(foldername,fname_CTD);
load(filepath_CTD);
CTD_shall = S12073;
CTD_shall.sal = CTD_shall.sal - 0.69; % correct for offset observed between sensors


%%

% convert times from UTC to EST (4 hr offset)
O2_deep.datetime_EST = O2_deep.datetime - 4/24;
O2_shall.datetime_EST = O2_shall.datetime - 4/24;


% calculate datetime from miniTDGP time
TDGP_deep.datetime = datetime(miniTDGP.yyyy,miniTDGP.mm,miniTDGP.dd,miniTDGP.HH,miniTDGP.MM,miniTDGP.SS);

% calculate O2 in umol/L from mg/L
MW_O2 = 31.998;  % molar weight of O2 in g/mol
mg_per_g = 1000; % milligrams per gram
umol_per_mol = 1e6;
O2_shall.DO_umolL = O2_shall.DO ./ MW_O2 ./ mg_per_g .* umol_per_mol;
O2_deep.DO_umolL = O2_deep.DO ./ MW_O2 ./ mg_per_g .* umol_per_mol;


figure(1)
clf; %clear existing plot
hold on;
%xl = [datetime(2022,8,1,16,0,0) datetime(2022,8,3,9,20,0)];
%xl = [datetime(2022,8,5,11,0,0) datetime(2022,8,12,11,0,0)];
%xl = [datetime(2022,9,12,18,0,0) datetime(2022,11,21,12,0,0)];
xl = [datetime(2022,9,12,18,0,0) datetime(2022,9,18,18,0,0)];

% oxygen plot
%subplot called as (number of rows),(number of columns), (subplot number)
subplot(5,1,1) 
hold on; box on;
plot(O2_deep.datetime_EST,O2_deep.DO_umolL)
plot(O2_shall.datetime_EST,O2_shall.DO_umolL)
%plot([datetime(2022,8,5,11,15,0) datetime(2022,8,12,11,0,0)], [215 215])
%xlim([datetime(2022,8,5,11,15,0) datetime(2022,8,12,11,0,0)])
%ylim([170 240])
%ylim([100 350]);
set(gca,'xticklabel',[]) % remove date label from all but bottom
%xlim([737838 737848])
ylabel("O_2 (\mumol/L)")
legend('deep','shall','location','northeast')
xlim(xl);
title('surface and benthic')


subplot(5,1,2)
hold on; box on;
plot(TDGP_deep.datetime,TDGP_deep.P);
ylabel('Gas pressure (mbar)')
set(gca,'xticklabel',[])
xlim(xl)


% temperature plot
% there is temperature from multiple sensors
% SO = Star-Oddi (CTD sensor)
% PME = Precision Measurement Engineering (O2 sensor)
% PO = Pro-Oceanus (TDGP sensor)
subplot(5,1,3)
hold on; box on;
plot(CTD_deep.sample_date,CTD_deep.temp); % star-oddi temperature
plot(CTD_shall.sample_date,CTD_shall.temp); % star-oddi temperature

plot(O2_deep.datetime_EST,O2_deep.T); % PME temperature
plot(O2_shall.datetime_EST,O2_shall.T); % PME temperature


legend('SO deep','SO shall','PME deep','PME shall','location','northeast')
%ylim([22.5 29])
%xlim([datetime(2022,8,5,11,15,0) datetime(2022,8,12,11,0,0)])
set(gca,'xticklabel',[]) % remove date label from all but bottom
ylabel("Temp (^oC)")
xlim(xl)


% salinity plot
subplot(5,1,4)
hold on; box on;
plot(CTD_deep.sample_date,CTD_deep.sal)
plot(CTD_shall.sample_date,CTD_shall.sal)
legend('deep','shall')

%ylim([20 30]) 
%ylim([24 30]);
%xlim([datetime(2022,8,5,11,15,0) datetime(2022,8,12,11,0,0)])
set(gca,'xticklabel',[]) % remove date label from all but bottom
ylabel('Sal (PSS)')
xlim(xl)



% depth plot
% note: for now I am applying a manual offset of 1 m to the depth data
% based on the depth we measured with the rope. 
% labeled as Approx Depth
% we need to verify the offset more accurately
subplot(5,1,5)
hold on; box on;
plot(CTD_deep.sample_date,CTD_deep.depth+1)
plot(CTD_shall.sample_date,CTD_shall.depth+1)
%xlim([datetime(2022,8,5,11,15,0) datetime(2022,8,12,11,0,0)])
ylabel('Approx. Depth (m)')
xlim(xl)

%ylim([0.85 2.25]);


%print -dpng -r300 'Data Sep 12 to Nov 21 Shallow Deep.png'
print -dpng -r300 'Data Sep 12 to Sep 18 Shallow Deep.png'