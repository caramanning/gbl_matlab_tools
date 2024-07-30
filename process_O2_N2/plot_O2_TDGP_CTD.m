%foldername = "G:\Shared drives\Gas Biogeochemistry Lab\projects\CIRCA 2022 seed grant\test deployment August 1 to 3";
foldername = "H:\Shared drives\Gas Biogeochemistry Lab\projects\CIRCA 2022 seed grant\deployment Aug 5";

fname_O2 =  'PME_O2_cc.mat';
filepath_O2 = fullfile(foldername,fname_O2);
load(filepath_O2);

fname_CTD =  'S12071.mat';
filepath_CTD = fullfile(foldername,fname_CTD);
load(filepath_CTD);

fname_TDGP = 'miniTDGP_20220805.mat';
filepath_TDGP = fullfile(foldername,fname_TDGP);
load(filepath_TDGP);

%%

% convert times from UTC to EST (4 hr offset)
PME_O2_cc.datetime_EST = PME_O2_cc.datetime - 4/24;

% calculate datetime from miniTDGP time
miniTDGP.datetime = datetime(miniTDGP.yyyy,miniTDGP.mm,miniTDGP.dd,miniTDGP.HH,miniTDGP.MM,miniTDGP.SS);

% calculate O2 in umol/L from mg/L
MW_O2 = 31.998;  % molar weight of O2 in g/mol
mg_per_g = 1000; % milligrams per gram
umol_per_mol = 1e6;
PME_O2_cc.DO_umolL = PME_O2_cc.DO ./ MW_O2 ./ mg_per_g .* umol_per_mol;

figure(1)
clf; %clear existing plot
hold on;
%xl = [datetime(2022,8,1,16,0,0) datetime(2022,8,3,9,20,0)];
xl = [datetime(2022,8,5,11,0,0) datetime(2022,8,12,11,0,0)];

% oxygen plot
%subplot called as (number of rows),(number of columns), (subplot number)
subplot(5,1,1) 
hold on; box on;
plot(PME_O2_cc.datetime_EST,PME_O2_cc.DO_umolL)
%plot([datetime(2022,8,5,11,15,0) datetime(2022,8,12,11,0,0)], [215 215])
%xlim([datetime(2022,8,5,11,15,0) datetime(2022,8,12,11,0,0)])
%ylim([170 240])
ylim([100 350]);
set(gca,'xticklabel',[]) % remove date label from all but bottom
%xlim([737838 737848])
ylabel("O_2 (\mumol/L)")
xlim(xl);


subplot(5,1,2)
hold on; box on;
plot(miniTDGP.datetime,miniTDGP.P);
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
plot(S12071.sample_date,S12071.temp); % star-oddi temperature
plot(PME_O2_cc.datetime_EST,PME_O2_cc.T); % PME temperature
plot(miniTDGP.datetime,miniTDGP.T); % Pro-Oceanus temperature
legend('SO T','PME T','PO T')
%ylim([22.5 29])
%xlim([datetime(2022,8,5,11,15,0) datetime(2022,8,12,11,0,0)])
set(gca,'xticklabel',[]) % remove date label from all but bottom
ylabel("Temp (^oC)")
xlim(xl)


% salinity plot
subplot(5,1,4)
hold on; box on;
plot(S12071.sample_date,S12071.sal)
%ylim([20 30]) 
ylim([24 30]);
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
plot(S12071.sample_date,S12071.depth+1)
%xlim([datetime(2022,8,5,11,15,0) datetime(2022,8,12,11,0,0)])
ylabel('Approx. Depth (m)')
xlim(xl)

%ylim([0.85 2.25]);


print -dpng -r300 'Data Aug 5 to 12.png'