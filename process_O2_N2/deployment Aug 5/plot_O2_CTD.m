load PME_O2_cc.mat;
load('S12071.mat');
% add in loading of the Pro-Oceanus data

%T3 = outerjoin(PME_O2_cc(1),PME_O2_cc(2),'MergeKeys', true);
%%

% convert times from UTC to EST (4 hr offset)
PME_O2_cc.datetime_EST = PME_O2_cc.datetime - 4/24;

% calculate O2 in umol/L from mg/L
MW_O2 = 31.998; % molar weight of O2 in g/mol
mg_per_g = 1000; % milligrams per gram
umol_per_mol = 1e6;
PME_O2_cc.DO_umolL = PME_O2_cc.DO ./ MW_O2 ./ mg_per_g .* umol_per_mol;

figure(1)
clf; hold on;

% oxygen plot
subplot(4,1,1)
hold on; box on;
plot(PME_O2_cc.datetime_EST,PME_O2_cc.DO_umolL)
plot([datetime(2022,8,5,11,15,0) datetime(2022,8,12,11,0,0)], [215 215])
xlim([datetime(2022,8,5,11,15,0) datetime(2022,8,12,11,0,0)])
ylim([100 350])
set(gca,'xticklabel',[]) % remove date label from all but bottom
%xlim([737838 737848])
ylabel("O_2 (\mumol/L)")

% temperature plot
% there is temperature from multiple sensors
% SO = Star-Oddi (CTD sensor)
% PME = Precision Measurement Engineering (O2 sensor)
% PO = Pro-Oceanus (TDGP sensor)
subplot(4,1,2)
hold on; box on;
plot(S12071.sample_date,S12071.temp); % star-oddi temperature
plot(PME_O2_cc.datetime_EST,PME_O2_cc.T)
legend('SO T','PME T')
ylim([22.5 29])
xlim([datetime(2022,8,5,11,15,0) datetime(2022,8,12,11,0,0)])
set(gca,'xticklabel',[]) % remove date label from all but bottom
ylabel("Temp (^oC)")

% salinity plot
subplot(4,1,3)
hold on; box on;
plot(S12071.sample_date,S12071.sal)
ylim([24 30]) 
xlim([datetime(2022,8,5,11,15,0) datetime(2022,8,12,11,0,0)])
set(gca,'xticklabel',[]) % remove date label from all but bottom
ylabel('Sal (PSS)')

% depth plot
% note: for now I am applying a manual offset of 1 m to the depth data
% based on the depth we measured with the rope. 
% labeled as Approx Depth
% we need to verify the offset more accurately
subplot(4,1,4)
plot(S12071.sample_date,S12071.depth+1)
xlim([datetime(2022,8,5,11,15,0) datetime(2022,8,12,11,0,0)])
ylabel('Approx. Depth (m)')
ylim([0.85 2.25]);

print -dpng -r300 'Data Week 1 umolL.png'