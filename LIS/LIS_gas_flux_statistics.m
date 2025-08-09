% stnlist from west to east
stnlist = {'EXR1', 'EXRX', 'MID3', 'MID4', 'MID5', 'WLIS', 'WLI6'};

km_transect = [0    1.4482    5.4270    9.0377   13.0335   16.2552   18.4318];
%%
% load August data
load LIS_gas_flux_Aug.mat;
gfA = LIS_gas_flux_Aug;
gfA.station_str = string(gfA.station);

sl = string(stnlist);
gfA_avg = table;

gfA_avg.station = sl'; % create table for average flux data
gfA_avg.lat = nan(numel(stnlist),1);
gfA_avg.lon = nan(numel(stnlist),1);
gfA_avg.n2o_nmolkg = nan(numel(stnlist),1);
gfA_avg.n2o_std_nmolkg = nan(numel(stnlist),1);
gfA_avg.n2o_eq_nmolkg = nan(numel(stnlist),1);
gfA_avg.ch4_nmolkg = nan(numel(stnlist),1);
gfA_avg.ch4_std_nmolkg = nan(numel(stnlist),1);
gfA_avg.ch4_eq_nmolkg = nan(numel(stnlist),1);
gfA_avg.F_CH4_15 = nan(numel(stnlist),1);
gfA_avg.F_CH4_15_std = nan(numel(stnlist),1);
gfA_avg.F_N2O_15 = nan(numel(stnlist),1);
gfA_avg.F_N2O_15_std = nan(numel(stnlist),1);



for i = 1:length(stnlist);
    substr = stnlist(i);
    strfindResult = strfind(gfA.station_str, substr);
    indices = find(~cellfun('isempty', strfindResult));
    gfA_avg.lat(i) = mean(gfA.lat(indices));
    gfA_avg.lon(i) = mean(gfA.lon(indices));  
    gfA_avg.n2o_nmolkg(i) = mean(gfA.n2o_nmolkg(indices));  
    % if there are multiple samples for the station, take standard
    % deviation of repeat measurements
    % else, take the standard deviation from duplicates
    if numel(indices) > 1
        gfA_avg.n2o_std_nmolkg(i) = std(gfA.n2o_nmolkg(indices)); 
    else
        gfA_avg.n2o_std_nmolkg(i) = mean(gfA.n2o_std_nmolkg(indices)); 
    end;    

    gfA_avg.n2o_eq_nmolkg(i) = mean(gfA.n2o_eq_nmolkg(indices)); 

    gfA_avg.ch4_nmolkg(i) = mean(gfA.ch4_nmolkg(indices));  
    % if there are multiple samples for the station, take standard
    % deviation of repeat measurements
    % else, take the standard deviation from duplicates
    if numel(indices) > 1
        gfA_avg.ch4_std_nmolkg(i) = std(gfA.ch4_nmolkg(indices)); 
    else
        gfA_avg.ch4_std_nmolkg(i) = mean(gfA.ch4_std_nmolkg(indices)); 
    end;

    gfA_avg.ch4_eq_nmolkg(i) = mean(gfA.ch4_eq_nmolkg(indices));       
    gfA_avg.F_N2O_15(i) = mean(gfA.F_N2O_15(indices));  
    gfA_avg.F_N2O_15_std(i) = std(gfA.F_N2O_15(indices)); 
    gfA_avg.F_CH4_15(i) = mean(gfA.F_CH4_15(indices));  
    gfA_avg.F_CH4_15_std(i) = std(gfA.F_CH4_15(indices)); 
    gfA_avg.k_wt_15_N2O(i) = mean(gfA.k_wt_15_N2O(indices));  
    gfA_avg.k_wt_15_CH4(i) = mean(gfA.k_wt_15_CH4(indices));     
end;

gfA_avg.DCH4 = (gfA_avg.ch4_nmolkg - gfA_avg.ch4_eq_nmolkg)./gfA_avg.ch4_eq_nmolkg .* 100;
gfA_avg.DN2O = (gfA_avg.n2o_nmolkg - gfA_avg.n2o_eq_nmolkg)./gfA_avg.n2o_eq_nmolkg .* 100;

save LIS_gas_flux_Aug_avg.mat gfA_avg gfA;
%%

%%

% load Oct data
load LIS_gas_flux_Oct.mat;
gfO = LIS_gas_flux_Oct;
gfO.station_str = string(gfO.station);

sl = string(stnlist);
gfO_avg = table;

gfO_avg.station = sl'; % create table for average flux data
gfO_avg.lat = nan(numel(stnlist),1);
gfO_avg.lon = nan(numel(stnlist),1);
gfO_avg.n2o_nmolkg = nan(numel(stnlist),1);
gfO_avg.n2o_std_nmolkg = nan(numel(stnlist),1);
gfO_avg.n2o_eq_nmolkg = nan(numel(stnlist),1);
gfO_avg.ch4_nmolkg = nan(numel(stnlist),1);
gfO_avg.ch4_std_nmolkg = nan(numel(stnlist),1);
gfO_avg.ch4_eq_nmolkg = nan(numel(stnlist),1);
gfO_avg.F_CH4_15 = nan(numel(stnlist),1);
gfO_avg.F_CH4_15_std = nan(numel(stnlist),1);
gfO_avg.F_N2O_15 = nan(numel(stnlist),1);
gfO_avg.F_N2O_15_std = nan(numel(stnlist),1);



for i = 1:length(stnlist);
    substr = stnlist(i);
    strfindResult = strfind(gfO.station_str, substr);
    indices = find(~cellfun('isempty', strfindResult));
    gfO_avg.lat(i) = mean(gfO.lat(indices));
    gfO_avg.lon(i) = mean(gfO.lon(indices));  
    gfO_avg.n2o_nmolkg(i) = mean(gfO.n2o_nmolkg(indices));  
    % if there are multiple samples for the station, take standard
    % deviation of repeat measurements
    % else, take the standard deviation from duplicates
    if numel(indices) > 1
        gfO_avg.n2o_std_nmolkg(i) = std(gfO.n2o_nmolkg(indices)); 
    else
        gfO_avg.n2o_std_nmolkg(i) = mean(gfO.n2o_std_nmolkg(indices)); 
    end;     
    gfO_avg.n2o_eq_nmolkg(i) = mean(gfO.n2o_eq_nmolkg(indices));       
    gfO_avg.ch4_nmolkg(i) = mean(gfO.ch4_nmolkg(indices));  
    % if there are multiple samples for the station, take standard
    % deviation of repeat measurements
    % else, take the standard deviation from duplicates
    if numel(indices) > 1
        gfO_avg.ch4_std_nmolkg(i) = std(gfO.ch4_nmolkg(indices)); 
    else
        gfO_avg.ch4_std_nmolkg(i) = mean(gfO.ch4_std_nmolkg(indices)); 
    end;

    gfO_avg.ch4_eq_nmolkg(i) = mean(gfO.ch4_eq_nmolkg(indices));       
    gfO_avg.F_N2O_15(i) = mean(gfO.F_N2O_15(indices));  
    gfO_avg.F_N2O_15_std(i) = std(gfO.F_N2O_15(indices)); 
    gfO_avg.F_CH4_15(i) = mean(gfO.F_CH4_15(indices));  
    gfO_avg.F_CH4_15_std(i) = std(gfO.F_CH4_15(indices)); 
    gfO_avg.k_wt_15_N2O(i) = mean(gfO.k_wt_15_N2O(indices));  
    gfO_avg.k_wt_15_CH4(i) = mean(gfO.k_wt_15_CH4(indices));     
end;

gfO_avg.DCH4 = (gfO_avg.ch4_nmolkg - gfO_avg.ch4_eq_nmolkg)./gfO_avg.ch4_eq_nmolkg .* 100;
gfO_avg.DN2O = (gfO_avg.n2o_nmolkg - gfO_avg.n2o_eq_nmolkg)./gfO_avg.n2o_eq_nmolkg .* 100;

save LIS_gas_flux_Oct_avg.mat gfO_avg gfO;

% load May data
load LIS_gas_flux_May.mat;
gfM = LIS_gas_flux_May;
gfM.station_str = string(gfM.station);

sl = string(stnlist);
gfM_avg = table;

gfM_avg.station = sl'; % create table for average flux data
gfM_avg.lat = nan(numel(stnlist),1);
gfM_avg.lon = nan(numel(stnlist),1);
gfM_avg.n2o_nmolkg = nan(numel(stnlist),1);
    % if there are multiple samples for the station, take standard
    % deviation of repeat measurements
    % else, take the standard deviation from duplicates
    if numel(indices) > 1
        gfM_avg.n2o_std_nmolkg(i) = std(gfM.n2o_nmolkg(indices)); 
    else
        gfM_avg.n2o_std_nmolkg(i) = mean(gfM.n2o_std_nmolkg(indices)); 
    end;    
gfM_avg.n2o_eq_nmolkg = nan(numel(stnlist),1);
gfM_avg.ch4_nmolkg = nan(numel(stnlist),1);
    % if there are multiple samples for the station, take standard
    % deviation of repeat measurements
    % else, take the standard deviation from duplicates
    if numel(indices) > 1
        gfM_avg.ch4_std_nmolkg(i) = std(gfM.ch4_nmolkg(indices)); 
    else
        gfM_avg.ch4_std_nmolkg(i) = mean(gfM.ch4_std_nmolkg(indices)); 
    end;

gfM_avg.ch4_eq_nmolkg = nan(numel(stnlist),1);
gfM_avg.F_CH4_15 = nan(numel(stnlist),1);
gfM_avg.F_CH4_15_std = nan(numel(stnlist),1);
gfM_avg.F_N2O_15 = nan(numel(stnlist),1);
gfM_avg.F_N2O_15_std = nan(numel(stnlist),1);


for i = 1:length(stnlist);
    substr = stnlist(i);
    strfindResult = strfind(gfM.station_str, substr);
    indices = find(~cellfun('isempty', strfindResult));
    gfM_avg.lat(i) = mean(gfM.lat(indices));
    gfM_avg.lon(i) = mean(gfM.lon(indices));  
    gfM_avg.n2o_nmolkg(i) = mean(gfM.n2o_nmolkg(indices));  
    gfM_avg.n2o_std_nmolkg(i) = std(gfM.n2o_nmolkg(indices)); 
    gfM_avg.n2o_eq_nmolkg(i) = mean(gfM.n2o_eq_nmolkg(indices));    
    gfM_avg.ch4_nmolkg(i) = mean(gfM.ch4_nmolkg(indices));  
    gfM_avg.ch4_std_nmolkg(i) = std(gfM.ch4_nmolkg(indices)); 
    gfM_avg.ch4_eq_nmolkg(i) = mean(gfM.ch4_eq_nmolkg(indices));      
    gfM_avg.F_N2O_15(i) = mean(gfM.F_N2O_15(indices));  
    gfM_avg.F_N2O_15_std(i) = std(gfM.F_N2O_15(indices)); 
    gfM_avg.F_CH4_15(i) = mean(gfM.F_CH4_15(indices));  
    gfM_avg.F_CH4_15_std(i) = std(gfM.F_CH4_15(indices)); 
    gfM_avg.k_wt_15_N2O(i) = mean(gfM.k_wt_15_N2O(indices));  
    gfM_avg.k_wt_15_CH4(i) = mean(gfM.k_wt_15_CH4(indices));     
end;

gfM_avg.DCH4 = (gfM_avg.ch4_nmolkg - gfM_avg.ch4_eq_nmolkg)./gfM_avg.ch4_eq_nmolkg .* 100;
gfM_avg.DN2O = (gfM_avg.n2o_nmolkg - gfM_avg.n2o_eq_nmolkg)./gfM_avg.n2o_eq_nmolkg .* 100;

save LIS_gas_flux_May_avg.mat gfM_avg gfM;

%%
[median(gfA_avg.DCH4) min(gfA_avg.DCH4) max(gfA_avg.DCH4)]
[median(gfO_avg.DCH4) min(gfO_avg.DCH4) max(gfO_avg.DCH4)]
[median(gfM_avg.DCH4) min(gfM_avg.DCH4) max(gfM_avg.DCH4)]

[median(gfA_avg.F_CH4_15) min(gfA_avg.F_CH4_15) max(gfA_avg.F_CH4_15)]
[median(gfO_avg.F_CH4_15) min(gfO_avg.F_CH4_15) max(gfO_avg.F_CH4_15)]
[median(gfM_avg.F_CH4_15) min(gfM_avg.F_CH4_15) max(gfM_avg.F_CH4_15)]

%%
[median(gfA_avg.DN2O) min(gfA_avg.DN2O) max(gfA_avg.DN2O)]
[median(gfO_avg.DN2O) min(gfO_avg.DN2O) max(gfO_avg.DN2O)]
[median(gfM_avg.DN2O) min(gfM_avg.DN2O) max(gfM_avg.DN2O)]

[median(gfA_avg.F_N2O_15) min(gfA_avg.F_N2O_15) max(gfA_avg.F_N2O_15)]
[median(gfO_avg.F_N2O_15) min(gfO_avg.F_N2O_15) max(gfO_avg.F_N2O_15)]
[median(gfM_avg.F_N2O_15) min(gfM_avg.F_N2O_15) max(gfM_avg.F_N2O_15)]


%%
[mean(gfA_avg.F_N2O_15) std(gfA_avg.F_N2O_15)
mean(gfO_avg.F_N2O_15) std(gfO_avg.F_N2O_15)
mean(gfM_avg.F_N2O_15) std(gfM_avg.F_N2O_15)]


[mean(gfA_avg.DN2O) std(gfA_avg.DN2O)
mean(gfO_avg.DN2O) std(gfO_avg.DN2O)
mean(gfM_avg.DN2O) std(gfM_avg.DN2O)]

%%
[mean(gfA_avg.F_CH4_15) std(gfA_avg.F_CH4_15)
mean(gfO_avg.F_CH4_15) std(gfO_avg.F_CH4_15)
mean(gfM_avg.F_CH4_15) std(gfM_avg.F_CH4_15)]


[mean(gfA_avg.DCH4) std(gfA_avg.DCH4)
mean(gfO_avg.DCH4) std(gfO_avg.DCH4)
mean(gfM_avg.DCH4) std(gfM_avg.DCH4)]


%%

figure(1)
clf; 
subplot(2,2,1)
hold on; box on;
set(gca,'tickdir','out');
set(gca,'xlim',[-0.5 19]);
%set(gca,'xticklabel',stnlist);
set(gca,'xtick',km_transect)

%plot(gfA_avg.ch4_nmolkg - gfA_avg.ch4_eq_nmolkg,'-o');
%plot(gfO_avg.ch4_nmolkg - gfO_avg.ch4_eq_nmolkg,'-o');
plot(km_transect, gfA_avg.DCH4,'-o');
plot(km_transect, gfO_avg.DCH4,'-o');
plot(km_transect, gfM_avg.DCH4,'-o');

ax1.XAxisLocation = 'bottom';
ax1.YAxisLocation = 'left';
ax1.Box = 'off'; % Optional: remove the box to avoid double ticks
ax1.XTick = 0:5:15;         % Bottom ticks

legend('Aug','Oct','May');
%ylabel('\DeltaCH_4 (nmol kg^{-1})')
ylabel('\DeltaCH_4 (%)');
xlabel('Transect distance [km]');

subplot(2,2,2)
hold on; box on;
set(gca,'tickdir','out');
%set(gca,'xlim',[-0.2 7.2]);
set(gca,'xticklabel',stnlist);
ylabel('CH_4 flux (\mumol m^{-2} d^{-1})')

plot(gfA_avg.F_CH4_15,'-o');
plot(gfO_avg.F_CH4_15,'-o');
plot(gfM_avg.F_CH4_15,'-o');
%legend('Aug','Oct','May');

subplot(2,2,3)
hold on; box on;
set(gca,'tickdir','out');
%set(gca,'xlim',[-0.2 7.2]);
set(gca,'xticklabel',stnlist);

plot(gfA_avg.DN2O,'-o');
plot(gfO_avg.DN2O,'-o');
plot(gfM_avg.DN2O,'-o');
%legend('Aug','Oct','May');
%ylabel('N_2O (nmol kg^{-1})')
ylabel('\DeltaN_2O (%)');


subplot(2,2,4)
hold on; box on;
set(gca,'tickdir','out');
%set(gca,'xlim',[-0.2 7.2]);
set(gca,'xticklabel',stnlist);
ylabel('N_2O flux (\mumol m^{-2} d^{-1})')

plot(gfA_avg.F_N2O_15,'-o');
plot(gfO_avg.F_N2O_15,'-o');
plot(gfM_avg.F_N2O_15,'-o');
%legend('Aug','Oct','May');


%%
plot(gfO.F_CH4_15,'-o');
plot(gfM.F_CH4_15,'-o');

% make a plot of station vs flux and station vs sat anomaly
% average the results for repeat stations and show stdev