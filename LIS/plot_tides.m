% cruise 1: August 2-3, 2023 = julian day 214.5, year 2023.588
% cruise 2: October 19-20, 2023 = julian day 292.5, year 2023.801 
% cruise 3: May 22-23, 2024 = julian day 142.5, year 2024.390

UTC_to_local = -4/24;

load Bridgeport_Aug2023;
BM = Bridgeport_Aug2023;

BM.Date.Format = 'yyyy/MM/dd HH:mm';
BM.TimeGMT.Format = 'yyyy/MM/dd HH:mm';
BM.Datetime = BM.Date + timeofday(BM.TimeGMT);
BM.Datetime_local = BM.Datetime + UTC_to_local;

Bridgeport_Aug2023 = BM;
save Bridgeport_Aug2023.mat Bridgeport_Aug2023;

%%

load KingsPoint_May2024;
KM = KingsPoint_May2024;

KM.Date.Format = 'yyyy/MM/dd HH:mm';
KM.TimeGMT.Format = 'yyyy/MM/dd HH:mm';
KM.Datetime = KM.Date + timeofday(KM.TimeGMT);
KM.Datetime_local = KM.Datetime + UTC_to_local;

KingsPoint_May2024 = KM;
save KingsPoint_May2024.mat KingsPoint_May2024;
%%

figure(1);
clf; hold on;
plot(BM.Datetime_local,BM.Verifiedft);
plot(KM.Datetime_local,KM.Verifiedft);