% load mims data
ws = 5; % window size: number of points to include when calculating movmean, movstd etc.
file_path = 'data/Membrane Test Data/';
mat_fn = '20251111_1231_mims.mat';

load([file_path mat_fn])

% if only want to select part of the data, use this line, otherwise comment
% out
%mims = mims(500:end,:);

% background corrected values
% calculated at start in case the bkg masses change later
mims.m22b = (mims.m22 - mims.m23);
mims.m36b = (mims.m36 - mims.m33p5);
mims.m38b = (mims.m38 - mims.m33p5);
mims.m84b = (mims.m84 - mims.m88);

% create new ratios if not already defined in file
% ratios to 84
mims.r2284 = mims.m22 ./ mims.m84;

% ratios to 84 with bkg subtraction
mims.r2284b = mims.m22b ./ mims.m84b;

% ratios to 36
mims.r8436 = mims.m84 ./ mims.m36;
mims.r2236 = mims.m22 ./ mims.m36;
mims.r3836 = mims.m38 ./ mims.m36;

% ratios to 36 with bkg subtraction
mims.r8436b = mims.m84b ./ mims.m36b;
mims.r2236b = mims.m22b ./ mims.m36b;
mims.r3836b = mims.m38b ./ mims.m36b;

% ratios to 38
mims.r8438 = mims.m84 ./ mims.m38;
mims.r2238 = mims.m22 ./ mims.m38;
mims.r3638 = mims.m36 ./ mims.m38;

% ratios to 38 with bkg subtraction
mims.r8438b = mims.m84b ./ mims.m38b;
mims.r2238b = mims.m22b ./ mims.m38b;
mims.r3638b = mims.m36b ./ mims.m38b;

% moving means
mims.mm_m22 = movmean(mims.m22, ws);
mims.mm_m36 = movmean(mims.m36, ws);
mims.mm_m38 = movmean(mims.m38, ws);
mims.mm_m84 = movmean(mims.m84, ws);

mims.mm_r2238 = movmean(mims.r2238, ws);
mims.mm_r3638 = movmean(mims.r3638, ws);
mims.mm_r8438 = movmean(mims.r8438, ws);

mims.mm_r2284 = movmean(mims.r2284, ws);

mims.mm_r2236 = movmean(mims.r2236, ws);
mims.mm_r3836 = movmean(mims.r3836, ws);
mims.mm_r8436 = movmean(mims.r8436, ws);

% moving means with bkg subtraction
mims.mm_m22 = movmean(mims.m22b, ws);
mims.mm_m36 = movmean(mims.m36b, ws);
mims.mm_m38 = movmean(mims.m38b, ws);
mims.mm_m84 = movmean(mims.m84b, ws);

mims.mm_r2238b = movmean(mims.r2238b, ws);
mims.mm_r3638b = movmean(mims.r3638b, ws);
mims.mm_r8438b = movmean(mims.r8438b, ws);

mims.mm_r2284b = movmean(mims.r2284b, ws);

mims.mm_r2236b = movmean(mims.r2236b, ws);
mims.mm_r3836b = movmean(mims.r3836b, ws);
mims.mm_r8436b = movmean(mims.r8436b, ws);

% coefficient of variance based on moving mean and moving std
mims.cv_m22 = movstd(mims.m22, ws) ./ mims.mm_m22;
mims.cv_m36 = movstd(mims.m36 ,ws) ./ mims.mm_m36;
mims.cv_m38 = movstd(mims.m38, ws)  ./ mims.mm_m38;
mims.cv_m84 = movstd(mims.m84, ws)  ./ mims.mm_m84;

mims.cv_r2238 = movstd(mims.r2238, ws) ./ mims.mm_r2238;
mims.cv_r3638 = movstd(mims.r3638 ,ws) ./ mims.mm_r3638;
mims.cv_r8438 = movstd(mims.r8438, ws)  ./ mims.mm_r8438;
mims.cv_r2284 = movstd(mims.r2284, ws)  ./ mims.mm_r2284;

mims.cv_r2236 = movstd(mims.r2236, ws) ./ mims.mm_r2236;
mims.cv_r3836 = movmean(mims.r3836, ws) ./ mims.mm_r3836;
mims.cv_r8436 = movmean(mims.r8436, ws) ./ mims.mm_r8436;

% coefficient of variance w/ bkg subtraction
mims.cv_m22b = movstd(mims.m22b, ws) ./ mims.mm_m22b;
mims.cv_m36b = movstd(mims.m36b, ws) ./ mims.mm_m36b;
mims.cv_m38b = movstd(mims.m38b, ws)  ./ mims.mm_m38b;
mims.cv_m84b = movstd(mims.m84b, ws)  ./ mims.mm_m84b;

mims.cv_r2238b = movstd(mims.r2238b, ws) ./ mims.mm_r2238b;
mims.cv_r3638b = movstd(mims.r3638b ,ws) ./ mims.mm_r3638b;
mims.cv_r8438b = movstd(mims.r8438b, ws)  ./ mims.mm_r8438b;
mims.cv_r2284b = movstd(mims.r2284b, ws)  ./ mims.mm_r2284b;

mims.cv_r2236b = movstd(mims.r2236b, ws) ./ mims.mm_r2236b;
mims.cv_r3836b = movmean(mims.r3836b, ws) ./ mims.mm_r3836b;
mims.cv_r8436b = movmean(mims.r8436b, ws) ./ mims.mm_r8436b;

%%
% plot of ion currents and pressure, rows are:
% 22
% 84
% 36
% 38 
% TP
figure(1)
clf; 

tiledlayout(5,1, 'TileSpacing', 'tight', 'Padding', 'tight');

nexttile;
hold on; box on;
plot(mims.dt,mims.m22,'r');
%plot(mims.dt,mims.m23,'r');
plot(mims.dt,mims.m22-mims.m23,'k');
plot(mims.dt, mims.mm_m22,'c');
legend('22','22-23', '22 movmean', 'location','southeast');
ylabel('ion current (A)');
axis tight;
xticklabels([]);


nexttile;
hold on; box on;
plot(mims.dt,mims.m84,'r');
%plot(mims.dt,mims.m88,'b');
plot(mims.dt,mims.m84b,'k');
plot(mims.dt, mims.mm_m84,'c');
legend('84','84-88', '84 movmean','location','southeast');
ylabel('ion current (A)');
axis tight;
xticklabels([]);

nexttile;
hold on; box on;
plot(mims.dt,mims.m36,'r');
%plot(mims.dt,mims.m33p5,'b');
plot(mims.dt,mims.m36b,'k');
plot(mims.dt, mims.mm_m36,'c');
legend('36','36-33.5', '36 movmean','location','southeast');
ylabel('ion current (A)');
axis tight;
xticklabels([]);


nexttile;
hold on; box on;
plot(mims.dt,mims.m38,'r');
%plot(mims.dt,mims.m33p5,'b');
plot(mims.dt,mims.m38b,'k');
plot(mims.dt, mims.mm_m38,'c');
legend('38','38-33.5', '38 movmean', 'location','southeast');
ylabel('ion current (A)');
axis tight;
xticklabels([]);

nexttile;
hold on; box on;
plot(mims.dt,mims.tp,'r')
legend('TP','location','southeast');
ylabel('pressure (mbar)');
axis tight;
%xticklabels([]);

%%
% plot of ratios and pressure, rows are:
% 22
% 84
% 36
% 38 
% TP
figure(2)
clf; 
tiledlayout(5,1, 'TileSpacing', 'tight', 'Padding', 'tight');

nexttile;
hold on; box on;
plot(mims.dt,mims.r2238,'r');
plot(mims.dt, mims.mm_r2238,'k');
legend('22/38','22/38 mm', 'location','southeast');
ylabel('ratio 22/38');
axis tight;
xticklabels([]);


nexttile;
hold on; box on;
plot(mims.dt,mims.r8438,'r');
plot(mims.dt, mims.mm_r8438,'k');
legend('84/38','84/38 mm','location','southeast');
ylabel('ratio 84/38');
axis tight;
xticklabels([]);

nexttile;
hold on; box on;
plot(mims.dt,mims.r3638,'r');
plot(mims.dt, mims.mm_r3638,'k');
legend('36/38','36/38 mm','location','southeast');
ylabel('ratio 36/38');
axis tight;
xticklabels([]);


nexttile;
hold on; box on;
plot(mims.dt,mims.r2284,'r');
plot(mims.dt, mims.mm_r2284,'k');
legend('22/84','22/84 mm', 'location','southeast');
ylabel('ratio 22/84');
axis tight;
xticklabels([]);

nexttile;
hold on; box on;
plot(mims.dt,mims.tp,'r');
legend('TP', 'location','southeast');
ylabel('pressure (mbar)');
axis tight;
%xticklabels([]);
%%
nexttile;
hold on; box on;
plot(mims.dt,mims.tp,'r')
legend('TP','location','southeast');
ylabel('pressure (mbar)');
axis tight;
%xticklabels([]);

%%
figure(2)
clf; 

nexttile;
hold on; box on;
plot(mims.dt,mims.r2236,'b');
plot(mims.dt,mims.r2238,'r');
plot(mims.dt,mims.r2236b,'k');
plot(mims.dt,mims.r2238b,'g');

legend('22/36', '22/38', '22/36b', '22/38b');

%%

nexttile;
hold on; box on;
%plot(mims.dt,mims.m84,'b');
plot(mims.dt,mims.m84,'r');
plot(mims.dt,mims.m88,'b');
plot(mims.dt,mims.m84-mims.m88,'k');
legend('84','88','84-88');

nexttile;
hold on; box on;
plot(mims.dt,mims.m36,'r');
plot(mims.dt,mims.m33p5,'b');
plot(mims.dt,mims.m36-mims.m33p5,'k');
legend('36','33.5','36-33.5');


nexttile;
hold on; box on;
plot(mims.dt,mims.m38,'r');
plot(mims.dt,mims.m36,'b');
plot(mims.dt,mims.m33p5,'g');
legend('38','36','38-33.5');



%%

n = 172:185; % 2 min average
n = 839:852; % 2 min average
n = 777:790;
n = 777:782;

%%

% Calculate the coefficients of variation for the selected ranges
cv_m38 = std(mims.m38(n)) ./ mean(mims.m38(n))
cv_m36 = std(mims.m36(n)) ./ mean(mims.m36(n))
cv_m84 = std(mims.m84(n)) ./ mean(mims.m84(n))
cv_m22 = std(mims.m22(n)) ./ mean(mims.m22(n))

