knots_to_meters_per_second = 0.514444;
load exrx2023-2024_met.mat;
exrxMet.windSpd_Kts(exrxMet.windSpd_Kts==0) = NaN;
exrxMet.windSpd_MS = exrxMet.windSpd_Kts.* knots_to_meters_per_second;

mbar_to_atm = 1./1013.25; % convert pressure data from mbar to atm for use in gas sol eqns

exrxMet.baroPress_atm = exrxMet.baroPress_Avg .* mbar_to_atm;

exrxMet_adj = exrxMet;


save exrxMet_adj.mat exrxMet_adj;
%%
load wlis2023-2024_met.mat;
A=find(wlisMet.windSpd_Kts>50);
wlisMet.windSpd_Kts(A) = nan;
wlisMet.windSpd_Kts(wlisMet.windSpd_Kts==0) = NaN;
wlisMet.windSpd_MS = wlisMet.windSpd_Kts.* knots_to_meters_per_second;


load clis2022-2024_met.mat;
clisMet.windSpd_Kts(clisMet.windSpd_Kts==0) = NaN;
clisMet.windSpd_MS = clisMet.windSpd_Kts.* knots_to_meters_per_second;

clisMet_adj = clisMet;
save clisMet_adj.mat clisMet_adj;

%%

% cruise dates
% Aug 2-3 2023 - EXRX and CLIS
% Oct 19-20 2023 - EXRX and CLIS
% May 22-23, 2024 - EXRX and CLIS only

d1 = datetime(2023,9,10,0,0,0);
d2 = datetime(2023,10,31,0,0,0);
Ai = (exrxMet.TIMESTAMP>=d1) & (exrxMet.TIMESTAMP<=d2) & ~isnan(exrxMet.windSpd_MS);
A = find(Ai);
if(isempty(A))
    disp('no values for exrx')
end;

Bi = (wlisMet.TIMESTAMP>=d1) & (wlisMet.TIMESTAMP<=d2) & ~isnan(wlisMet.windSpd_MS);
B = find(Bi);
if(isempty(B))
    disp('no values for wlis')
end;

Ci = (clisMet.TIMESTAMP>=d1) & (clisMet.TIMESTAMP<=d2)  & ~isnan(clisMet.windSpd_MS);
C = find(Ci);
if(isempty(C))
    disp('no values for clis')
end;
%%
% note
figure(1)
clf; hold on;
plot(exrxMet.TIMESTAMP(A),exrxMet.windSpd_MS(A));
%plot(wlisMet.TIMESTAMP(B),wlisMet.windSpd_MS(B));
plot(clisMet.TIMESTAMP(C),clisMet.windSpd_MS(C));
%legend('exrx','wlis','clis')
ylim([0,15]);


%%
dl = datetime(2022,1,1,0,0,0):15/60/24:datetime(2025,1,1,0,0,0);
dl = datetime(2023,9,1,0,0,0):15/60/24:datetime(2023,11,1,0,0,0);

% find the unique values in exrx
[unique_times_exrx, ~, idx] = unique(exrxMet.TIMESTAMP);  % idx gives the group index for each key
num_keys = numel(unique_times_exrx);
avg_wspd_exrx = accumarray(idx(:), exrxMet.windSpd_MS(:), [], @mean);

i_exrx = interp1(unique_times_exrx,avg_wspd_exrx,dl); % interpolate to fill in the gaps

wspd_exrx_dl = nan(size(dl));

for i = 1:length(dl)
    A = find(unique_times_exrx == dl(i));
    if ~isempty(A)
        wspd_exrx_dl(i) = avg_wspd_exrx(A);
    end;
end;

% get a list of the nan values
for i = 1:length(i_exrx)
    A = find(unique_times_exrx == dl(i));
    if ~isempty(A)
        wspd_exrx_dl(i) = avg_wspd_exrx(A);
    end;
end;

%%
figure(12)
clf; hold on;
plot(dl,wspd_exrx_dl);

%%
[unique_times_wlis, ~, idx] = unique(wlisMet.TIMESTAMP);  % idx gives the group index for each key
num_keys = numel(unique_times_wlis);
avg_wspd_wlis = accumarray(idx(:), wlisMet.windSpd_MS(:), [], @mean);

i_wlis = interp1(unique_times_wlis,avg_wspd_wlis,dl); % interpolate to fill in the 

[unique_times_clis, ~, idx] = unique(clisMet.TIMESTAMP);  % idx gives the group index for each key
num_keys = numel(unique_times_clis);
avg_wspd_clis = accumarray(idx(:), clisMet.windSpd_MS(:), [], @mean);

i_clis = interp1(unique_times_clis,avg_wspd_clis,dl);

wspd_clis_dl = nan(size(dl));

for i = 1:length(dl)
    A = find(unique_times_clis == dl(i));
    if ~isempty(A)
        wspd_clis_dl(i) = avg_wspd_clis(A);
    end;
end;

%%
x = i_exrx;  % data
win = 25;    % size of the moving window (must be odd)

% Compute moving statistics
mov_mean = movmean(x, win, 'omitnan');
mov_std = movstd(x, win, 'omitnan');

% Find spikes
spike_idx = x > mov_mean + 2 * mov_std; % find spikes that are more than 2x the running stdev

% Option 1: Replace spikes with NaN
x_cleaned = x;
x_cleaned(spike_idx) = NaN;

% Option 2: Replace spikes with moving mean (smoother option)
x_smoothed = x;
x_smoothed(spike_idx) = mov_mean(spike_idx);

figure(3)
clf; hold on;
plot(dl,i_exrx,'.');
plot(dl,x_smoothed);

%%
% Linear fit (degree 1)
pi = ~isnan(wspd_clis_dl) & ~isnan(wspd_exrx_dl);
p = polyfit(wspd_clis_dl(pi), wspd_exrx_dl(pi), 1);

% p(1) = slope, p(2) = intercept
slope = p(1);
intercept = p(2);

figure(4)
clf; hold on;
scatter(wspd_clis_dl,wspd_exrx_dl);
plot(wspd_clis_dl(pi),wspd_clis_dl(pi)*slope + intercept,'.');
ylabel('exrx');
xlabel('clis');

figure(5)
clf; hold on;
plot(dl,wspd_exrx_dl);
plot(dl,wspd_clis_dl);
legend('exrx','clis');

%%
figure(2)
clf; hold on;
plot(dl,i_exrx);
plot(dl,i_wlis);
plot(dl,i_clis);

figure(3)
clf; hold on;
plot(dl,i_wlis./i_exrx);
plot(dl,i_wlis./i_clis);
ylim([0 50]);

%%
A = i_wlis./i_exrx;
B = i_wlis./i_clis;

%%

mean(A<=5)

mean(Ai);
figure(3)
clf;
subplot(1,2,1)
hold on;
histogram(A,'binedges',[0:0.05:5]);

subplot(1,2,2)
hold on;
histogram(A);

%nanmean(A(1000:2000))

%%

