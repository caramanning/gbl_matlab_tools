load exrx2023-2024_met.mat;
exrxMet.windSpd_Kts(exrxMet.windSpd_Kts==0) = NaN;

load wlis2023-2024_met.mat;
A=find(wlisMet.windSpd_Kts>50);
wlisMet.windSpd_Kts(A) = nan;
wlisMet.windSpd_Kts(wlisMet.windSpd_Kts==0) = NaN;


load clis2022-2024_met.mat;
clisMet.windSpd_Kts(clisMet.windSpd_Kts==0) = NaN;

%%

% cruise dates
% Aug 2-3 2023 - EXRX and CLIS
% Oct 19-20 2023 - EXRX and CLIS
% May 22-23, 2024 - 

d1 = datetime(2023,9,10,0,0,0);
d2 = datetime(2023,10,31,0,0,0);
Ai = (exrxMet.TIMESTAMP>=d1) & (exrxMet.TIMESTAMP<=d2) & ~isnan(exrxMet.windSpd_Kts);
A = find(Ai);
if(isempty(A))
    disp('no values for exrx')
end;

Bi = (wlisMet.TIMESTAMP>=d1) & (wlisMet.TIMESTAMP<=d2) & ~isnan(wlisMet.windSpd_Kts);
B = find(Bi);
if(isempty(B))
    disp('no values for wlis')
end;

Ci = (clisMet.TIMESTAMP>=d1) & (clisMet.TIMESTAMP<=d2)  & ~isnan(clisMet.windSpd_Kts);
C = find(Ci);
if(isempty(C))
    disp('no values for clis')
end;
%%
% note
figure(1)
clf; hold on;
plot(exrxMet.TIMESTAMP(A),exrxMet.windSpd_Kts(A));
%plot(wlisMet.TIMESTAMP(B),wlisMet.windSpd_Kts(B));
plot(clisMet.TIMESTAMP(C),clisMet.windSpd_Kts(C));
%legend('exrx','wlis','clis')
ylim([0,40]);

%%
dl = datetime(2022,1,1,0,0,0):15/60/24:datetime(2025,1,1,0,0,0);

% find the unique values in exrx
[unique_times, ~, idx] = unique(exrxMet.TIMESTAMP);  % idx gives the group index for each key
num_keys = numel(unique_times);
avg_wspd_exrx = accumarray(idx(:), exrxMet.windSpd_Kts(:), [], @mean);

i_exrx = interp1(unique_times,avg_wspd_exrx,dl);

[unique_times, ~, idx] = unique(wlisMet.TIMESTAMP);  % idx gives the group index for each key
num_keys = numel(unique_times);
avg_wspd_wlis = accumarray(idx(:), wlisMet.windSpd_Kts(:), [], @mean);

i_wlis = interp1(unique_times,avg_wspd_wlis,dl);

[unique_times, ~, idx] = unique(clisMet.TIMESTAMP);  % idx gives the group index for each key
num_keys = numel(unique_times);
avg_wspd_clis = accumarray(idx(:), clisMet.windSpd_Kts(:), [], @mean);

i_clis = interp1(unique_times,avg_wspd_clis,dl);

%%
x = i_exrx;  % data
win = 21;    % size of the moving window (should be odd)

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
plot(dl,i_exrx);
plot(dl,x_smoothed);

%%


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

