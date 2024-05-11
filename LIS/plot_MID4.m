load LISAug23_CH4N2O_CTD.mat;
LIS = LISAug23_CH4N2O_CTD;

load LISOct23_CH4N2O_CTD.mat;
LISO = LISOct23_CH4N2O_CTD;

%
% we need to make a grid that is evenly spaced so that all the casts are
% interpolated onto the same spacing

dl = [0:0.5:20]'; % depth spacing for interpolation
x = [1 3 5 7 9 11 13 14]; % temporary x values for interpolation

n_stn = numel(x);

dl_grid = repmat(dl,1,n_stn);

x_grid = repmat(x,length(dl),1);

CH4i = nan(numel(dl),n_stn); % make a blank grid for storing CH4
N2Oi = nan(numel(dl),n_stn); % make a blank grid for storing CH4


stn = 'MID4-cast01';
A = find(LIS.Station==stn);
CH4i(:,1) = interp1(LIS.Depth(A),LIS.mean_CH4_nM(A),dl);
N2Oi(:,1) = interp1(LIS.Depth(A),LIS.mean_N2O_nM(A),dl);


stn = 'MID4-cast03';
A = find(LIS.Station==stn);
CH4i(:,2) = interp1(LIS.Depth(A),LIS.mean_CH4_nM(A),dl);
N2Oi(:,2) = interp1(LIS.Depth(A),LIS.mean_N2O_nM(A),dl);


stn = 'MID4-cast05';
A = find(LIS.Station==stn);
CH4i(:,3) = interp1(LIS.Depth(A),LIS.mean_CH4_nM(A),dl);
N2Oi(:,3) = interp1(LIS.Depth(A),LIS.mean_N2O_nM(A),dl);

stn = 'MID4-cast07';
A = find(LIS.Station==stn);
CH4i(:,4) = interp1(LIS.Depth(A),LIS.mean_CH4_nM(A),dl);
N2Oi(:,4) = interp1(LIS.Depth(A),LIS.mean_N2O_nM(A),dl);


stn = 'MID4-cast09';
A = find(LIS.Station==stn);
CH4i(:,5) = interp1(LIS.Depth(A),LIS.mean_CH4_nM(A),dl);
N2Oi(:,5) = interp1(LIS.Depth(A),LIS.mean_N2O_nM(A),dl);

stn = 'MID4-cast11';
A = find(LIS.Station==stn);
CH4i(:,6) = interp1(LIS.Depth(A),LIS.mean_CH4_nM(A),dl);
N2Oi(:,6) = interp1(LIS.Depth(A),LIS.mean_N2O_nM(A),dl);

stn = 'MID4-cast13';
A = find(LIS.Station==stn);
CH4i(:,7) = interp1(LIS.Depth(A),LIS.mean_CH4_nM(A),dl);
N2Oi(:,7) = interp1(LIS.Depth(A),LIS.mean_N2O_nM(A),dl);

stn = 'MID4-cast14';
A = find(LIS.Station==stn);
CH4i(:,8) = interp1(LIS.Depth(A),LIS.mean_CH4_nM(A),dl);
N2Oi(:,8) = interp1(LIS.Depth(A),LIS.mean_N2O_nM(A),dl);


figure(2)
clf; hold on;
box on;
set(gca,'tickdir','out');
clevel = [40:1:80];
C = contourf(x_grid,dl_grid,CH4i,clevel,'edgecolor','none');
xlabel('cast number');
ylabel('Depth (m)');
c = colorbar;
c.Label.String = 'CH_4 (nM)';
c.Label.FontSize = 16;

axis ij;

figure(3)
clf; hold on;
box on;
set(gca,'tickdir','out');
clevel = [8:0.1:12];
C = contourf(x_grid,dl_grid,N2Oi,clevel,'edgecolor','none');
xlabel('cast number');
ylabel('Depth (m)');
c = colorbar;
c.Label.String = 'N_2O (nM)';
c.Label.FontSize = 16;

axis ij;


