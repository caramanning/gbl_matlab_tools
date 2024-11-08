load LISAug23_CH4N2O_CTD.mat
LIS = LISAug23_CH4N2O_CTD;


yl = [0 20]; % y-axis limits
xl = [0 500]; % x-axis limits
figure(1)
clf; 

subplot(1,2,1)
hold on;
stn="EXCR-cast01";
A=find(LIS.Station==stn & isnumeric(LIS.Lon));

errorbar(LIS.mean_CH4_nM(A),LIS.Depth(A),LIS.std_CH4_nM(A),'horizontal','o-b','linewidth',2,'markerfacecolor','b');
xlabel('CH_4 (nM)');
ylabel('Depth (m)');
title(stn);
axis ij;
set(gca,'xlim',xl);
set(gca,'ylim',yl);

subplot(1,2,2)
hold on;
stn="EXCR-cast02";
A=find(LIS.Station==stn & isnumeric(LIS.Lon));

errorbar(LIS.mean_CH4_nM(A),LIS.Depth(A),LIS.std_CH4_nM(A),'horizontal','o-b','linewidth',2,'markerfacecolor','b');
xlabel('CH_4 (nM)');
ylabel('Depth (m)');
title(stn);
axis ij;
set(gca,'xlim',xl);
set(gca,'ylim',yl);


%%
dl = [0:0.5:20]'; % depths for interpolation
x = [1:8]; % temporary x values for interpolation
n_stn = numel(x);

CH4i = nan(length(dl),n_stn);

dl_grid = repmat(dl,1,n_stn);

x_grid = repmat(x,length(dl),1);

stn="MID4-cast01";
A=find(LIS.Station==stn & isnumeric(LIS.Lon));
CH4i(:,1) = interp1(LIS.Depth(A),LIS.mean_CH4_nM(A),dl);

stn="MID4-cast03";
A=find(LIS.Station==stn & isnumeric(LIS.Lon));
CH4i(:,2) = interp1(LIS.Depth(A),LIS.mean_CH4_nM(A),dl);

stn="MID4-cast05";
A=find(LIS.Station==stn & isnumeric(LIS.Lon));
CH4i(:,3) = interp1(LIS.Depth(A),LIS.mean_CH4_nM(A),dl);

stn="MID4-cast07";
A=find(LIS.Station==stn & isnumeric(LIS.Lon));
CH4i(:,4) = interp1(LIS.Depth(A),LIS.mean_CH4_nM(A),dl);

stn="MID4-cast09";
A=find(LIS.Station==stn & isnumeric(LIS.Lon));
CH4i(:,5) = interp1(LIS.Depth(A),LIS.mean_CH4_nM(A),dl);

stn="MID4-cast11";
A=find(LIS.Station==stn & isnumeric(LIS.Lon));
CH4i(:,6) = interp1(LIS.Depth(A),LIS.mean_CH4_nM(A),dl);

stn="MID4-cast13";
A=find(LIS.Station==stn & isnumeric(LIS.Lon));
CH4i(:,7) = interp1(LIS.Depth(A),LIS.mean_CH4_nM(A),dl);

stn="MID4-cast09";
A=find(LIS.Station==stn & isnumeric(LIS.Lon));
CH4i(:,8) = interp1(LIS.Depth(A),LIS.mean_CH4_nM(A),dl);



figure(2)
clf; hold on;
contourf(x_grid,dl_grid,CH4i,'edgecolor','none')
axis ij;
colorbar;

%%

N2Oi = nan(length(dl),n_stn);

dl_grid = repmat(dl,1,n_stn);

x_grid = repmat(x,length(dl),1);

stn="MID4-cast01";
A=find(LIS.Station==stn & isnumeric(LIS.Lon));
N2Oi(:,1) = interp1(LIS.Depth(A),LIS.mean_N2O_nM(A),dl);

stn="MID4-cast03";
A=find(LIS.Station==stn & isnumeric(LIS.Lon));
N2Oi(:,2) = interp1(LIS.Depth(A),LIS.mean_N2O_nM(A),dl);

stn="MID4-cast05";
A=find(LIS.Station==stn & isnumeric(LIS.Lon));
N2Oi(:,3) = interp1(LIS.Depth(A),LIS.mean_N2O_nM(A),dl);

stn="MID4-cast07";
A=find(LIS.Station==stn & isnumeric(LIS.Lon));
N2Oi(:,4) = interp1(LIS.Depth(A),LIS.mean_N2O_nM(A),dl);

stn="MID4-cast09";
A=find(LIS.Station==stn & isnumeric(LIS.Lon));
N2Oi(:,5) = interp1(LIS.Depth(A),LIS.mean_N2O_nM(A),dl);

stn="MID4-cast11";
A=find(LIS.Station==stn & isnumeric(LIS.Lon));
N2Oi(:,6) = interp1(LIS.Depth(A),LIS.mean_N2O_nM(A),dl);

stn="MID4-cast13";
A=find(LIS.Station==stn & isnumeric(LIS.Lon));
N2Oi(:,7) = interp1(LIS.Depth(A),LIS.mean_N2O_nM(A),dl);

stn="MID4-cast09";
A=find(LIS.Station==stn & isnumeric(LIS.Lon));
N2Oi(:,8) = interp1(LIS.Depth(A),LIS.mean_N2O_nM(A),dl);



figure(3)
clf; hold on;
contourf(x_grid,dl_grid,N2Oi,'edgecolor','none')
axis ij;
colorbar;
