% fill in the bottom water CTD data
LIS = LISAug23_CH4N2O_CTD;

load LISAug23_deep;
LISdeep = LISAug23_deep;

castlist = unique(LIS.CastNum);

%for i = 1:length(castlist)
for i = 1:length(castlist)
       
    A = find(LIS.CastNum == castlist(i) & isnan(LIS.NiskinNum));
   % B = isnan(LIS.NiskinNum);
   % C = intersect(A,B)
   B = find(LIS.CastNum == castlist(i));
   C = find(max(LIS.Depth(B)));
   LIS.NiskinNum(A) = 0;
   LIS.S(A) = LIS.S(B(C));
   LIS.T(A) = LIS.T(B(C));
   LIS.PTemp(A) = LIS.PTemp(B(C));
   LIS.Dens(A) = LIS.Dens(B(C));
   LIS.PDen(A) = LIS.PDen(B(C));
   LIS.O2_umolkg(A) = LIS.O2_umolkg(B(C));   
   LIS.O2_umolkg_A(A) = LIS.O2_umolkg_A(B(C));  
   LIS.Chl(A) = LIS.Chl(B(C)); 
   LIS.pH(A) = LIS.pH(B(C));  
   LIS.Lat(A) = LIS.Lat(B(C));  
   LIS.Lon(A) = LIS.Lon(B(C));  
   LIS.datetime(A) = LIS.datetime(B(C));  
   LIS.Station(A) = LIS.Station(B(C));  

   X = find(LISdeep.Cast==castlist(i));
   LIS.P(A) = LISdeep.Pmax(X)-0.2; % removing 0.2 m due to height of niskin bottle, to give roughly center
   LIS.Depth(A) = sw_dpth(LIS.P(A),41);

end;

%%

LISAug23_CH4N2O_CTD = LIS;

save LISAug23_CH4N2O_CTD.mat LISAug23_CH4N2O_CTD;