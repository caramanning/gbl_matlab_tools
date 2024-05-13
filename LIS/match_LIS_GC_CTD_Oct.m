load LISOct23_CH4N2O.mat
GC = LISOct23_CH4N2O;

%% rename all the cast numbers in the file
GC.CastNum = nan.*GC.mean_CH4_nM;

A=find(GC.Cast=='Ca2');
GC.CastNum(A) = 2;

A=find(GC.Cast=='Ca3');
GC.CastNum(A) = 3;

A=find(GC.Cast=='Ca5');
GC.CastNum(A) = 5;

A=find(GC.Cast=='Ca7');
GC.CastNum(A) = 7;

A=find(GC.Cast=='Ca9');
GC.CastNum(A) = 9;

A=find(GC.Cast=='Ca11');
GC.CastNum(A) = 11;

A=find(GC.Cast=='Ca13');
GC.CastNum(A) = 13;

A=find(GC.Cast=='Ca15');
GC.CastNum(A) = 15;

A=find(GC.Cast=='Ca16');
GC.CastNum(A) = 16;

A=find(GC.Cast=='Ca17');
GC.CastNum(A) = 17;

A=find(GC.Cast=='Ca18');
GC.CastNum(A) = 18;

A=find(GC.Cast=='Ca19');
GC.CastNum(A) = 19;

A=find(GC.Cast=='Ca20');
GC.CastNum(A) = 20;

A=find(GC.Cast=='Ca21');
GC.CastNum(A) = 21;

A=find(GC.Cast=='Ca22');
GC.CastNum(A) = 22;

%% rename all the numeric niskin numbers in the file
GC.NiskinNum = nan.*GC.mean_CH4_nM;


A=find(GC.Niskin=='Ni1');
GC.NiskinNum(A) = 1;

A=find(GC.Niskin=='Ni2');
GC.NiskinNum(A) = 2;

A=find(GC.Niskin=='Ni3');
GC.NiskinNum(A) = 3;

A=find(GC.Niskin=='Ni4');
GC.NiskinNum(A) = 4;

A=find(GC.Niskin=='Ni5');
GC.NiskinNum(A) = 5;

A=find(GC.Niskin=='Ni6');
GC.NiskinNum(A) = 6;

A=find(GC.Niskin=='Ni7');
GC.NiskinNum(A) = 7;

A=find(GC.Niskin=='Ni8');
GC.NiskinNum(A) = 8;

A=find(GC.Niskin=='Ni9');
GC.NiskinNum(A) = 9;

A=find(GC.Niskin=='Ni10');
GC.NiskinNum(A) = 10;

A=find(GC.Niskin=='Ni11');
GC.NiskinNum(A) = 11;

A=find(GC.Niskin=='Ni12');
GC.NiskinNum(A) = 12;


LISAug23_CH4N2O_CTD = GC;

save LISOct23_CH4N2O_CTD.mat LISAug23_CH4N2O_CTD ;

%%
load LISOct23btl.mat;
btl = LISOct23btl;

LIS = GC; % shorten name temporarily
LIS.Depth = nan.*LIS.CastNum;
LIS.P = nan.*LIS.CastNum;
LIS.S = nan.*LIS.CastNum;
LIS.T = nan.*LIS.CastNum;
LIS.PTemp = nan.*LIS.CastNum;
LIS.Dens = nan.*LIS.CastNum;
LIS.PDen = nan.*LIS.CastNum;
LIS.O2_umolkg = nan.*LIS.CastNum;
LIS.O2_umolkg_A = nan.*LIS.CastNum;
LIS.Chl = nan.*LIS.CastNum;
LIS.pH = nan.*LIS.CastNum;
LIS.Lat = nan.*LIS.CastNum;
LIS.Lon = nan.*LIS.CastNum;
%LIS.datetime = nan.*LIS.CastNum;

for i = 1:length(GC.Cast)
%for i = 1  
    A = find(btl.Cast==LIS.CastNum(i));
    B = find(btl.Niskin==LIS.NiskinNum(i));
    C = intersect(A,B);
    if ~isempty(C)
    LIS.Depth(i) = btl.Depth(C); 
    LIS.P(i) = btl.P(C);
    LIS.S(i) = btl.S(C);
    LIS.T(i) = btl.T(C);
    LIS.PTemp(i) = btl.PTemp(C);
    LIS.Dens(i) = btl.Dens(C);
    LIS.PDen(i) = btl.PDen(C);
    LIS.O2_umolkg(i) = btl.O2_umolkg(C);
    %LIS.O2_umolkg_A(i) = btl.O2_umolkg_A(C);
    LIS.Chl(i) = btl.Chl(C);
    LIS.pH(i) = btl.pH(C);
    LIS.Lat(i) = btl.Lat(C);
    LIS.Lon(i) = btl.Lon(C);
    LIS.datetime(i) = btl.datetime(C);
    LIS.Station(i) = btl.Station(C);
    end;
end

LISOct23_CH4N2O_CTD = LIS;

save LISOct23_CH4N2O_CTD.mat LISOct23_CH4N2O_CTD

%%





