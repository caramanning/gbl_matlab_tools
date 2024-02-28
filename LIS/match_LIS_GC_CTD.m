load LISAug23_CH4N2O.mat

load LISAug23btl.mat;
%% rename all the cast numbers in the file
LISAug23_CH4N2O.CastNum = nan.*LISAug23_CH4N2O.mean_CH4_nM;

A=find(LISAug23_CH4N2O.Cast=='Ca2');
LISAug23_CH4N2O.CastNum(A) = 2;

A=find(LISAug23_CH4N2O.Cast=='Ca3');
LISAug23_CH4N2O.CastNum(A) = 3;

A=find(LISAug23_CH4N2O.Cast=='Ca5');
LISAug23_CH4N2O.CastNum(A) = 5;

A=find(LISAug23_CH4N2O.Cast=='Ca7');
LISAug23_CH4N2O.CastNum(A) = 7;

A=find(LISAug23_CH4N2O.Cast=='Ca9');
LISAug23_CH4N2O.CastNum(A) = 9;

A=find(LISAug23_CH4N2O.Cast=='Ca11');
LISAug23_CH4N2O.CastNum(A) = 11;

A=find(LISAug23_CH4N2O.Cast=='Ca13');
LISAug23_CH4N2O.CastNum(A) = 13;

A=find(LISAug23_CH4N2O.Cast=='Ca15');
LISAug23_CH4N2O.CastNum(A) = 15;

A=find(LISAug23_CH4N2O.Cast=='Ca16');
LISAug23_CH4N2O.CastNum(A) = 16;

A=find(LISAug23_CH4N2O.Cast=='Ca17');
LISAug23_CH4N2O.CastNum(A) = 17;

A=find(LISAug23_CH4N2O.Cast=='Ca18');
LISAug23_CH4N2O.CastNum(A) = 18;

A=find(LISAug23_CH4N2O.Cast=='Ca19');
LISAug23_CH4N2O.CastNum(A) = 19;

A=find(LISAug23_CH4N2O.Cast=='Ca20');
LISAug23_CH4N2O.CastNum(A) = 20;

A=find(LISAug23_CH4N2O.Cast=='Ca21');
LISAug23_CH4N2O.CastNum(A) = 21;

A=find(LISAug23_CH4N2O.Cast=='Ca22');
LISAug23_CH4N2O.CastNum(A) = 22;

%% rename all the numeric niskin numbers in the file
LISAug23_CH4N2O.NiskinNum = nan.*LISAug23_CH4N2O.mean_CH4_nM;


A=find(LISAug23_CH4N2O.Niskin=='Ni1');
LISAug23_CH4N2O.NiskinNum(A) = 1;

A=find(LISAug23_CH4N2O.Niskin=='Ni2');
LISAug23_CH4N2O.NiskinNum(A) = 2;

A=find(LISAug23_CH4N2O.Niskin=='Ni3');
LISAug23_CH4N2O.NiskinNum(A) = 3;

A=find(LISAug23_CH4N2O.Niskin=='Ni4');
LISAug23_CH4N2O.NiskinNum(A) = 4;

A=find(LISAug23_CH4N2O.Niskin=='Ni5');
LISAug23_CH4N2O.NiskinNum(A) = 5;

A=find(LISAug23_CH4N2O.Niskin=='Ni6');
LISAug23_CH4N2O.NiskinNum(A) = 6;

A=find(LISAug23_CH4N2O.Niskin=='Ni7');
LISAug23_CH4N2O.NiskinNum(A) = 7;

A=find(LISAug23_CH4N2O.Niskin=='Ni8');
LISAug23_CH4N2O.NiskinNum(A) = 8;

A=find(LISAug23_CH4N2O.Niskin=='Ni9');
LISAug23_CH4N2O.NiskinNum(A) = 9;

A=find(LISAug23_CH4N2O.Niskin=='Ni10');
LISAug23_CH4N2O.NiskinNum(A) = 10;

A=find(LISAug23_CH4N2O.Niskin=='Ni11');
LISAug23_CH4N2O.NiskinNum(A) = 11;

A=find(LISAug23_CH4N2O.Niskin=='Ni12');
LISAug23_CH4N2O.NiskinNum(A) = 12;


LISAug23_CH4N2O_CTD = LISAug23_CH4N2O;

save LISAug23_CH4N2O_CTD.mat LISAug23_CH4N2O_CTD;

%%
load LISAug23btl.mat;
LIS = LISAug23_CH4N2O_CTD; % shorten name temporarily
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

for i = 1:length(LISAug23_CH4N2O_CTD.Cast)
%for i = 1  
    A = find(LISAug23btl.Cast==LIS.CastNum(i));
    B = find(LISAug23btl.Niskin==LIS.NiskinNum(i));
    C = intersect(A,B);
    if ~isempty(C)
    LIS.Depth(i) = LISAug23btl.Depth(C); 
    LIS.P(i) = LISAug23btl.P(C);
    LIS.S(i) = LISAug23btl.S(C);
    LIS.T(i) = LISAug23btl.T(C);
    LIS.PTemp(i) = LISAug23btl.PTemp(C);
    LIS.Dens(i) = LISAug23btl.Dens(C);
    LIS.PDen(i) = LISAug23btl.PDen(C);
    LIS.O2_umolkg(i) = LISAug23btl.O2_umolkg(C);
    LIS.O2_umolkg_A(i) = LISAug23btl.O2_umolkg_A(C);
    LIS.Chl(i) = LISAug23btl.Chl(C);
    LIS.pH(i) = LISAug23btl.pH(C);
    LIS.Lat(i) = LISAug23btl.Lat(C);
    LIS.Lon(i) = LISAug23btl.Lon(C);
    LIS.datetime(i) = LISAug23btl.datetime(C);
    LIS.Station(i) = LISAug23btl.Station(C);
    end;
end

LISAug23_CH4N2O_CTD = LIS;

save LISAug23_CH4N2O_CTD.mat LISAug23_CH4N2O_CTD

%%





