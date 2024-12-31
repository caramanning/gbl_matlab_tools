
%% Import all the data for each cast

load WLIScast01_btl.mat;
load WLIScast01_bl.mat;
WLIScast01_btl.Niskin = nan.*WLIScast01_btl.FireSeq;
WLIScast01_btl.Cast = nan.*WLIScast01_btl.FireSeq;

for i = 1:numel(WLIScast01_bl.FireSeq)
    ifs = find(WLIScast01_bl.FireSeq==WLIScast01_btl.FireSeq(i));
    WLIScast01_btl.Niskin(i) = WLIScast01_bl.Niskin(ifs);
    WLIScast01_btl.Cast(i) = 1;
end;
%%
load MID4cast01_btl.mat;
load MID4cast01_bl.mat;
MID4cast01_btl.Niskin = nan.*MID4cast01_btl.FireSeq;
MID4cast01_btl.Cast = nan.*MID4cast01_btl.FireSeq;

for i = 1:numel(MID4cast01_bl.FireSeq)
    ifs = find(MID4cast01_bl.FireSeq==MID4cast01_btl.FireSeq(i));
    MID4cast01_btl.Niskin(i) = MID4cast01_bl.Niskin(ifs);
    MID4cast01_btl.Cast(i) = 2;
end;


load EXRXcast01_btl.mat;
load EXRXcast01_bl.mat;
EXRXcast01_btl.Niskin = nan.*EXRXcast01_btl.FireSeq;
EXRXcast01_btl.Cast = nan.*EXRXcast01_btl.FireSeq;

for i = 1:numel(EXRXcast01_bl.FireSeq)
    ifs = find(EXRXcast01_bl.FireSeq==EXRXcast01_btl.FireSeq(i));
    EXRXcast01_btl.Niskin(i) = EXRXcast01_bl.Niskin(ifs);
    EXRXcast01_btl.Cast(i) = 3;
end;



load MID4cast02_btl.mat;
load MID4cast02_bl.mat;
MID4cast02_btl.Niskin = nan.*MID4cast02_btl.FireSeq;
MID4cast02_btl.Cast = nan.*MID4cast02_btl.FireSeq;

for i = 1:numel(MID4cast02_bl.FireSeq)
    ifs = find(MID4cast02_bl.FireSeq==MID4cast02_btl.FireSeq(i));
    MID4cast02_btl.Niskin(i) = MID4cast02_bl.Niskin(ifs);
    MID4cast02_btl.Cast(i) = 4;
end;
%%
load MID4cast03_btl.mat;
load MID4cast03_bl.mat;
MID4cast03_btl.Niskin = nan.*MID4cast03_btl.FireSeq;
MID4cast03_btl.Cast = nan.*MID4cast03_btl.FireSeq;

for i = 1:numel(MID4cast03_bl.FireSeq)
    ifs = find(MID4cast03_bl.FireSeq==MID4cast03_btl.FireSeq(i));
    MID4cast03_btl.Niskin(i) = MID4cast03_bl.Niskin(ifs);
    MID4cast03_btl.Cast(i) = 5;
end;

load MID4cast04_btl.mat;
load MID4cast04_bl.mat;
MID4cast04_btl.Niskin = nan.*MID4cast04_btl.FireSeq;
MID4cast04_btl.Cast = nan.*MID4cast04_btl.FireSeq;

for i = 1:numel(MID4cast04_bl.FireSeq)
    ifs = find(MID4cast04_bl.FireSeq==MID4cast04_btl.FireSeq(i));
    MID4cast04_btl.Niskin(i) = MID4cast04_bl.Niskin(ifs);
    MID4cast04_btl.Cast(i) = 6;
end;

load MID4cast05_btl.mat;
load MID4cast05_bl.mat;
MID4cast05_btl.Niskin = nan.*MID4cast05_btl.FireSeq;
MID4cast05_btl.Cast = nan.*MID4cast05_btl.FireSeq;

for i = 1:numel(MID4cast05_bl.FireSeq)
    ifs = find(MID4cast05_bl.FireSeq==MID4cast05_btl.FireSeq(i));
    MID4cast05_btl.Niskin(i) = MID4cast05_bl.Niskin(ifs);
    MID4cast05_btl.Cast(i) = 7;
end;

load MID4cast06_btl.mat;
load MID4cast06_bl.mat;
MID4cast06_btl.Niskin = nan.*MID4cast06_btl.FireSeq;
MID4cast06_btl.Cast = nan.*MID4cast06_btl.FireSeq;

for i = 1:numel(MID4cast06_bl.FireSeq)
    ifs = find(MID4cast06_bl.FireSeq==MID4cast06_btl.FireSeq(i));
    MID4cast06_btl.Niskin(i) = MID4cast06_bl.Niskin(ifs);
    MID4cast06_btl.Cast(i) = 8;
end;

load MID4cast07_btl.mat;
load MID4cast07_bl.mat;
MID4cast07_btl.Niskin = nan.*MID4cast07_btl.FireSeq;
MID4cast07_btl.Cast = nan.*MID4cast07_btl.FireSeq;

for i = 1:numel(MID4cast07_bl.FireSeq)
    ifs = find(MID4cast07_bl.FireSeq==MID4cast07_btl.FireSeq(i));
    MID4cast07_btl.Niskin(i) = MID4cast07_bl.Niskin(ifs);
    MID4cast07_btl.Cast(i) = 9;
end;

load MID4cast08_btl.mat;
load MID4cast08_bl.mat;
MID4cast08_btl.Niskin = nan.*MID4cast08_btl.FireSeq;
MID4cast08_btl.Cast = nan.*MID4cast08_btl.FireSeq;

for i = 1:numel(MID4cast08_bl.FireSeq)
    ifs = find(MID4cast08_bl.FireSeq==MID4cast08_btl.FireSeq(i));
    MID4cast08_btl.Niskin(i) = MID4cast08_bl.Niskin(ifs);
    MID4cast08_btl.Cast(i) = 10;
end;

load MID4cast09_btl.mat;
load MID4cast09_bl.mat;
MID4cast09_btl.Niskin = nan.*MID4cast09_btl.FireSeq;
MID4cast09_btl.Cast = nan.*MID4cast09_btl.FireSeq;

for i = 1:numel(MID4cast09_bl.FireSeq)
    ifs = find(MID4cast09_bl.FireSeq==MID4cast09_btl.FireSeq(i));
    MID4cast09_btl.Niskin(i) = MID4cast09_bl.Niskin(ifs);
    MID4cast09_btl.Cast(i) = 11;
end;

load MID4cast10_btl.mat;
load MID4cast10_bl.mat;
MID4cast10_btl.Niskin = nan.*MID4cast10_btl.FireSeq;
MID4cast10_btl.Cast = nan.*MID4cast10_btl.FireSeq;

for i = 1:numel(MID4cast10_bl.FireSeq)
    ifs = find(MID4cast10_bl.FireSeq==MID4cast10_btl.FireSeq(i));
    MID4cast10_btl.Niskin(i) = MID4cast10_bl.Niskin(ifs);
    MID4cast10_btl.Cast(i) = 12;
end;


load MID4cast11_btl.mat;
load MID4cast11_bl.mat;
MID4cast11_btl.Niskin = nan.*MID4cast11_btl.FireSeq;
MID4cast11_btl.Cast = nan.*MID4cast11_btl.FireSeq;

for i = 1:numel(MID4cast11_bl.FireSeq)
    ifs = find(MID4cast11_bl.FireSeq==MID4cast11_btl.FireSeq(i));
    MID4cast11_btl.Niskin(i) = MID4cast11_bl.Niskin(ifs);
    MID4cast11_btl.Cast(i) = 13;
end;

load MID4cast12_btl.mat;
load MID4cast12_bl.mat;
MID4cast12_btl.Niskin = nan.*MID4cast12_btl.FireSeq;
MID4cast12_btl.Cast = nan.*MID4cast12_btl.FireSeq;

for i = 1:numel(MID4cast12_bl.FireSeq)
    ifs = find(MID4cast12_bl.FireSeq==MID4cast12_btl.FireSeq(i));
    MID4cast12_btl.Niskin(i) = MID4cast12_bl.Niskin(ifs);
    MID4cast12_btl.Cast(i) = 14;
end;

load MID4cast13_btl.mat;
load MID4cast13_bl.mat;
MID4cast13_btl.Niskin = nan.*MID4cast13_btl.FireSeq;
MID4cast13_btl.Cast = nan.*MID4cast13_btl.FireSeq;

for i = 1:numel(MID4cast13_bl.FireSeq)
    ifs = find(MID4cast13_bl.FireSeq==MID4cast13_btl.FireSeq(i));
    MID4cast13_btl.Niskin(i) = MID4cast13_bl.Niskin(ifs);
    MID4cast13_btl.Cast(i) = 15;
end;


load EXR1cast01_btl.mat;
load EXR1cast01_bl.mat;
EXR1cast01_btl.Niskin = nan.*EXR1cast01_btl.FireSeq;
EXR1cast01_btl.Cast = nan.*EXR1cast01_btl.FireSeq;

for i = 1:numel(EXR1cast01_bl.FireSeq)
    ifs = find(EXR1cast01_bl.FireSeq==EXR1cast01_btl.FireSeq(i));
    EXR1cast01_btl.Niskin(i) = EXR1cast01_bl.Niskin(ifs);
    EXR1cast01_btl.Cast(i) = 16;
end;


load EXRXcast02_btl.mat;
load EXRXcast02_bl.mat;
EXRXcast02_btl.Niskin = nan.*EXRXcast02_btl.FireSeq;
EXRXcast02_btl.Cast = nan.*EXRXcast02_btl.FireSeq;

for i = 1:numel(EXRXcast02_bl.FireSeq)
    ifs = find(EXRXcast02_bl.FireSeq==EXRXcast02_btl.FireSeq(i));
    EXRXcast02_btl.Niskin(i) = EXRXcast02_bl.Niskin(ifs);
    EXRXcast02_btl.Cast(i) = 17;
end;


load MID3cast01_btl.mat;
load MID3cast01_bl.mat;
MID3cast01_btl.Niskin = nan.*MID3cast01_btl.FireSeq;
MID3cast01_btl.Cast = nan.*MID3cast01_btl.FireSeq;

for i = 1:numel(MID3cast01_bl.FireSeq)
    ifs = find(MID3cast01_bl.FireSeq==MID3cast01_btl.FireSeq(i));
    MID3cast01_btl.Niskin(i) = MID3cast01_bl.Niskin(ifs);
    MID3cast01_btl.Cast(i) = 18;
end;


load MID4cast14_btl.mat;
load MID4cast14_bl.mat;
MID4cast14_btl.Niskin = nan.*MID4cast14_btl.FireSeq;
MID4cast14_btl.Cast = nan.*MID4cast14_btl.FireSeq;

for i = 1:numel(MID4cast14_bl.FireSeq)
    ifs = find(MID4cast14_bl.FireSeq==MID4cast14_btl.FireSeq(i));
    MID4cast14_btl.Niskin(i) = MID4cast14_bl.Niskin(ifs);
    MID4cast14_btl.Cast(i) = 19;
end;

%%
load MID5cast01_btl.mat;
load MID5cast01_bl.mat;
MID5cast01_btl.Niskin = nan.*MID5cast01_btl.FireSeq;
MID5cast01_btl.Cast = nan.*MID5cast01_btl.FireSeq;

for i = 1:numel(MID5cast01_bl.FireSeq)
    ifs = find(MID5cast01_bl.FireSeq==MID5cast01_btl.FireSeq(i));
    MID5cast01_btl.Niskin(i) = MID5cast01_bl.Niskin(ifs);
    MID5cast01_btl.Cast(i) = 20;
end;

load WLIScast02_btl.mat;
load WLIScast02_bl.mat;
WLIScast02_btl.Niskin = nan.*WLIScast02_btl.FireSeq;
WLIScast02_btl.Cast = nan.*WLIScast02_btl.FireSeq;

for i = 1:numel(WLIScast02_bl.FireSeq)
    ifs = find(WLIScast02_bl.FireSeq==WLIScast02_btl.FireSeq(i));
    WLIScast02_btl.Niskin(i) = WLIScast02_bl.Niskin(ifs);
    WLIScast02_btl.Cast(i) = 21;
end;


load WLI6cast01_btl.mat;
load WLI6cast01_bl.mat;
WLI6cast01_btl.Niskin = nan.*WLI6cast01_btl.FireSeq;
WLI6cast01_btl.Cast = nan.*WLI6cast01_btl.FireSeq;

for i = 1:numel(WLI6cast01_bl.FireSeq)
    ifs = find(WLI6cast01_bl.FireSeq==WLI6cast01_btl.FireSeq(i));
    WLI6cast01_btl.Niskin(i) = WLI6cast01_bl.Niskin(ifs);
    WLI6cast01_btl.Cast(i) = 22;
end;


load WLI7cast01_btl.mat;
load WLI7cast01_bl.mat;
WLI7cast01_btl.Niskin = nan.*WLI7cast01_btl.FireSeq;
WLI7cast01_btl.Cast = nan.*WLI7cast01_btl.FireSeq;
 
for i = 1:numel(WLI7cast01_bl.FireSeq)
     ifs = find(WLI7cast01_bl.FireSeq==WLI7cast01_btl.FireSeq(i));
     WLI7cast01_btl.Niskin(i) = WLI7cast01_bl.Niskin(ifs);
     WLI7cast01_btl.Cast(i) = 23;
 end;


load ARTGcast01_btl.mat;
load ARTGcast01_bl.mat;
ARTGcast01_btl.Niskin = nan.*ARTGcast01_btl.FireSeq;
ARTGcast01_btl.Cast = nan.*ARTGcast01_btl.FireSeq;
 
for i = 1:numel(ARTGcast01_bl.FireSeq)
     ifs = find(ARTGcast01_bl.FireSeq==ARTGcast01_btl.FireSeq(i));
     ARTGcast01_btl.Niskin(i) = ARTGcast01_bl.Niskin(ifs);
     ARTGcast01_btl.Cast(i) = 24;
 end;


load CLIScast01_btl.mat;
load CLIScast01_bl.mat;
CLIScast01_btl.Niskin = nan.*CLIScast01_btl.FireSeq;
CLIScast01_btl.Cast = nan.*CLIScast01_btl.FireSeq;
 
for i = 1:numel(CLIScast01_bl.FireSeq)
     ifs = find(CLIScast01_bl.FireSeq==CLIScast01_btl.FireSeq(i));
     CLIScast01_btl.Niskin(i) = CLIScast01_bl.Niskin(ifs);
     CLIScast01_btl.Cast(i) = 25;
 end;



%%


Cast = [WLIScast01_btl.Cast; MID4cast01_btl.Cast; EXRXcast01_btl.Cast; MID4cast02_btl.Cast; ...
    MID4cast03_btl.Cast; MID4cast04_btl.Cast; MID4cast05_btl.Cast; MID4cast06_btl.Cast; ...
    MID4cast07_btl.Cast; MID4cast08_btl.Cast; MID4cast09_btl.Cast; MID4cast10_btl.Cast; ...
    MID4cast11_btl.Cast; MID4cast12_btl.Cast; MID4cast13_btl.Cast; EXR1cast01_btl.Cast; ...
    EXRXcast02_btl.Cast; MID3cast01_btl.Cast; MID4cast14_btl.Cast; MID5cast01_btl.Cast; ...
    WLIScast02_btl.Cast; WLI6cast01_btl.Cast; WLI7cast01_btl.Cast; ARTGcast01_btl.Cast; CLIScast01_btl.Cast; ...
    ];

Niskin = [WLIScast01_btl.Niskin; MID4cast01_btl.Niskin; EXRXcast01_btl.Niskin; MID4cast02_btl.Niskin; ...
    MID4cast03_btl.Niskin; MID4cast04_btl.Niskin; MID4cast05_btl.Niskin; MID4cast06_btl.Niskin; ...
    MID4cast07_btl.Niskin; MID4cast08_btl.Niskin; MID4cast09_btl.Niskin; MID4cast10_btl.Niskin; ...
    MID4cast11_btl.Niskin; MID4cast12_btl.Niskin; MID4cast13_btl.Niskin; EXR1cast01_btl.Niskin; ...
    EXRXcast02_btl.Niskin; MID3cast01_btl.Niskin; MID4cast14_btl.Niskin; MID5cast01_btl.Niskin; ...
    WLIScast02_btl.Niskin; WLI6cast01_btl.Niskin; WLI7cast01_btl.Niskin; ARTGcast01_btl.Niskin; CLIScast01_btl.Niskin; ...
    ];


%%
P = [WLIScast01_btl.PrDM; MID4cast01_btl.PrDM; EXRXcast01_btl.PrDM; MID4cast02_btl.PrDM; ...
    MID4cast03_btl.PrDM; MID4cast04_btl.PrDM; MID4cast05_btl.PrDM; MID4cast06_btl.PrDM; ...
    MID4cast07_btl.PrDM; MID4cast08_btl.PrDM; MID4cast09_btl.PrDM; MID4cast10_btl.PrDM; ...
    MID4cast11_btl.PrDM; MID4cast12_btl.PrDM; MID4cast13_btl.PrDM; EXR1cast01_btl.PrDM; ...
    EXRXcast02_btl.PrDM; MID3cast01_btl.PrDM; MID4cast14_btl.PrDM; MID5cast01_btl.PrDM; ...
    WLIScast02_btl.PrDM; WLI6cast01_btl.PrDM; WLI7cast01_btl.PrDM; ARTGcast01_btl.PrDM; CLIScast01_btl.PrDM; ...
    ];

S = [WLIScast01_btl.Sal00; MID4cast01_btl.Sal00; EXRXcast01_btl.Sal00; MID4cast02_btl.Sal00; ...
    MID4cast03_btl.Sal00; MID4cast04_btl.Sal00; MID4cast05_btl.Sal00; MID4cast06_btl.Sal00; ...
    MID4cast07_btl.Sal00; MID4cast08_btl.Sal00; MID4cast09_btl.Sal00; MID4cast10_btl.Sal00; ...
    MID4cast11_btl.Sal00; MID4cast12_btl.Sal00; MID4cast13_btl.Sal00; EXR1cast01_btl.Sal00; ...
    EXRXcast02_btl.Sal00; MID3cast01_btl.Sal00; MID4cast14_btl.Sal00; MID5cast01_btl.Sal00; ...
    WLIScast02_btl.Sal00; WLI6cast01_btl.Sal00; WLI7cast01_btl.Sal00; ARTGcast01_btl.Sal00; CLIScast01_btl.Sal00; ...
    ];

T = [WLIScast01_btl.T090C; MID4cast01_btl.T090C; EXRXcast01_btl.T090C; MID4cast02_btl.T090C; ...
    MID4cast03_btl.T090C; MID4cast04_btl.T090C; MID4cast05_btl.T090C; MID4cast06_btl.T090C; ...
    MID4cast07_btl.T090C; MID4cast08_btl.T090C; MID4cast09_btl.T090C; MID4cast10_btl.T090C; ...
    MID4cast11_btl.T090C; MID4cast12_btl.T090C; MID4cast13_btl.T090C; EXR1cast01_btl.T090C; ...
    EXRXcast02_btl.T090C; MID3cast01_btl.T090C; MID4cast14_btl.T090C; MID5cast01_btl.T090C; ...
    WLIScast02_btl.T090C; WLI6cast01_btl.T090C; WLI7cast01_btl.T090C; ARTGcast01_btl.T090C; CLIScast01_btl.T090C; ...
    ];

%%
Depth = sw_dpth(P,41);

PTemp = sw_ptmp(S,T,P,0);


Dens = sw_dens(S,T,P);

PDen = sw_pden(S,T,P,0) - 1000;

O2_umolL = [WLIScast01_btl.O2_umolL; MID4cast01_btl.O2_umolL; EXRXcast01_btl.O2_umolL; MID4cast02_btl.O2_umolL; ...
    MID4cast03_btl.O2_umolL; MID4cast04_btl.O2_umolL; MID4cast05_btl.O2_umolL; MID4cast06_btl.O2_umolL; ...
    MID4cast07_btl.O2_umolL; MID4cast08_btl.O2_umolL; MID4cast09_btl.O2_umolL; MID4cast10_btl.O2_umolL; ...
    MID4cast11_btl.O2_umolL; MID4cast12_btl.O2_umolL; MID4cast13_btl.O2_umolL; EXR1cast01_btl.O2_umolL; ...
    EXRXcast02_btl.O2_umolL; MID3cast01_btl.O2_umolL; MID4cast14_btl.O2_umolL; MID5cast01_btl.O2_umolL; ...
    WLIScast02_btl.O2_umolL; WLI6cast01_btl.O2_umolL; WLI7cast01_btl.O2_umolL; ARTGcast01_btl.O2_umolL; CLIScast01_btl.O2_umolL; ...
    ];

O2_umolkg = O2_umolL .* (1000 ./ Dens);


Chl = [WLIScast01_btl.Fluor; MID4cast01_btl.Fluor; EXRXcast01_btl.Fluor; MID4cast02_btl.Fluor; ...
    MID4cast03_btl.Fluor; MID4cast04_btl.Fluor; MID4cast05_btl.Fluor; MID4cast06_btl.Fluor; ...
    MID4cast07_btl.Fluor; MID4cast08_btl.Fluor; MID4cast09_btl.Fluor; MID4cast10_btl.Fluor; ...
    MID4cast11_btl.Fluor; MID4cast12_btl.Fluor; MID4cast13_btl.Fluor; EXR1cast01_btl.Fluor; ...
    EXRXcast02_btl.Fluor; MID3cast01_btl.Fluor; MID4cast14_btl.Fluor; MID5cast01_btl.Fluor; ...
    WLIScast02_btl.Fluor; WLI6cast01_btl.Fluor; WLI7cast01_btl.Fluor; ARTGcast01_btl.Fluor; CLIScast01_btl.Fluor; ...
    ];

pH = [WLIScast01_btl.pH; MID4cast01_btl.pH; EXRXcast01_btl.pH; MID4cast02_btl.pH; ...
    MID4cast03_btl.pH; MID4cast04_btl.pH; MID4cast05_btl.pH; MID4cast06_btl.pH; ...
    MID4cast07_btl.pH; MID4cast08_btl.pH; MID4cast09_btl.pH; MID4cast10_btl.pH; ...
    MID4cast11_btl.pH; MID4cast12_btl.pH; MID4cast13_btl.pH; EXR1cast01_btl.pH; ...
    EXRXcast02_btl.pH; MID3cast01_btl.pH; MID4cast14_btl.pH; MID5cast01_btl.pH; ...
    WLIScast02_btl.pH; WLI6cast01_btl.pH; WLI7cast01_btl.pH; ARTGcast01_btl.pH; CLIScast01_btl.pH; ...
    ];



%%
btl = table(Cast);
btl.Niskin = Niskin;
btl.Depth = Depth;
btl.P = P;
btl.S = S;
btl.T = T;
btl.PTemp = PTemp;
btl.Dens = Dens;
btl.PDen = PDen;
btl.O2_umolkg = O2_umolkg;
%btl.O2_umolkg_A = O2_umolkg_A;
btl.Chl = Chl;
btl.pH = pH;


%%

load LISMay2024CastData.mat;
LIS = LISMay2024CastData;
infmt = 'MMM dd yyyy HH:mm:ss';
LIS.datetime=datetime(LIS.DateTime_UTC,'inputFormat',infmt);


%btl.Station = nan.*btl.Cast;

for i = 1:length(btl.Cast)
    A=find(LIS.Cast==btl.Cast(i));
    btl.Station(i) = LIS.Station(A);
    btl.Lat(i) = LIS.Lat(A);   
    btl.Lon(i) = LIS.Lon(A);   
    btl.datetime(i) = LIS.datetime(A);   
    btl.CH4N2Ocast(i) = LIS.CH4N2O(A);    
end;

btl = movevars(btl, "Station", "Before", "Cast");
%%
LISMay24btl = btl;

save LISMay24btl.mat LISMay24btl;
