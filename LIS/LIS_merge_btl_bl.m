
%%

%%

load WLIScast01_btl.mat;
load WLIScast01_fireseq.mat;
WLIScast01_btl.Niskin = nan.*WLIScast01_btl.FireSeq;
WLIScast01_btl.Cast = nan.*WLIScast01_btl.FireSeq;

for i = 1:numel(WLIScast01_fireseq.FireSeq)
    ifs = find(WLIScast01_fireseq.FireSeq==WLIScast01_btl.FireSeq(i));
    WLIScast01_btl.Niskin(i) = WLIScast01_fireseq.Niskin(ifs);
    WLIScast01_btl.Cast(i) = 1;
end;

load MID4cast01_btl.mat;
load MID4cast01_fireseq.mat;
MID4cast01_btl.Niskin = nan.*MID4cast01_btl.FireSeq;
MID4cast01_btl.Cast = nan.*MID4cast01_btl.FireSeq;

for i = 1:numel(MID4cast01_fireseq.FireSeq)
    ifs = find(MID4cast01_fireseq.FireSeq==MID4cast01_btl.FireSeq(i));
    MID4cast01_btl.Niskin(i) = MID4cast01_fireseq.Niskin(ifs);
    MID4cast01_btl.Cast(i) = 2;
end;


load EXCRcast01_btl.mat;
load EXCRcast01_fireseq.mat;
EXCRcast01_btl.Niskin = nan.*EXCRcast01_btl.FireSeq;
EXCRcast01_btl.Cast = nan.*EXCRcast01_btl.FireSeq;

for i = 1:numel(EXCRcast01_fireseq.FireSeq)
    ifs = find(EXCRcast01_fireseq.FireSeq==EXCRcast01_btl.FireSeq(i));
    EXCRcast01_btl.Niskin(i) = EXCRcast01_fireseq.Niskin(ifs);
    EXCRcast01_btl.Cast(i) = 3;
end;



load MID4cast02_btl.mat;
load MID4cast02_fireseq.mat;
MID4cast02_btl.Niskin = nan.*MID4cast02_btl.FireSeq;
MID4cast02_btl.Cast = nan.*MID4cast02_btl.FireSeq;

for i = 1:numel(MID4cast02_fireseq.FireSeq)
    ifs = find(MID4cast02_fireseq.FireSeq==MID4cast02_btl.FireSeq(i));
    MID4cast02_btl.Niskin(i) = MID4cast02_fireseq.Niskin(ifs);
    MID4cast02_btl.Cast(i) = 4;
end;

load MID4cast03_btl.mat;
load MID4cast03_fireseq.mat;
MID4cast03_btl.Niskin = nan.*MID4cast03_btl.FireSeq;
MID4cast03_btl.Cast = nan.*MID4cast03_btl.FireSeq;

for i = 1:numel(MID4cast03_fireseq.FireSeq)
    ifs = find(MID4cast03_fireseq.FireSeq==MID4cast03_btl.FireSeq(i));
    MID4cast03_btl.Niskin(i) = MID4cast03_fireseq.Niskin(ifs);
    MID4cast03_btl.Cast(i) = 5;
end;

load MID4cast04_btl.mat;
load MID4cast04_fireseq.mat;
MID4cast04_btl.Niskin = nan.*MID4cast04_btl.FireSeq;
MID4cast04_btl.Cast = nan.*MID4cast04_btl.FireSeq;

for i = 1:numel(MID4cast04_fireseq.FireSeq)
    ifs = find(MID4cast04_fireseq.FireSeq==MID4cast04_btl.FireSeq(i));
    MID4cast04_btl.Niskin(i) = MID4cast04_fireseq.Niskin(ifs);
    MID4cast04_btl.Cast(i) = 6;
end;

load MID4cast05_btl.mat;
load MID4cast05_fireseq.mat;
MID4cast05_btl.Niskin = nan.*MID4cast05_btl.FireSeq;
MID4cast05_btl.Cast = nan.*MID4cast05_btl.FireSeq;

for i = 1:numel(MID4cast05_fireseq.FireSeq)
    ifs = find(MID4cast05_fireseq.FireSeq==MID4cast05_btl.FireSeq(i));
    MID4cast05_btl.Niskin(i) = MID4cast05_fireseq.Niskin(ifs);
    MID4cast05_btl.Cast(i) = 7;
end;

load MID4cast06_btl.mat;
load MID4cast06_fireseq.mat;
MID4cast06_btl.Niskin = nan.*MID4cast06_btl.FireSeq;
MID4cast06_btl.Cast = nan.*MID4cast06_btl.FireSeq;

for i = 1:numel(MID4cast06_fireseq.FireSeq)
    ifs = find(MID4cast06_fireseq.FireSeq==MID4cast06_btl.FireSeq(i));
    MID4cast06_btl.Niskin(i) = MID4cast06_fireseq.Niskin(ifs);
    MID4cast06_btl.Cast(i) = 8;
end;

load MID4cast07_btl.mat;
load MID4cast07_fireseq.mat;
MID4cast07_btl.Niskin = nan.*MID4cast07_btl.FireSeq;
MID4cast07_btl.Cast = nan.*MID4cast07_btl.FireSeq;

for i = 1:numel(MID4cast07_fireseq.FireSeq)
    ifs = find(MID4cast07_fireseq.FireSeq==MID4cast07_btl.FireSeq(i));
    MID4cast07_btl.Niskin(i) = MID4cast07_fireseq.Niskin(ifs);
    MID4cast07_btl.Cast(i) = 9;
end;

load MID4cast08_btl.mat;
load MID4cast08_fireseq.mat;
MID4cast08_btl.Niskin = nan.*MID4cast08_btl.FireSeq;
MID4cast08_btl.Cast = nan.*MID4cast08_btl.FireSeq;

for i = 1:numel(MID4cast08_fireseq.FireSeq)
    ifs = find(MID4cast08_fireseq.FireSeq==MID4cast08_btl.FireSeq(i));
    MID4cast08_btl.Niskin(i) = MID4cast08_fireseq.Niskin(ifs);
    MID4cast08_btl.Cast(i) = 10;
end;

load MID4cast09_btl.mat;
load MID4cast09_fireseq.mat;
MID4cast09_btl.Niskin = nan.*MID4cast09_btl.FireSeq;
MID4cast09_btl.Cast = nan.*MID4cast09_btl.FireSeq;

for i = 1:numel(MID4cast09_fireseq.FireSeq)
    ifs = find(MID4cast09_fireseq.FireSeq==MID4cast09_btl.FireSeq(i));
    MID4cast09_btl.Niskin(i) = MID4cast09_fireseq.Niskin(ifs);
    MID4cast09_btl.Cast(i) = 11;
end;

load MID4cast10_btl.mat;
load MID4cast10_fireseq.mat;
MID4cast10_btl.Niskin = nan.*MID4cast10_btl.FireSeq;
MID4cast10_btl.Cast = nan.*MID4cast10_btl.FireSeq;

for i = 1:numel(MID4cast10_fireseq.FireSeq)
    ifs = find(MID4cast10_fireseq.FireSeq==MID4cast10_btl.FireSeq(i));
    MID4cast10_btl.Niskin(i) = MID4cast10_fireseq.Niskin(ifs);
    MID4cast10_btl.Cast(i) = 12;
end;


load MID4cast11_btl.mat;
load MID4cast11_fireseq.mat;
MID4cast11_btl.Niskin = nan.*MID4cast11_btl.FireSeq;
MID4cast11_btl.Cast = nan.*MID4cast11_btl.FireSeq;

for i = 1:numel(MID4cast11_fireseq.FireSeq)
    ifs = find(MID4cast11_fireseq.FireSeq==MID4cast11_btl.FireSeq(i));
    MID4cast11_btl.Niskin(i) = MID4cast11_fireseq.Niskin(ifs);
    MID4cast11_btl.Cast(i) = 13;
end;

load MID4cast12_btl.mat;
load MID4cast12_fireseq.mat;
MID4cast12_btl.Niskin = nan.*MID4cast12_btl.FireSeq;
MID4cast12_btl.Cast = nan.*MID4cast12_btl.FireSeq;

for i = 1:numel(MID4cast12_fireseq.FireSeq)
    ifs = find(MID4cast12_fireseq.FireSeq==MID4cast12_btl.FireSeq(i));
    MID4cast12_btl.Niskin(i) = MID4cast12_fireseq.Niskin(ifs);
    MID4cast12_btl.Cast(i) = 14;
end;

load MID4cast13_btl.mat;
load MID4cast13_fireseq.mat;
MID4cast13_btl.Niskin = nan.*MID4cast13_btl.FireSeq;
MID4cast13_btl.Cast = nan.*MID4cast13_btl.FireSeq;

for i = 1:numel(MID4cast13_fireseq.FireSeq)
    ifs = find(MID4cast13_fireseq.FireSeq==MID4cast13_btl.FireSeq(i));
    MID4cast13_btl.Niskin(i) = MID4cast13_fireseq.Niskin(ifs);
    MID4cast13_btl.Cast(i) = 15;
end;


load EXR1cast01_btl.mat;
load EXR1cast01_fireseq.mat;
EXR1cast01_btl.Niskin = nan.*EXR1cast01_btl.FireSeq;
EXR1cast01_btl.Cast = nan.*EXR1cast01_btl.FireSeq;

for i = 1:numel(EXR1cast01_fireseq.FireSeq)
    ifs = find(EXR1cast01_fireseq.FireSeq==EXR1cast01_btl.FireSeq(i));
    EXR1cast01_btl.Niskin(i) = EXR1cast01_fireseq.Niskin(ifs);
    EXR1cast01_btl.Cast(i) = 16;
end;


load EXCRcast02_btl.mat;
load EXCRcast02_fireseq.mat;
EXCRcast02_btl.Niskin = nan.*EXCRcast02_btl.FireSeq;
EXCRcast02_btl.Cast = nan.*EXCRcast02_btl.FireSeq;

for i = 1:numel(EXCRcast02_fireseq.FireSeq)
    ifs = find(EXCRcast02_fireseq.FireSeq==EXCRcast02_btl.FireSeq(i));
    EXCRcast02_btl.Niskin(i) = EXCRcast02_fireseq.Niskin(ifs);
    EXCRcast02_btl.Cast(i) = 17;
end;


load MID3cast01_btl.mat;
load MID3cast01_fireseq.mat;
MID3cast01_btl.Niskin = nan.*MID3cast01_btl.FireSeq;
MID3cast01_btl.Cast = nan.*MID3cast01_btl.FireSeq;

for i = 1:numel(MID3cast01_fireseq.FireSeq)
    ifs = find(MID3cast01_fireseq.FireSeq==MID3cast01_btl.FireSeq(i));
    MID3cast01_btl.Niskin(i) = MID3cast01_fireseq.Niskin(ifs);
    MID3cast01_btl.Cast(i) = 18;
end;


load MID4cast14_btl.mat;
load MID4cast14_fireseq.mat;
MID4cast14_btl.Niskin = nan.*MID4cast14_btl.FireSeq;
MID4cast14_btl.Cast = nan.*MID4cast14_btl.FireSeq;

for i = 1:numel(MID4cast14_fireseq.FireSeq)
    ifs = find(MID4cast14_fireseq.FireSeq==MID4cast14_btl.FireSeq(i));
    MID4cast14_btl.Niskin(i) = MID4cast14_fireseq.Niskin(ifs);
    MID4cast14_btl.Cast(i) = 19;
end;


load MID5cast01_btl.mat;
load MID5cast01_fireseq.mat;
MID5cast01_btl.Niskin = nan.*MID5cast01_btl.FireSeq;
MID5cast01_btl.Cast = nan.*MID5cast01_btl.FireSeq;

for i = 1:numel(MID5cast01_fireseq.FireSeq)
    ifs = find(MID5cast01_fireseq.FireSeq==MID5cast01_btl.FireSeq(i));
    MID5cast01_btl.Niskin(i) = MID5cast01_fireseq.Niskin(ifs);
    MID5cast01_btl.Cast(i) = 20;
end;

load WLIScast02_btl.mat;
load WLIScast02_fireseq.mat;
WLIScast02_btl.Niskin = nan.*WLIScast02_btl.FireSeq;
WLIScast02_btl.Cast = nan.*WLIScast02_btl.FireSeq;

for i = 1:numel(WLIScast02_fireseq.FireSeq)
    ifs = find(WLIScast02_fireseq.FireSeq==WLIScast02_btl.FireSeq(i));
    WLIScast02_btl.Niskin(i) = WLIScast02_fireseq.Niskin(ifs);
    WLIScast02_btl.Cast(i) = 21;
end;


load WLI6cast01_btl.mat;
load WLI6cast01_fireseq.mat;
WLI6cast01_btl.Niskin = nan.*WLI6cast01_btl.FireSeq;
WLI6cast01_btl.Cast = nan.*WLI6cast01_btl.FireSeq;

for i = 1:numel(WLI6cast01_fireseq.FireSeq)
    ifs = find(WLI6cast01_fireseq.FireSeq==WLI6cast01_btl.FireSeq(i));
    WLI6cast01_btl.Niskin(i) = WLI6cast01_fireseq.Niskin(ifs);
    WLI6cast01_btl.Cast(i) = 22;
end;


load WLI7cast01_btl.mat;
load WLI7cast01_fireseq.mat;
WLI7cast01_btl.Niskin = nan.*WLI7cast01_btl.FireSeq;
WLI7cast01_btl.Cast = nan.*WLI7cast01_btl.FireSeq;

for i = 1:numel(WLI7cast01_fireseq.FireSeq)
    ifs = find(WLI7cast01_fireseq.FireSeq==WLI7cast01_btl.FireSeq(i));
    WLI7cast01_btl.Niskin(i) = WLI7cast01_fireseq.Niskin(ifs);
    WLI7cast01_btl.Cast(i) = 23;
end;


%%


Cast = [WLIScast01_btl.Cast; MID4cast01_btl.Cast; EXCRcast01_btl.Cast; MID4cast02_btl.Cast; ...
    MID4cast03_btl.Cast; MID4cast04_btl.Cast; MID4cast05_btl.Cast; MID4cast06_btl.Cast; ...
    MID4cast07_btl.Cast; MID4cast08_btl.Cast; MID4cast09_btl.Cast; MID4cast10_btl.Cast; ...
    MID4cast11_btl.Cast; MID4cast12_btl.Cast; MID4cast13_btl.Cast; EXR1cast01_btl.Cast; ...
    EXCRcast02_btl.Cast; MID3cast01_btl.Cast; MID4cast14_btl.Cast; MID5cast01_btl.Cast; ...
    WLIScast02_btl.Cast; WLI6cast01_btl.Cast; WLI7cast01_btl.Cast; ...
    ];

Niskin = [WLIScast01_btl.Niskin; MID4cast01_btl.Niskin; EXCRcast01_btl.Niskin; MID4cast02_btl.Niskin; ...
    MID4cast03_btl.Niskin; MID4cast04_btl.Niskin; MID4cast05_btl.Niskin; MID4cast06_btl.Niskin; ...
    MID4cast07_btl.Niskin; MID4cast08_btl.Niskin; MID4cast09_btl.Niskin; MID4cast10_btl.Niskin; ...
    MID4cast11_btl.Niskin; MID4cast12_btl.Niskin; MID4cast13_btl.Niskin; EXR1cast01_btl.Niskin; ...
    EXCRcast02_btl.Niskin; MID3cast01_btl.Niskin; MID4cast14_btl.Niskin; MID5cast01_btl.Niskin; ...
    WLIScast02_btl.Niskin; WLI6cast01_btl.Niskin; WLI7cast01_btl.Niskin; ...
    ];

Depth = [WLIScast01_btl.Depth; MID4cast01_btl.Depth; EXCRcast01_btl.Depth; MID4cast02_btl.Depth; ...
    MID4cast03_btl.Depth; MID4cast04_btl.Depth; MID4cast05_btl.Depth; MID4cast06_btl.Depth; ...
    MID4cast07_btl.Depth; MID4cast08_btl.Depth; MID4cast09_btl.Depth; MID4cast10_btl.Depth; ...
    MID4cast11_btl.Depth; MID4cast12_btl.Depth; MID4cast13_btl.Depth; EXR1cast01_btl.Depth; ...
    EXCRcast02_btl.Depth; MID3cast01_btl.Depth; MID4cast14_btl.Depth; MID5cast01_btl.Depth; ...
    WLIScast02_btl.Depth; WLI6cast01_btl.Depth; WLI7cast01_btl.Depth; ...
    ];


P = [WLIScast01_btl.P; MID4cast01_btl.P; EXCRcast01_btl.P; MID4cast02_btl.P; ...
    MID4cast03_btl.P; MID4cast04_btl.P; MID4cast05_btl.P; MID4cast06_btl.P; ...
    MID4cast07_btl.P; MID4cast08_btl.P; MID4cast09_btl.P; MID4cast10_btl.P; ...
    MID4cast11_btl.P; MID4cast12_btl.P; MID4cast13_btl.P; EXR1cast01_btl.P; ...
    EXCRcast02_btl.P; MID3cast01_btl.P; MID4cast14_btl.P; MID5cast01_btl.P; ...
    WLIScast02_btl.P; WLI6cast01_btl.P; WLI7cast01_btl.P; ...
    ];

S = [WLIScast01_btl.S; MID4cast01_btl.S; EXCRcast01_btl.S; MID4cast02_btl.S; ...
    MID4cast03_btl.S; MID4cast04_btl.S; MID4cast05_btl.S; MID4cast06_btl.S; ...
    MID4cast07_btl.S; MID4cast08_btl.S; MID4cast09_btl.S; MID4cast10_btl.S; ...
    MID4cast11_btl.S; MID4cast12_btl.S; MID4cast13_btl.S; EXR1cast01_btl.S; ...
    EXCRcast02_btl.S; MID3cast01_btl.S; MID4cast14_btl.S; MID5cast01_btl.S; ...
    WLIScast02_btl.S; WLI6cast01_btl.S; WLI7cast01_btl.S; ...
    ];

T = [WLIScast01_btl.T; MID4cast01_btl.T; EXCRcast01_btl.T; MID4cast02_btl.T; ...
    MID4cast03_btl.T; MID4cast04_btl.T; MID4cast05_btl.T; MID4cast06_btl.T; ...
    MID4cast07_btl.T; MID4cast08_btl.T; MID4cast09_btl.T; MID4cast10_btl.T; ...
    MID4cast11_btl.T; MID4cast12_btl.T; MID4cast13_btl.T; EXR1cast01_btl.T; ...
    EXCRcast02_btl.T; MID3cast01_btl.T; MID4cast14_btl.T; MID5cast01_btl.T; ...
    WLIScast02_btl.T; WLI6cast01_btl.T; WLI7cast01_btl.T; ...
    ];

PTemp = [WLIScast01_btl.PTemp; MID4cast01_btl.PTemp; EXCRcast01_btl.PTemp; MID4cast02_btl.PTemp; ...
    MID4cast03_btl.PTemp; MID4cast04_btl.PTemp; MID4cast05_btl.PTemp; MID4cast06_btl.PTemp; ...
    MID4cast07_btl.PTemp; MID4cast08_btl.PTemp; MID4cast09_btl.PTemp; MID4cast10_btl.PTemp; ...
    MID4cast11_btl.PTemp; MID4cast12_btl.PTemp; MID4cast13_btl.PTemp; EXR1cast01_btl.PTemp; ...
    EXCRcast02_btl.PTemp; MID3cast01_btl.PTemp; MID4cast14_btl.PTemp; MID5cast01_btl.PTemp; ...
    WLIScast02_btl.PTemp; WLI6cast01_btl.PTemp; WLI7cast01_btl.PTemp; ...
    ];

Dens = [WLIScast01_btl.Dens; MID4cast01_btl.Dens; EXCRcast01_btl.Dens; MID4cast02_btl.Dens; ...
    MID4cast03_btl.Dens; MID4cast04_btl.Dens; MID4cast05_btl.Dens; MID4cast06_btl.Dens; ...
    MID4cast07_btl.Dens; MID4cast08_btl.Dens; MID4cast09_btl.Dens; MID4cast10_btl.Dens; ...
    MID4cast11_btl.Dens; MID4cast12_btl.Dens; MID4cast13_btl.Dens; EXR1cast01_btl.Dens; ...
    EXCRcast02_btl.Dens; MID3cast01_btl.Dens; MID4cast14_btl.Dens; MID5cast01_btl.Dens; ...
    WLIScast02_btl.Dens; WLI6cast01_btl.Dens; WLI7cast01_btl.Dens; ...
    ];

PDen = [WLIScast01_btl.PDen; MID4cast01_btl.PDen; EXCRcast01_btl.PDen; MID4cast02_btl.PDen; ...
    MID4cast03_btl.PDen; MID4cast04_btl.PDen; MID4cast05_btl.PDen; MID4cast06_btl.PDen; ...
    MID4cast07_btl.PDen; MID4cast08_btl.PDen; MID4cast09_btl.PDen; MID4cast10_btl.PDen; ...
    MID4cast11_btl.PDen; MID4cast12_btl.PDen; MID4cast13_btl.PDen; EXR1cast01_btl.PDen; ...
    EXCRcast02_btl.PDen; MID3cast01_btl.PDen; MID4cast14_btl.PDen; MID5cast01_btl.PDen; ...
    WLIScast02_btl.PDen; WLI6cast01_btl.PDen; WLI7cast01_btl.PDen; ...
    ];

O2_umolkg = [WLIScast01_btl.O2_umolkg; MID4cast01_btl.O2_umolkg; EXCRcast01_btl.O2_umolkg; MID4cast02_btl.O2_umolkg; ...
    MID4cast03_btl.O2_umolkg; MID4cast04_btl.O2_umolkg; MID4cast05_btl.O2_umolkg; MID4cast06_btl.O2_umolkg; ...
    MID4cast07_btl.O2_umolkg; MID4cast08_btl.O2_umolkg; MID4cast09_btl.O2_umolkg; MID4cast10_btl.O2_umolkg; ...
    MID4cast11_btl.O2_umolkg; MID4cast12_btl.O2_umolkg; MID4cast13_btl.O2_umolkg; EXR1cast01_btl.O2_umolkg; ...
    EXCRcast02_btl.O2_umolkg; MID3cast01_btl.O2_umolkg; MID4cast14_btl.O2_umolkg; MID5cast01_btl.O2_umolkg; ...
    WLIScast02_btl.O2_umolkg; WLI6cast01_btl.O2_umolkg; WLI7cast01_btl.O2_umolkg; ...
    ];

O2_umolkg_A = [WLIScast01_btl.O2_umolkg_A; MID4cast01_btl.O2_umolkg_A; EXCRcast01_btl.O2_umolkg_A; MID4cast02_btl.O2_umolkg_A; ...
    MID4cast03_btl.O2_umolkg_A; MID4cast04_btl.O2_umolkg_A; MID4cast05_btl.O2_umolkg_A; MID4cast06_btl.O2_umolkg_A; ...
    MID4cast07_btl.O2_umolkg_A; MID4cast08_btl.O2_umolkg_A; MID4cast09_btl.O2_umolkg_A; MID4cast10_btl.O2_umolkg_A; ...
    MID4cast11_btl.O2_umolkg_A; MID4cast12_btl.O2_umolkg_A; MID4cast13_btl.O2_umolkg_A; EXR1cast01_btl.O2_umolkg_A; ...
    EXCRcast02_btl.O2_umolkg_A; MID3cast01_btl.O2_umolkg_A; MID4cast14_btl.O2_umolkg_A; MID5cast01_btl.O2_umolkg_A; ...
    WLIScast02_btl.O2_umolkg_A; WLI6cast01_btl.O2_umolkg_A; WLI7cast01_btl.O2_umolkg_A; ...
    ];

Chl = [WLIScast01_btl.Chl; MID4cast01_btl.Chl; EXCRcast01_btl.Chl; MID4cast02_btl.Chl; ...
    MID4cast03_btl.Chl; MID4cast04_btl.Chl; MID4cast05_btl.Chl; MID4cast06_btl.Chl; ...
    MID4cast07_btl.Chl; MID4cast08_btl.Chl; MID4cast09_btl.Chl; MID4cast10_btl.Chl; ...
    MID4cast11_btl.Chl; MID4cast12_btl.Chl; MID4cast13_btl.Chl; EXR1cast01_btl.Chl; ...
    EXCRcast02_btl.Chl; MID3cast01_btl.Chl; MID4cast14_btl.Chl; MID5cast01_btl.Chl; ...
    WLIScast02_btl.Chl; WLI6cast01_btl.Chl; WLI7cast01_btl.Chl; ...
    ];

pH = [WLIScast01_btl.pH; MID4cast01_btl.pH; EXCRcast01_btl.pH; MID4cast02_btl.pH; ...
    MID4cast03_btl.pH; MID4cast04_btl.pH; MID4cast05_btl.pH; MID4cast06_btl.pH; ...
    MID4cast07_btl.pH; MID4cast08_btl.pH; MID4cast09_btl.pH; MID4cast10_btl.pH; ...
    MID4cast11_btl.pH; MID4cast12_btl.pH; MID4cast13_btl.pH; EXR1cast01_btl.pH; ...
    EXCRcast02_btl.pH; MID3cast01_btl.pH; MID4cast14_btl.pH; MID5cast01_btl.pH; ...
    WLIScast02_btl.pH; WLI6cast01_btl.pH; WLI7cast01_btl.pH; ...
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
btl.O2_umolkg_A = O2_umolkg_A;
btl.Chl = Chl;
btl.pH = pH;


%%

load LISAug2023CastData.mat;

infmt = 'MMM dd yyyy HH:mm:ss';
LISAug2023CastData.datetime=datetime(LISAug2023CastData.DateTime_UTC,'inputFormat',infmt);


%btl.Station = nan.*btl.Cast;

for i = 1:length(btl.Cast)
    A=find(LISAug2023CastData.Cast==btl.Cast(i));
    btl.Station(i) = LISAug2023CastData.Station(A);
    btl.Lat(i) = LISAug2023CastData.Lat(A);   
    btl.Lon(i) = LISAug2023CastData.Lon(A);   
    btl.datetime(i) = LISAug2023CastData.datetime(A);   
    btl.CH4N2Ocast(i) = LISAug2023CastData.CH4N2O(A);    
end;

btl = movevars(btl, "Station", "Before", "Cast");
%%
LISAug23btl = btl;

save LISAug23btl.mat LISAug23btl;

%%

A = join(EXCRcast01_btl,EXCRcast02_btl);

% next steps
% make table with cast number, lat, lon and cast time in UTC


%%

Month = 'EXCR-cast01,EXCR-cast02,03,04,05,06,07,08,09,10,11,12';

numFiles = 12;
for n = 1:numFiles
   randomData = rand(n);
   currentFile = sprintf('myfile%d.mat',n);
   load(currentFile)
end



%%

monthsArray = strsplit(Month,',');
for i = 1:numel(monthsArray)
Varnames{i} = matlab.lang.makeValidName(strcat('Indiv_Reg_',monthsArray{i}));
myStruct.(Varnames{i}) = randi(20,1,1);
end
myStruct.(Varnames{1,1}) % should give you a value of a random number
myStruct.Indiv_Reg_01 % same result above


%%
A = ['EXCR-cast01'
    'EXCR-cast02']

