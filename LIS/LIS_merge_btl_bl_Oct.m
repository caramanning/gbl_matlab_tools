
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


load EXRXcast01_btl.mat;
load EXRXcast01_fireseq.mat;
EXRXcast01_btl.Niskin = nan.*EXRXcast01_btl.FireSeq;
EXRXcast01_btl.Cast = nan.*EXRXcast01_btl.FireSeq;

for i = 1:numel(EXRXcast01_fireseq.FireSeq)
    ifs = find(EXRXcast01_fireseq.FireSeq==EXRXcast01_btl.FireSeq(i));
    EXRXcast01_btl.Niskin(i) = EXRXcast01_fireseq.Niskin(ifs);
    EXRXcast01_btl.Cast(i) = 3;
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


load EXCR1cast01_btl.mat;
load EXCR1cast01_fireseq.mat;
EXCR1cast01_btl.Niskin = nan.*EXCR1cast01_btl.FireSeq;
EXCR1cast01_btl.Cast = nan.*EXCR1cast01_btl.FireSeq;

for i = 1:numel(EXCR1cast01_fireseq.FireSeq)
    ifs = find(EXCR1cast01_fireseq.FireSeq==EXCR1cast01_btl.FireSeq(i));
    EXCR1cast01_btl.Niskin(i) = EXCR1cast01_fireseq.Niskin(ifs);
    EXCR1cast01_btl.Cast(i) = 16;
end;


load EXRXcast02_btl.mat;
load EXRXcast02_fireseq.mat;
EXRXcast02_btl.Niskin = nan.*EXRXcast02_btl.FireSeq;
EXRXcast02_btl.Cast = nan.*EXRXcast02_btl.FireSeq;

for i = 1:numel(EXRXcast02_fireseq.FireSeq)
    ifs = find(EXRXcast02_fireseq.FireSeq==EXRXcast02_btl.FireSeq(i));
    EXRXcast02_btl.Niskin(i) = EXRXcast02_fireseq.Niskin(ifs);
    EXRXcast02_btl.Cast(i) = 17;
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

%
% not sampled in October
% load WLI7cast01_btl.mat;
% load WLI7cast01_fireseq.mat;
% WLI7cast01_btl.Niskin = nan.*WLI7cast01_btl.FireSeq;
% WLI7cast01_btl.Cast = nan.*WLI7cast01_btl.FireSeq;
% 
% for i = 1:numel(WLI7cast01_fireseq.FireSeq)
%     ifs = find(WLI7cast01_fireseq.FireSeq==WLI7cast01_btl.FireSeq(i));
%     WLI7cast01_btl.Niskin(i) = WLI7cast01_fireseq.Niskin(ifs);
%     WLI7cast01_btl.Cast(i) = 23;
% end;


%%


Cast = [WLIScast01_btl.Cast; MID4cast01_btl.Cast; EXRXcast01_btl.Cast; MID4cast02_btl.Cast; ...
    MID4cast03_btl.Cast; MID4cast04_btl.Cast; MID4cast05_btl.Cast; MID4cast06_btl.Cast; ...
    MID4cast07_btl.Cast; MID4cast08_btl.Cast; MID4cast09_btl.Cast; MID4cast10_btl.Cast; ...
    MID4cast11_btl.Cast; MID4cast12_btl.Cast; MID4cast13_btl.Cast; EXCR1cast01_btl.Cast; ...
    EXRXcast02_btl.Cast; MID3cast01_btl.Cast; MID4cast14_btl.Cast; MID5cast01_btl.Cast; ...
    WLIScast02_btl.Cast; WLI6cast01_btl.Cast; ...
    ];

Niskin = [WLIScast01_btl.Niskin; MID4cast01_btl.Niskin; EXRXcast01_btl.Niskin; MID4cast02_btl.Niskin; ...
    MID4cast03_btl.Niskin; MID4cast04_btl.Niskin; MID4cast05_btl.Niskin; MID4cast06_btl.Niskin; ...
    MID4cast07_btl.Niskin; MID4cast08_btl.Niskin; MID4cast09_btl.Niskin; MID4cast10_btl.Niskin; ...
    MID4cast11_btl.Niskin; MID4cast12_btl.Niskin; MID4cast13_btl.Niskin; EXCR1cast01_btl.Niskin; ...
    EXRXcast02_btl.Niskin; MID3cast01_btl.Niskin; MID4cast14_btl.Niskin; MID5cast01_btl.Niskin; ...
    WLIScast02_btl.Niskin; WLI6cast01_btl.Niskin; ...
    ];

Depth = [WLIScast01_btl.DepSM; MID4cast01_btl.DepSM; EXRXcast01_btl.DepSM; MID4cast02_btl.DepSM; ...
    MID4cast03_btl.DepSM; MID4cast04_btl.DepSM; MID4cast05_btl.DepSM; MID4cast06_btl.DepSM; ...
    MID4cast07_btl.DepSM; MID4cast08_btl.DepSM; MID4cast09_btl.DepSM; MID4cast10_btl.DepSM; ...
    MID4cast11_btl.DepSM; MID4cast12_btl.DepSM; MID4cast13_btl.DepSM; EXCR1cast01_btl.DepSM; ...
    EXRXcast02_btl.DepSM; MID3cast01_btl.DepSM; MID4cast14_btl.DepSM; MID5cast01_btl.DepSM; ...
    WLIScast02_btl.DepSM; WLI6cast01_btl.DepSM; ...
    ];


P = [WLIScast01_btl.PrDM; MID4cast01_btl.PrDM; EXRXcast01_btl.PrDM; MID4cast02_btl.PrDM; ...
    MID4cast03_btl.PrDM; MID4cast04_btl.PrDM; MID4cast05_btl.PrDM; MID4cast06_btl.PrDM; ...
    MID4cast07_btl.PrDM; MID4cast08_btl.PrDM; MID4cast09_btl.PrDM; MID4cast10_btl.PrDM; ...
    MID4cast11_btl.PrDM; MID4cast12_btl.PrDM; MID4cast13_btl.PrDM; EXCR1cast01_btl.PrDM; ...
    EXRXcast02_btl.PrDM; MID3cast01_btl.PrDM; MID4cast14_btl.PrDM; MID5cast01_btl.PrDM; ...
    WLIScast02_btl.PrDM; WLI6cast01_btl.PrDM; ...
    ];

S = [WLIScast01_btl.Sal00; MID4cast01_btl.Sal00; EXRXcast01_btl.Sal00; MID4cast02_btl.Sal00; ...
    MID4cast03_btl.Sal00; MID4cast04_btl.Sal00; MID4cast05_btl.Sal00; MID4cast06_btl.Sal00; ...
    MID4cast07_btl.Sal00; MID4cast08_btl.Sal00; MID4cast09_btl.Sal00; MID4cast10_btl.Sal00; ...
    MID4cast11_btl.Sal00; MID4cast12_btl.Sal00; MID4cast13_btl.Sal00; EXCR1cast01_btl.Sal00; ...
    EXRXcast02_btl.Sal00; MID3cast01_btl.Sal00; MID4cast14_btl.Sal00; MID5cast01_btl.Sal00; ...
    WLIScast02_btl.Sal00; WLI6cast01_btl.Sal00; ...
    ];

T = [WLIScast01_btl.T090C; MID4cast01_btl.T090C; EXRXcast01_btl.T090C; MID4cast02_btl.T090C; ...
    MID4cast03_btl.T090C; MID4cast04_btl.T090C; MID4cast05_btl.T090C; MID4cast06_btl.T090C; ...
    MID4cast07_btl.T090C; MID4cast08_btl.T090C; MID4cast09_btl.T090C; MID4cast10_btl.T090C; ...
    MID4cast11_btl.T090C; MID4cast12_btl.T090C; MID4cast13_btl.T090C; EXCR1cast01_btl.T090C; ...
    EXRXcast02_btl.T090C; MID3cast01_btl.T090C; MID4cast14_btl.T090C; MID5cast01_btl.T090C; ...
    WLIScast02_btl.T090C; WLI6cast01_btl.T090C; ...
    ];

PTemp = [WLIScast01_btl.Potemp090C; MID4cast01_btl.Potemp090C; EXRXcast01_btl.Potemp090C; MID4cast02_btl.Potemp090C; ...
    MID4cast03_btl.Potemp090C; MID4cast04_btl.Potemp090C; MID4cast05_btl.Potemp090C; MID4cast06_btl.Potemp090C; ...
    MID4cast07_btl.Potemp090C; MID4cast08_btl.Potemp090C; MID4cast09_btl.Potemp090C; MID4cast10_btl.Potemp090C; ...
    MID4cast11_btl.Potemp090C; MID4cast12_btl.Potemp090C; MID4cast13_btl.Potemp090C; EXCR1cast01_btl.Potemp090C; ...
    EXRXcast02_btl.Potemp090C; MID3cast01_btl.Potemp090C; MID4cast14_btl.Potemp090C; MID5cast01_btl.Potemp090C; ...
    WLIScast02_btl.Potemp090C; WLI6cast01_btl.Potemp090C; ...
    ];

Dens = [WLIScast01_btl.Density00; MID4cast01_btl.Density00; EXRXcast01_btl.Density00; MID4cast02_btl.Density00; ...
    MID4cast03_btl.Density00; MID4cast04_btl.Density00; MID4cast05_btl.Density00; MID4cast06_btl.Density00; ...
    MID4cast07_btl.Density00; MID4cast08_btl.Density00; MID4cast09_btl.Density00; MID4cast10_btl.Density00; ...
    MID4cast11_btl.Density00; MID4cast12_btl.Density00; MID4cast13_btl.Density00; EXCR1cast01_btl.Density00; ...
    EXRXcast02_btl.Density00; MID3cast01_btl.Density00; MID4cast14_btl.Density00; MID5cast01_btl.Density00; ...
    WLIScast02_btl.Density00; WLI6cast01_btl.Density00; ...
    ];

PDen = [WLIScast01_btl.Sigma_00; MID4cast01_btl.Sigma_00; EXRXcast01_btl.Sigma_00; MID4cast02_btl.Sigma_00; ...
    MID4cast03_btl.Sigma_00; MID4cast04_btl.Sigma_00; MID4cast05_btl.Sigma_00; MID4cast06_btl.Sigma_00; ...
    MID4cast07_btl.Sigma_00; MID4cast08_btl.Sigma_00; MID4cast09_btl.Sigma_00; MID4cast10_btl.Sigma_00; ...
    MID4cast11_btl.Sigma_00; MID4cast12_btl.Sigma_00; MID4cast13_btl.Sigma_00; EXCR1cast01_btl.Sigma_00; ...
    EXRXcast02_btl.Sigma_00; MID3cast01_btl.Sigma_00; MID4cast14_btl.Sigma_00; MID5cast01_btl.Sigma_00; ...
    WLIScast02_btl.Sigma_00; WLI6cast01_btl.Sigma_00; ...
    ];

O2_umolkg = [WLIScast01_btl.Sbox0Mmkg; MID4cast01_btl.Sbox0Mmkg; EXRXcast01_btl.Sbox0Mmkg; MID4cast02_btl.Sbox0Mmkg; ...
    MID4cast03_btl.Sbox0Mmkg; MID4cast04_btl.Sbox0Mmkg; MID4cast05_btl.Sbox0Mmkg; MID4cast06_btl.Sbox0Mmkg; ...
    MID4cast07_btl.Sbox0Mmkg; MID4cast08_btl.Sbox0Mmkg; MID4cast09_btl.Sbox0Mmkg; MID4cast10_btl.Sbox0Mmkg; ...
    MID4cast11_btl.Sbox0Mmkg; MID4cast12_btl.Sbox0Mmkg; MID4cast13_btl.Sbox0Mmkg; EXCR1cast01_btl.Sbox0Mmkg; ...
    EXRXcast02_btl.Sbox0Mmkg; MID3cast01_btl.Sbox0Mmkg; MID4cast14_btl.Sbox0Mmkg; MID5cast01_btl.Sbox0Mmkg; ...
    WLIScast02_btl.Sbox0Mmkg; WLI6cast01_btl.Sbox0Mmkg; ...
    ];

% O2_umolkg_A = [WLIScast01_btl.O2_umolkg_A; MID4cast01_btl.O2_umolkg_A; EXRXcast01_btl.O2_umolkg_A; MID4cast02_btl.O2_umolkg_A; ...
%     MID4cast03_btl.O2_umolkg_A; MID4cast04_btl.O2_umolkg_A; MID4cast05_btl.O2_umolkg_A; MID4cast06_btl.O2_umolkg_A; ...
%     MID4cast07_btl.O2_umolkg_A; MID4cast08_btl.O2_umolkg_A; MID4cast09_btl.O2_umolkg_A; MID4cast10_btl.O2_umolkg_A; ...
%     MID4cast11_btl.O2_umolkg_A; MID4cast12_btl.O2_umolkg_A; MID4cast13_btl.O2_umolkg_A; EXCR1cast01_btl.O2_umolkg_A; ...
%     EXRXcast02_btl.O2_umolkg_A; MID3cast01_btl.O2_umolkg_A; MID4cast14_btl.O2_umolkg_A; MID5cast01_btl.O2_umolkg_A; ...
%     WLIScast02_btl.O2_umolkg_A; WLI6cast01_btl.O2_umolkg_A; ...
%     ];

Chl = [WLIScast01_btl.Wetstar; MID4cast01_btl.Wetstar; EXRXcast01_btl.Wetstar; MID4cast02_btl.Wetstar; ...
    MID4cast03_btl.Wetstar; MID4cast04_btl.Wetstar; MID4cast05_btl.Wetstar; MID4cast06_btl.Wetstar; ...
    MID4cast07_btl.Wetstar; MID4cast08_btl.Wetstar; MID4cast09_btl.Wetstar; MID4cast10_btl.Wetstar; ...
    MID4cast11_btl.Wetstar; MID4cast12_btl.Wetstar; MID4cast13_btl.Wetstar; EXCR1cast01_btl.Wetstar; ...
    EXRXcast02_btl.Wetstar; MID3cast01_btl.Wetstar; MID4cast14_btl.Wetstar; MID5cast01_btl.Wetstar; ...
    WLIScast02_btl.Wetstar; WLI6cast01_btl.Wetstar; ...
    ];

pH = [WLIScast01_btl.pH; MID4cast01_btl.pH; EXRXcast01_btl.pH; MID4cast02_btl.pH; ...
    MID4cast03_btl.pH; MID4cast04_btl.pH; MID4cast05_btl.pH; MID4cast06_btl.pH; ...
    MID4cast07_btl.pH; MID4cast08_btl.pH; MID4cast09_btl.pH; MID4cast10_btl.pH; ...
    MID4cast11_btl.pH; MID4cast12_btl.pH; MID4cast13_btl.pH; EXCR1cast01_btl.pH; ...
    EXRXcast02_btl.pH; MID3cast01_btl.pH; MID4cast14_btl.pH; MID5cast01_btl.pH; ...
    WLIScast02_btl.pH; WLI6cast01_btl.pH; ...
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

load LISOct2023CastData.mat;
LIS = LISOct2023CastData;
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
LISOct23btl = btl;

save LISOct23btl.mat LISOct23btl;

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

