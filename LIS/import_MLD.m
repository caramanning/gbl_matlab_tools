load LISAug23_CH4N2O_CTD.mat;

sl= ["EXRX-cast01"
"EXRX-cast02"
"EXR1-cast01"
"MID3-cast01"
"MID4-cast01"
"MID4-cast02"
"MID4-cast03"
"MID4-cast04"
"MID4-cast05"
"MID4-cast06"
"MID4-cast07"
"MID4-cast08"
"MID4-cast09"
"MID4-cast10"
"MID4-cast11"
"MID4-cast12"
"MID4-cast13"
"MID4-cast14"
"MID5-cast01"
"WLI6-cast01"
"WLI7-cast01"
"WLIS-cast01"
"WLIS-cast02"];

mld = [8
5
5
8
4
5
4
4
4
2
3
5
6
8
7
10
5
7
8
5
4
3
5];

d=LISAug23_CH4N2O_CTD;
d.mld = nan.*d.CH4_mean_nmolkg;
for i = 1:length(d.Station)
    A=strcmp(sl,string(d.Station(i)));
    Ai = find(A);
    if ~isempty(Ai)
        d.mld(i) = mld(Ai);
    end;
end;

LISAug23_CH4N2O_CTD = d;
%%
save LISAug23_CH4N2O_CTD.mat LISAug23_CH4N2O_CTD;

%%

sl = ["EXR1-cast01"
"EXRX-cast01"
"EXRX-cast02"
"MID3-cast01"
"MID4-cast01"
"MID4-cast02"
"MID4-cast03"
"MID4-cast04"
"MID4-cast05"
"MID4-cast06"
"MID4-cast07"
"MID4-cast08"
"MID4-cast09"
"MID4-cast10"
"MID4-cast11"
"MID4-cast12"
"MID4-cast13"
"MID4-cast14"
"MID5-cast01"
"WLI6-cast01"
"WLIS-cast01"
"WLIS-cast02"];

mld = [7
14
9
7
6
8
8
7
5
7
6
7
9
9
9
8
7
9
7
7
7
7
];

load LISOct23_CH4N2O_CTD.mat;
%%
d=LISOct23_CH4N2O_CTD;
d.mld = nan.*d.CH4_mean_nmolkg;
for i = 1:length(d.Station)
    A=strcmp(sl,string(d.Station(i)));
    Ai = find(A);
    if ~isempty(Ai)
        d.mld(i) = mld(Ai);
    end;
end;

LISOct23_CH4N2O_CTD = d;

save LISOct23_CH4N2O_CTD.mat LISOct23_CH4N2O_CTD;

%%

mld = [4
3
3
4
7
4
3
4
4
2
4
3
6
4
4
7
5
4
3
5
5
3
2
2
4
];

sl = ["ARTG-cast01"
"CLIS-cast01"
"EXR1-cast01"
"EXRX-cast01"
"EXRX-cast02"
"MID3-cast01"
"MID4-cast01"
"MID4-cast02"
"MID4-cast03"
"MID4-cast04"
"MID4-cast05"
"MID4-cast06"
"MID4-cast07"
"MID4-cast08"
"MID4-cast09"
"MID4-cast10"
"MID4-cast11"
"MID4-cast12"
"MID4-cast13"
"MID4-cast14"
"MID5-cast01"
"WLI6-cast01"
"WLI7-cast01"
"WLIS-cast01"
"WLIS-cast02"
];


load LISMay24_CH4N2O_CTD.mat;
%%
d=LISMay24_CH4N2O_CTD;
d.mld = nan.*d.CH4_mean_nmolkg;
for i = 1:length(d.Station)
    A=strcmp(sl,string(d.Station(i)));
    Ai = find(A);
    if ~isempty(Ai)
        d.mld(i) = mld(Ai);
    end;
end;

LISMay24_CH4N2O_CTD = d;
%%
save LISMay24_CH4N2O_CTD.mat LISMay24_CH4N2O_CTD;
