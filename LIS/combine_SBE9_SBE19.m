
% this script will allow us to get the in situ T and S for each sample
load LISAug23btl.mat;

% deep cast / cast number / pressure dbar / pressure - 0.2 m dbar
deepC = [1	18.86	18.66
2	16.688	16.488
3	17.389	17.189
4	18.287	18.087
5	18.241	18.041
6	19.128	18.928
7	19.447	19.247
8	18.972	18.772
9	16.909	16.709
10	15.843	15.643
11	17.691	17.491
12	19.115	18.915
13	19.638	19.438
14	18.783	18.583
15	17.184	16.984
16	28.372	28.172
17	17.816	17.616
18	12.915	12.715
19	17.057	16.857
20	18.292	18.092
21	21.389	21.189
22	21.919	21.719
23	22.783	22.583];


btl = LISAug23btl;

% add in the Hand CTD data
sl = numel(btl.Cast);
uC = unique(btl.Cast);
C = nan.*uC;
pdeep = nan.*uC;

for i = 1:length(uC)
    A = find(btl.Cast==uC(i));
    B = find(btl.Niskin==12); % Niskin 12 is always the deepest in August;
    C(i) = intersect(A,B)
    D = find(deepC(:,1)==uC(i));
    pdeep(i) = deepC(D,3);
   % btl.Cast(sl+i) = btl.Cast(C)
   % btl.
end;

si= sl+1:sl+numel(uC);
btl.Station(si) = btl.Station(C);
btl.Cast(si) = btl.Cast(C);
btl.Niskin(si) = 0;

btl.S(si) = btl.S(C);
btl.T(si) = btl.T(C);
btl.P(si) = pdeep;
btl.Depth(si) = sw_dpth(btl.P(si),40.97);
btl.PDen(si) = sw_pden(btl.S(si),btl.T(si),btl.P(si),0)-1000;
btl.PTemp(si) = sw_ptmp(btl.S(si),btl.T(si),btl.P(si),0);

btl.Dens(si) = btl.Dens(C);
btl.O2_umolkg(si) = btl.O2_umolkg(C);
btl.O2_umolkg_A(si) = btl.O2_umolkg_A(C);
btl.Chl(si) = btl.Chl(C);
btl.pH(si) = btl.pH(C);
btl.Lat(si) = btl.Lat(C);
btl.Lon(si) = btl.Lon(C);
btl.datetime(si) = btl.datetime(C);
btl.CH4N2Ocast(si) = btl.CH4N2Ocast(C);

LISAug23btl_combo = btl;
save LISAug23btl_combo.mat LISAug23btl_combo;