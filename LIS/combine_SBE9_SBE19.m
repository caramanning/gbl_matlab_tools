
% this script will allow us to get the in situ T and S for each sample
%load LISAug23btl.mat;
%btl = LISAug23btl;
% load LISOct23btl.mat;
% btl = LISOct23btl;
% 
% load LISOct2023SID.mat;
% SID = LISOct2023SID;

load LISMay24btl.mat;
btl = LISMay24btl;

load LISMay2024SID.mat;
SID = LISMay2024SID;
%%
%deepest_Niskin = 12; % in August niskin 12 is always deepest but this changes in future cruises

% deep cast / cast number / max pressure dbar / max pressure - 0.2 dbar
% deepC = [1	18.86	18.66
% 2	16.688	16.488
% 3	17.389	17.189
% 4	18.287	18.087
% 5	18.241	18.041
% 6	19.128	18.928
% 7	19.447	19.247
% 8	18.972	18.772
% 9	16.909	16.709
% 10	15.843	15.643
% 11	17.691	17.491
% 12	19.115	18.915
% 13	19.638	19.438
% 14	18.783	18.583
% 15	17.184	16.984
% 16	28.372	28.172
% 17	17.816	17.616
% 18	12.915	12.715
% 19	17.057	16.857
% 20	18.292	18.092
% 21	21.389	21.189
% 22	21.919	21.719
% 23	22.783	22.583];

% deep cast data from October
% cast / max pressure dbar / maxpressure - 0.2 dbar
% deepC = [1	19.872	19.672
% 2	16.97	16.77
% 3	16.672	16.472
% 4	16.748	16.548
% 5	17.238	17.038
% 6	18.716	18.516
% 7	19.173	18.973
% 8	17.906	17.706
% 9	17.169	16.969
% 10	16.734	16.534
% 11	17.304	17.104
% 12	18.291	18.091
% 13	19.127	18.927
% 14	18.733	18.533
% 15	17.988	17.788
% 16	31.795	31.595
% 17	22.929	22.729
% 18	12.864	12.664
% 19	17.061	16.861
% 20	16.9	16.7
% 21	19.065	18.865
% 22	20.305	20.105
% ];


% MAY 2024 deep cast data
% cast number / max pressure dbar / max pressure - 0.2 dbar
deepC = [1	18.855	18.655
2	17.559	17.359
3	22.01	21.81
4	18.533	18.333
5	18.946	18.746
6	18.85	18.65
7	17.748	17.548
8	16.79	16.59
9	17.06	16.86
10	18.301	18.101
11	19.037	18.837
12	18.726	18.526
13	17.863	17.663
14	16.887	16.687
15	16.601	16.401
16	29.632	29.432
17	18.784	18.584
18	13.995	13.795
19	18.74	18.54
20	18.731	18.531
21	20.775	20.575
22	21.268	21.068
23	23.292	23.092
24	30.866	30.666
25	27.233	27.033
];
%%
% add in the deep cast Station, Cast, Niskin, Pressure
sl = numel(btl.Cast);
si= sl+1:sl+numel(deepC(:,1));

btl.Cast(si) = deepC(:,1);
btl.Niskin(si) = 0;
btl.P(si) = deepC(:,3);


% now run a script to add in the sample IDs
btl.SID1 = nan .* btl.Cast;
btl.SID2 = nan .* btl.Cast;

for i = 1:length(btl.Cast)
    A = find(SID.Cast == btl.Cast(i));
    B = find(SID.Niskin == btl.Niskin(i));
    C = intersect(A,B);
    if numel(C) == 2
        btl.SID1(i) = C(1);
        btl.SID2(i) = C(2);
    elseif numel(C) > 2
        disp(['more than 2 bottles matching at index ', num2str(i)]);
    end;
end;
%%

% add in the Hand CTD data
sl = numel(btl.Cast);
uC = unique(btl.Cast);
C = nan.*uC;
pdeep = nan.*uC;

% find rows with niskin 0
Ni0 = find(btl.Niskin == 0);
E = nan(numel(pdeep(:,1)),1);
%%
for i = 1:length(Ni0)
%for i = 4
       
    %A = find(btl.Cast==btl.Cast(Ni0(i)) & ~isnan(btl.SID1) & btl.Niskin>0);
    A = ~isnan(btl.SID1);
    B = btl.Niskin>0; 
    C = btl.Cast==btl.Cast(Ni0(i));
    D = A & B & C;

    if ~isempty(find(D))
        e = find(D);
        f = find(max(btl.Depth(e)));
        E(i) = e(f);
    end;

    %B = find(btl.Niskin==deepest_Niskin);
  %  B = find(max(btl.Depth(A)));
  %  C(i) = A(B);
  %  D = find(deepC(:,1)==uC(i));
  %  pdeep(i) = deepC(D,3);
end;


%%
% now run the script to fill in all of the data for the deep casts for which we actually
% collected samples
A = ~isnan(E);

btl.Station(si(A)) = btl.Station(E(A));

btl.S(si(A)) = btl.S(E(A));
btl.T(si(A)) = btl.T(E(A));
%btl.P(si(A)) = pdeep;
btl.Depth(si(A)) = sw_dpth(btl.P(si(A)),40.97);
btl.PDen(si(A)) = sw_pden(btl.S(si(A)),btl.T(si(A)),btl.P(si(A)),0)-1000;
btl.PTemp(si(A)) = sw_ptmp(btl.S(si(A)),btl.T(si(A)),btl.P(si(A)),0);

btl.Dens(si(A)) = btl.Dens(E(A));
btl.O2_umolkg(si(A)) = btl.O2_umolkg(E(A));
%btl.O2_umolkg_A(si(A)) = btl.O2_umolkg_A(E(A));
btl.Chl(si(A)) = btl.Chl(E(A));
btl.pH(si(A)) = btl.pH(E(A));
btl.Lat(si(A)) = btl.Lat(E(A));
btl.Lon(si(A)) = btl.Lon(E(A));
btl.datetime(si(A)) = btl.datetime(E(A));
btl.CH4N2Ocast(si(A)) = btl.CH4N2Ocast(E(A));

%%
no_sample_casts=find(btl.Depth==0);
btl([no_sample_casts],:) = [];


%%
LISMay24btl_combo = btl;
save LISMay24btl_combo.mat LISMay24btl_combo;

%%

%%
LISOct23btl_combo = btl;
save LISOct23btl_combo.mat LISOct23btl_combo;

%%
LISAug23btl_combo = btl;
save LISAug23btl_combo.mat LISAug23btl_combo;