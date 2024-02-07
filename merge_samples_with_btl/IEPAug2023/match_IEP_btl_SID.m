% load MAT files created with import_IEP_SID.m and import_IEP_btl.m
load IEP_btl.mat;
load IEP_SID.mat;

match_sb = nan.*IEP_SID.Niskin; % create blank variable for storing matching indices

% find the index (row) in IEP_btl that matches every sample in IEP_SID
for i = 1:length(IEP_SID.Station)  
    match_s = find(IEP_btl.Grid == IEP_SID.Station(i)); % matching station name
    match_b = find(IEP_btl.Bottle == IEP_SID.Niskin(i)); % matching bottle
    match_sb(i) = intersect(match_s,match_b); % matching both station and bottle
end
%%
% calculate potential temperature, which will be used to calculate
% solubility later
IEP_btl.PotentialT = sw_ptmp(IEP_btl.SalinityPSU,IEP_btl.TemperatureITS90DegC,IEP_btl.Pressuredb,0);

% add pressure, depth, T, S, density variables to IEP_SID
IEP_SID.P = IEP_btl.Pressuredb(match_sb);
IEP_SID.depth = IEP_btl.Depthm(match_sb);
IEP_SID.T = IEP_btl.TemperatureITS90DegC(match_sb);
IEP_SID.S = IEP_btl.SalinityPSU(match_sb);
IEP_SID.dens = sw_dens(IEP_SID.S,IEP_SID.T,IEP_SID.P);

IEP_SID_TS = IEP_SID;

save IEP_SID_TS.mat IEP_SID_TS;

% CODE TO ADD LATER: CONVERSION TO umolkg
% O2 conversion from ml/l to umol/kg
% O2_umolkg = (O2_mll * 44.660)/((1000 + SIG0)/1000)




