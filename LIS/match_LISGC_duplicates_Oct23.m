load LISOct23GC.mat
LIS = LISOct23GC;

u_SID = unique(LIS.SID); % find unique sample IDs
u_SID(isnan(u_SID)) = []; % remove non-numeric values

LIS.Cast(LIS.Cast=='CA2') = 'Ca2';
LIS.Cast(LIS.Cast=='CA3') = 'Ca3';
LIS.Cast(LIS.Cast=='CA5') = 'Ca5';
LIS.Cast(LIS.Cast=='CA7') = 'Ca7';
LIS.Cast(LIS.Cast=='CA9') = 'Ca9';
LIS.Cast(LIS.Cast=='CA11') = 'Ca11';
LIS.Cast(LIS.Cast=='CA13') = 'Ca13';
LIS.Cast(LIS.Cast=='CA15') = 'Ca15';
LIS.Cast(LIS.Cast=='CA16') = 'Ca16';
LIS.Cast(LIS.Cast=='CA17') = 'Ca17';
LIS.Cast(LIS.Cast=='CA18') = 'Ca18';
LIS.Cast(LIS.Cast=='CA19') = 'Ca19';
LIS.Cast(LIS.Cast=='CA20') = 'Ca20';
LIS.Cast(LIS.Cast=='CA21') = 'Ca21';
LIS.Cast(LIS.Cast=='CA22') = 'Ca22';
LIS.Cast(LIS.Cast=='HAND') = 'hand';

LIS.Niskin(LIS.Niskin=='NI1') = 'Ni1';
LIS.Niskin(LIS.Niskin=='NI2') = 'Ni2';
LIS.Niskin(LIS.Niskin=='NI3') = 'Ni3';
LIS.Niskin(LIS.Niskin=='NI4') = 'Ni4';
LIS.Niskin(LIS.Niskin=='NI5') = 'Ni5';
LIS.Niskin(LIS.Niskin=='NI6') = 'Ni6';
LIS.Niskin(LIS.Niskin=='NI7') = 'Ni7';
LIS.Niskin(LIS.Niskin=='NI8') = 'Ni8';
LIS.Niskin(LIS.Niskin=='NI9') = 'Ni9';
LIS.Niskin(LIS.Niskin=='NI10') = 'Ni10';
LIS.Niskin(LIS.Niskin=='NI11') = 'Ni11';
LIS.Niskin(LIS.Niskin=='NI12') = 'Ni12';
LIS.Niskin(LIS.Niskin=='HAND') = 'hand';



%%

% run a loop to find the duplicate sample IDs
match_SID = nan(numel(u_SID),2);nan

for i = 1:length(u_SID)
    if ~isnan(u_SID(i)) % check that the SID is a non-nan value
        i_u_SID=find(LIS.SID == u_SID(i)); % find row of SID in table
        match_Cast=find(LIS.Cast(i_u_SID) == LIS.Cast); % matching Cast
        match_Niskin = find(LIS.Niskin(i_u_SID) == LIS.Niskin); % matching Niskin
        match_CastNiskin = intersect(match_Cast,match_Niskin); %matching both Cast and Niskin
        match_CastNiskin = sort(match_CastNiskin); % sort values low to high
    
        % now save the values into a row
        if numel(match_CastNiskin)==1 % if there is only one matching value
            match_SID(i,1) = LIS.SID(match_CastNiskin); % 
        elseif numel(match_CastNiskin)==2 % if there is two matching values
            match_SID(i,:) = LIS.SID(match_CastNiskin);
        elseif numel(match_CastNiskin)>2 % if there is two matching values
            disp(['three matching values for index ',num2str(i)]);   % this will come out if there are NaN values      
        end;
    end;
end;

%%
% c2 is the unique values in column 1
% i2 is the index of the first unique value of every pair of replicates
[c2, i2] = unique(match_SID(:,1));

u_SID = match_SID(i2,:);
u_SID(isnan(u_SID(:,1)),:) = []; % remove nan rows


%% calculate mean & std
a=LIS;
usn=u_SID; % unique sample numbers

% preallocate variables
%Cruise = nan(length(usn),1);
%Cast = nan(length(usn),1);
%Station = nan(length(usn),1);
%Niskin = nan(length(usn),1);
mean_CH4_nM = nan(length(usn),1);
std_CH4_nM = nan(length(usn),1);
mean_N2O_nM = nan(length(usn),1);
std_N2O_nM = nan(length(usn),1);

for i=1:length(usn)
    A=find(a.SID==usn(i,1)); % get index for duplicate A
    B=find(a.SID==usn(i,2)); % get index for duplicate B
    Cruise(i) = a.Cruise(A);
    Cast(i) = a.Cast(A);
    Station(i) = a.Station(A);
    Niskin(i) = a.Niskin(A);
    mean_CH4_nM(i) = nanmean(a.CH4nM([A,B]));
    std_CH4_nM(i)=nanstd(a.CH4nM([A,B]));
    mean_N2O_nM(i) = nanmean(a.N2OnM([A,B]));
    std_N2O_nM(i)=nanstd(a.N2OnM([A,B]));

end

Cruise = Cruise';
a_m = table(Cruise);
a_m.Cast = Cast';
a_m.Niskin = Niskin';
a_m.mean_CH4_nM = mean_CH4_nM;
a_m.std_CH4_nM = std_CH4_nM;
a_m.mean_N2O_nM = mean_N2O_nM;
a_m.std_N2O_nM = std_N2O_nM;

%%
LISOct23_CH4N2O = a_m;

save LISOct23_CH4N2O.mat LISOct23_CH4N2O;

%%
%a.
%std_ch4=a.std_ch4';
%%
for i=1:numel(usn)
    A=find(a.sample_num==usn(i));
    a.std_n2o(A)=nanstd(a.n2o(A));
end
a.std_n2o=a.std_n2o';

for i=1:numel(usn)
    A=find(a.sample_num==usn(i));
    a.ch4(A)=nanmean(a.ch4(A));
end

for i=1:numel(usn)
    A=find(a.sample_num==usn(i));
    a.n2o(A)=nanmean(a.n2o(A));
end


%%
Cruise = LIS.Cruise(dup_SID(:,1));
LIS_CH4N2O = table(Cruise);
%%
LIS_CH4N2O.Cast = LIS.Cast(dup_SID(:,1));
LIS_CH4N2O.Niskin = LIS.Niskin(dup_SID(:,1));
LIS_CH4N2O.CH4_nM = mean(LIS.CH4nM(dup_SID),2);
LIS_CH4N2O.CH4_nM_std = std(LIS.CH4nM(dup_SID,:),0,1);

%%
std(LIS.CH4nM(dup_SID(:,1),:),1,2)

%%

% now figure out code to plot all of the 