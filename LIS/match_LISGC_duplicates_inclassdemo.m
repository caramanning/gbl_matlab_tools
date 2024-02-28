% load data
load LISAug23GC.mat

% find unique sample IDs (their value)
u_SID = unique(LISAug23GC.SID);
u_SID(isnan(u_SID)) = [];

match_SID = nan(length(u_SID),2); %create variable to store matching SID
%%
% for each unique sample ID, find the other sample IDs that match the cast
% and niskin

for i = 1:length(u_SID)
%for i = 158
    i_u_SID=find(LISAug23GC.SID == u_SID(i)); % find row of SID in table
    match_Cast=find(LISAug23GC.Cast(i_u_SID) == LISAug23GC.Cast); % matching Cast
    match_Niskin = find(LISAug23GC.Niskin(i_u_SID) == LISAug23GC.Niskin); % matching Niskin
    match_CastNiskin = intersect(match_Cast,match_Niskin); %matching both Cast and Niskin
    match_CastNiskin = sort(match_CastNiskin); % sort values low to high, check whether this is needed

    % now save the values for match_CastNiskin to a variable
    if numel(match_CastNiskin)==1
        match_SID(i,1) = LISAug23GC.SID(match_CastNiskin);
    elseif numel(match_CastNiskin)==2 % if there are 2 matching values
        match_SID(i,:) = LISAug23GC.SID(match_CastNiskin);
    elseif numel(match_CastNiskin)>2
        disp(['more than 2 matching values for SID = ', num2str(u_SID(i))]);
    end

end;

%%

[c2, i2] = unique(match_SID(:,1));

u_SID = match_SID(i2,:);
u_SID(isnan(u_SID(:,1)),:) = []; % remove nan rows


%% calculate mean & std
a=LISAug23GC;
usn=u_SID; % unique sample numbers
%%
% preallocate variables
mean_CH4_nM = nan(length(usn),1);
std_CH4_nM = nan(length(usn),1);
mean_N2O_nM = nan(length(usn),1);
std_N2O_nM = nan(length(usn),1);

for i=1:length(usn)
%for i = 76
    A=find(a.SID==usn(i,1)); % get index for duplicate A
    B=find(a.SID==usn(i,2)); % get index for duplicate B
    Cruise(i) = a.Cruise(A);
    Cast(i) = a.Cast(A);
    Station(i) = a.Station(A);
    Niskin(i) = a.Niskin(A);
    mean_CH4_nM(i) = mean(a.CH4nM([A,B]));
    std_CH4_nM(i) = std(a.CH4nM([A,B]));
    mean_N2O_nM(i) = mean(a.N2OnM([A,B]));
    std_N2O_nM(i) = std(a.N2OnM([A,B]));
end;

Cruise = Cruise';
a_m = table(Cruise);

%%
a_m.Cast = Cast';
a_m.Niskin = Niskin';
a_m.mean_CH4_nM = mean_CH4_nM;
a_m.std_CH4_nM = std_CH4_nM;
a_m.mean_N2O_nM = mean_N2O_nM;
a_m.std_N2O_nM = std_N2O_nM;

%%

LISAug23_CH4N2O = a_m;

save LISAug23_CH4N2O_2.mat LISAug23_CH4N2O;


%%

A = find(a_m.Cast=='Ca22');

a_m.mean_CH4_nM(A)
