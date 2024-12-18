%load LISAug23GC.mat
load LISAug23GC_v2.mat
%%
% REMOVE THE ROWS WITH BAD DATA using the note 'delete'
toDelete = find(LISAug23GC.note=='delete');

%toDelete = Tnew.Age < 30;
LISAug23GC(toDelete,:) = [];

%%
u_SID = unique(LISAug23GC.SID); % find unique sample IDs
u_SID(isnan(u_SID)) = []; % remove non-numeric values
%%

% run a loop to find the duplicate sample IDs
match_SID = nan(numel(u_SID),2);

for i = 1:length(u_SID)
    if ~isnan(u_SID(i)) % check that the SID is a non-nan value
        i_u_SID=find(LISAug23GC.SID == u_SID(i)); % find row of SID in table
        match_Cast=find(LISAug23GC.Cast(i_u_SID) == LISAug23GC.Cast); % matching Cast
        match_Niskin = find(LISAug23GC.Niskin(i_u_SID) == LISAug23GC.Niskin); % matching Niskin
        match_CastNiskin = intersect(match_Cast,match_Niskin); %matching both Cast and Niskin
        match_CastNiskin = sort(match_CastNiskin); % sort values low to high
    
        % now save the values into a row
        if numel(match_CastNiskin)==1 % if there is only one matching value
            match_SID(i,1) = LISAug23GC.SID(match_CastNiskin); % 
        elseif numel(match_CastNiskin)==2 % if there is two matching values
            match_SID(i,:) = LISAug23GC.SID(match_CastNiskin);
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
a=LISAug23GC;
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
    mean_CH4_nM(i) = mean(a.CH4_nM([A,B]),'omitnan');
    std_CH4_nM(i) = std(a.CH4_nM([A,B]),'omitnan');
    mean_N2O_nM(i) = mean(a.N2O_nM([A,B]),'omitnan');
    std_N2O_nM(i) = std(a.N2O_nM([A,B]),'omitnan');

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
LISAug23_CH4N2O = a_m;

save LISAug23_CH4N2O.mat LISAug23_CH4N2O;

% %%
% %a.
% %std_ch4=a.std_ch4';
% %%
% for i=1:numel(usn)
%     A=find(a.sample_num==usn(i));
%     a.std_n2o(A)=nanstd(a.n2o(A));
% end
% a.std_n2o=a.std_n2o';
% 
% for i=1:numel(usn)
%     A=find(a.sample_num==usn(i));
%     a.ch4(A)=nanmean(a.ch4(A));
% end
% 
% for i=1:numel(usn)
%     A=find(a.sample_num==usn(i));
%     a.n2o(A)=nanmean(a.n2o(A));
% end
% 
% 
% %%
% Cruise = LISAug23GC.Cruise(dup_SID(:,1));
% LISAug23_CH4N2O = table(Cruise);
% %%
% LISAug23_CH4N2O.Cast = LISAug23GC.Cast(dup_SID(:,1));
% LISAug23_CH4N2O.Niskin = LISAug23GC.Niskin(dup_SID(:,1));
% LISAug23_CH4N2O.CH4_nM = mean(LISAug23GC.CH4nM(dup_SID),2);
% LISAug23_CH4N2O.CH4_nM_std = std(LISAug23GC.CH4nM(dup_SID,:),0,1);
% 
% %%
% std(LISAug23GC.CH4nM(dup_SID(:,1),:),1,2)
% 
% %%
% 
% % now figure out code to plot all of the 