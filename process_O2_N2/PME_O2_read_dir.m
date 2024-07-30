%PME_O2_read_dir
%directory of file
fdir = "deployment Sep 12\7450-383325\";
cd(fdir)

if ~exist('mat', 'dir')
       mkdir('mat')
end

foutdir = strcat(fdir,'mat\');
dirs = dir;
%fdir = "deployment Sep 12\7450-383325\"; % return to main directory
%dirs = dir(fdir); %directory structure

% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 5);
disp(opts)

% Specify range and delimiter
opts.DataLines = [4, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["unixtime", "BV", "T", "DO", "Q"];
opts.VariableTypes = ["double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

for i = 1:length(dirs)
    % if the file name is more than 2 characters then import and save
    if length(dirs(i).name)>3
        fpath = strcat(fdir,dirs(i).name); %file path
        fname = strcat('PME_O2_',dirs(i).name(1:10),'-',dirs(i).name(12:end-4),'.mat'); %file name to export
        % Import and then save the data
        PME_O2 = readtable(fpath, opts);
        PME_O2.datetime = datetime(PME_O2.unixtime, 'ConvertFrom', 'posixtime');
        floc = strcat('mat\',fname); 
        save(floc,'PME_O2'); 
       % save fname PME_O2;
    end
end

% Clear temporary variables
clear opts

figure(1)
clf; 
subplot(2,1,1)
hold on;
%title(PME_O2.datetime(1))
plot(PME_O2.datetime,PME_O2.T);
datetick;

subplot(2,1,2)
hold on;
plot(PME_O2.datetime,PME_O2.DO);
datetick;
