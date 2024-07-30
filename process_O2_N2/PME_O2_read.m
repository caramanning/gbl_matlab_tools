%PME_read_O2
%directory of file
%fdir = "deployment Aug 5\7450-383325\2022-08-12 202900Z.txt";
%fdir = "test deployment August 1 to 3\2022-08-02 172900Z.txt";
fdir = "lab calibration Aug 4\7450-364333\2022-08-04 012300Z.txt";
%fdir = "deployment Aug 12\2022-08-04 013900Z.txt";

% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 5);

% Specify range and delimiter
opts.DataLines = [4, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["unixtime", "BV", "T", "DO", "Q"];
opts.VariableTypes = ["double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
PME_O2 = readtable(fdir, opts);

% Clear temporary variables
clear opts


PME_O2.datetime = datetime(PME_O2.unixtime, 'ConvertFrom', 'posixtime');

figure(1)
clf; 
subplot(2,1,1)
hold on;
plot(PME_O2.datetime,PME_O2.T);
datetick;

subplot(2,1,2)
hold on;
plot(PME_O2.datetime,PME_O2.DO);
datetick;

save PME_O2_2022-08-04-012300Z.mat PME_O2;
