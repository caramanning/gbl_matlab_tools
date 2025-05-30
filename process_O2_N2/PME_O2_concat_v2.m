% 
%fdir = "deployment Sep 12\7450-364333\mat";
%cd 'deployment Sep 12\7450-364333\mat'
%fname = strcat('PME_cc_364333.mat'); %file name to export

cd 'deployment Sep 12\7450-383325\mat'
fname = strcat('PME_cc_383325.mat'); %file name to export

d=dir('*.mat');  % get the list of files
%PME_O2_cc=[]; % start w/ an empty array
%A = [];
for i=1:length(d)
    load(d(i).name)
    if i == 1
        A = PME_O2;
    else
        A = [A; PME_O2];
    end
end
PME_O2_cc = A;

save(fname,'PME_O2_cc');

% 
% %%
% PME_O2_cc=[PME_O2_cc; load(d(i).name)];   % read/concatenate into x
% end
% 
% save PME_concat_364333.mat PME_O2_cc
% 
% %%
% %his script does a manual concatenation of the daily files
% % you will have to edit to update the filenames
% % plan to update code to automate in the future
% 
% % load PME_O2_2022-08-01-172900Z.mat
% % A = PME_O2;
% % 
% % load PME_O2_2022-08-02-172900Z.mat
% % A = [A; PME_O2];
% % 
% % PME_O2_cc = A;
% % 
% % save PME_O2_cc.mat PME_O2_cc;
% 
% load PME_O2_2022-08-12-153300Z.mat
% A = PME_O2;
% 
% 
% load PME_O2_2022-08-13-153300Z.mat
% A = [A; PME_O2];
% 
% load PME_O2_2022-08-14-153300Z.mat
% A = [A; PME_O2];
% 
% load PME_O2_2022-08-15-153300Z.mat
% A = [A; PME_O2];
% 
% load PME_O2_2022-08-16-153300Z.mat
% A = [A; PME_O2];
% 
% load PME_O2_2022-08-17-153300Z.mat
% A = [A; PME_O2];
%  
% load PME_O2_2022-08-18-153300Z.mat
% A = [A; PME_O2];
% 
% load PME_O2_2022-08-19-153300Z.mat
% A = [A; PME_O2];
% 
% 
% % load PME_O2_2022-08-04-202900Z.mat
% % A = PME_O2;
% % 
% % load PME_O2_2022-08-05-202900Z.mat
% % A = [A; PME_O2];
% % 
% % load PME_O2_2022-08-06-202900Z.mat
% % A = [A; PME_O2];
% % 
% % load PME_O2_2022-08-07-202900Z.mat
% % A = [A; PME_O2];
% % 
% % load PME_O2_2022-08-08-202900Z.mat
% % A = [A; PME_O2];
% % 
% % load PME_O2_2022-08-09-202900Z.mat
% % A = [A; PME_O2];
% % 
% % load PME_O2_2022-08-10-202900Z.mat
% % A = [A; PME_O2];
% % 
% % load PME_O2_2022-08-11-202900Z.mat
% % A = [A; PME_O2];
% 
% PME_O2_cc = A;
% 
% save PME_O2_cc.mat PME_O2_cc;

 