% 1) get the list of sample IDs from a GC run xlsx sheet and copy below in the
% SID field
% 2) run the code below to extract the in situ T, S and Density
% 3) transfer those back to the excel sheet
% 4) ensure the summary sheet is updated for the new concentrations

%load LISAug23btl_combo.mat;
%btl = LISAug23btl_combo;

load LISOct23btl_combo.mat;
btl = LISOct23btl_combo;

load LISMay24btl_combo.mat;
btl = LISMay24btl_combo;

% find sample IDs for a list
SID = [11
37
29
7
47
49
27
45
1
];

A = nan.*SID;

for i = 1:length(SID)
    a = find(btl.SID1 == SID(i) | btl.SID2 == SID(i));
    if numel(a) == 1
        A(i) = a;
    else
        disp(['error at index ', num2str(i)]);
    end;
end;

TSD = nan(numel(SID),3);

TSD(:,1) = btl.T(A);
TSD(:,2) = btl.S(A);
TSD(:,3) = btl.Dens(A);
