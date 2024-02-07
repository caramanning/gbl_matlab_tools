% sample indices from the run
SID = [384
364
321
372
766
313
667
672
684
683
677
288
344
373
405
316
282
437
515
479
426
519
396
428
407
423
400
421];

%%

SID = [778
756
734
781
765
681
782
664
679
663
350
391
394
336
354
362
170
351
359
439
279
342
314
293
369
312
];
%%
run_data = table(SID); % create a table
%%
load IEP_SID_TS.mat

% now find the row in IEP_SID_TS that matches the SID of each sample in the
% run

% note: some of the sample IDs around 766 seem to be matching with multiple
% niskin bottles - need to check how the data was sorted
IEP_SID_TS_index = nan.*run_data.SID;
for i = 1:length(run_data.SID)
    IEP_SID_TS_index_i = find(IEP_SID_TS.SID == run_data.SID(i));
    if numel(IEP_SID_TS_index_i) ~= 1
       disp(['for index ', num2str(i), ' zero or multiple values matching SID'])
    else
       IEP_SID_TS_index(i) = IEP_SID_TS_index_i;
    end;
end
%%
run_data.T = IEP_SID_TS.T(IEP_SID_TS_index);
run_data.S = IEP_SID_TS.S(IEP_SID_TS_index);
run_data.dens = IEP_SID_TS.dens(IEP_SID_TS_index);
run_data.Station = IEP_SID_TS.Station(IEP_SID_TS_index);
run_data.Niskin = IEP_SID_TS.Niskin(IEP_SID_TS_index);

%%
writetable(run_data,'run_data.csv')

%% now calculate the exact density for the air-eq samples

T_eq = [20.95
20.95
20.95
5
5
5
];

S = 0.*T_eq;

P = 0.*T_eq;

dens = sw_dens(S,T_eq,P);
