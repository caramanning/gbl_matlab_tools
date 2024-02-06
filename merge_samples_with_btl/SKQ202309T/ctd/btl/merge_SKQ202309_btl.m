%merge_SKQ202309_btl

% load CTD files one by one and merge into a master table

for ca = 1:7
    % load CTD file and add in a column for cast number
    if ca==1
        load skq202309t001.mat;
        skq202309t001.Cast = ca.*ones(size(skq202309t001.Bottle));   

    elseif ca==2
        load skq202309t002.mat;
        skq202309t002.Cast = ca.*ones(size(skq202309t002.Bottle));

    elseif ca==3
        load skq202309t003.mat;
        skq202309t003.Cast = ca.*ones(size(skq202309t003.Bottle));

    elseif ca==4
        load skq202309t004.mat;
        skq202309t004.Cast = ca.*ones(size(skq202309t004.Bottle));   

    elseif ca==5
        load skq202309t005.mat;
        skq202309t005.Cast = ca.*ones(size(skq202309t005.Bottle)); 

    elseif ca==6
        load skq202309t006.mat;
        skq202309t006.Cast = ca.*ones(size(skq202309t006.Bottle)); 

    elseif ca==7
        load skq202309t007.mat; 
        skq202309t007.Cast = ca.*ones(size(skq202309t007.Bottle));             
    end;

end;

skq202309 = [skq202309t001
    skq202309t002
    skq202309t003
    skq202309t004
    skq202309t005
    skq202309t006
    skq202309t007];
%%
save skq202309.mat skq202309;
