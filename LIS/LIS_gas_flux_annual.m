cruisedates = [datenum(2023,8,2) datenum(2023,10,19) datenum(2024,5,22)];
fCH4 = [154 133 52];
fCH4med = [69 72 59];

fN2O = [2.5 3.3 4.8];
fN2Omed = [1.9 2.0 5.0];

cdf = [cruisedates cruisedates(1)+365];
fCH4f = [fCH4 fCH4(1)];
fN2Of = [fN2O fN2O(1)];

fCH4medf = [fCH4med fCH4med(1)];
fN2Omedf = [fN2Omed fN2Omed(1)];

%%

cd = [cdf(1):1:cdf(end)];

fCH4i = interp1(cdf,fCH4f,cd);
fN2Oi = interp1(cdf,fN2Of,cd);

fCH4medi = interp1(cdf,fCH4medf,cd);
fN2Omedi = interp1(cdf,fN2Omedf,cd);


mean(fCH4i)
mean(fN2Oi)

mean(fCH4medi)
mean(fN2Omedi)
