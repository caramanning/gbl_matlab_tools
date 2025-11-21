data = tdmsread('20250826_1824.tdms');
dataT = data{1};

%%

figure(1)
clf; hold on;
plot(dataT.TimeStamps,dataT.AI0);
plot(dataT.TimeStamps,dataT.AI1);
plot(dataT.TimeStamps,dataT.AI2);
legend('electronics','inside box','ambient room','location','northeast');