% T-S plot demo
% ABOUT: T-S plots of different years illustrating the different water
% masses
clear,clc,close all
base_dir='archived_data';
cd(base_dir)
load laurier2015.mat
load laurier2016.mat
load laurier2017.mat
load laurier2018.mat
load jois2015.mat
load jois2016.mat
load jois2017.mat
load jois2018.mat;
load amundsen2015.mat
load amundsen2016.mat
load amundsen2017.mat
load amundsen2018.mat

cd ..\
%%

% figure out appropriate ranges for the plots
theta=[-3 10];
s=[22 36];
smin=min(s)-0.01.*min(s);
smax=max(s)+0.01.*max(s);
thetamin=min(theta)-0.1*max(theta);
thetamax=max(theta)+0.1*max(theta);
xdim=round((smax-smin)./0.1+1);
ydim=round((thetamax-thetamin)+1);
dens=zeros(ydim,xdim);
thetai=((1:ydim)-1)*1+thetamin;
si=((1:xdim)-1)*0.1+smin;
disp(xdim);disp(ydim);
for j=1:ydim
    for i=1:xdim
        dens(j,i)=sw_dens(si(i),thetai(j),0);
    end
end

dens=dens-1000;

%%
% create structs of different variables for iteration
temp.y5=cat(1,laurier2015.temp,jois2015.temp);
temp.y5=cat(1,temp.y5,amundsen2015.temp);
temp.y6=cat(1,laurier2016.temp,jois2016.temp);
temp.y6=cat(1,temp.y6,amundsen2016.temp);
temp.y7=cat(1,laurier2017.temp,jois2017.temp);
temp.y7=cat(1,temp.y7,amundsen2017.temp);
temp.y8=cat(1,laurier2018.temp,jois2018.temp);


sal.y5=cat(1,laurier2015.sal,jois2015.sal);
sal.y5=cat(1,sal.y5,amundsen2015.sal);
sal.y6=cat(1,laurier2016.sal,jois2016.sal);
sal.y6=cat(1,sal.y6,amundsen2016.sal);
sal.y7=cat(1,laurier2017.sal,jois2017.sal);
sal.y7=cat(1,sal.y7,amundsen2017.sal);
sal.y8=cat(1,laurier2018.sal,jois2018.sal);

% pressure
press.y5=cat(1,laurier2015.press,jois2015.press,amundsen2015.press);
press.y6=cat(1,laurier2016.press,jois2016.press,amundsen2016.press);
press.y7=cat(1,laurier2017.press,jois2017.press,amundsen2017.press);
press.y8=cat(1,laurier2018.press,jois2018.press);


pt.y5 = sw_ptmp(sal.y5,temp.y5,press.y5,zeros(length(sal.y5),1));
pt.y6 = sw_ptmp(sal.y6,temp.y6,press.y6,zeros(length(sal.y6),1));
pt.y7 = sw_ptmp(sal.y7,temp.y7,press.y7,zeros(length(sal.y7),1));
pt.y8 = sw_ptmp(sal.y8,temp.y8,press.y8,zeros(length(sal.y8),1));

ch4.y5=cat(1,laurier2015.ch4_nmolkg,jois2015.ch4_nmolkg);
ch4.y5=cat(1,ch4.y5,amundsen2015.ch4_nmolkg);
ch4.y6=cat(1,laurier2016.ch4_nmolkg,jois2016.ch4_nmolkg);
ch4.y6=cat(1,ch4.y6,amundsen2016.ch4_nmolkg);
ch4.y7=cat(1,laurier2017.ch4_nmolkg,jois2017.ch4_nmolkg);
ch4.y7=cat(1,ch4.y7,amundsen2017.ch4_nmolkg);
ch4.y8=cat(1,laurier2018.ch4_nmolkg,jois2018.ch4_nmolkg);


n2o.y5=cat(1,laurier2015.n2o_nmolkg,jois2015.n2o_nmolkg);
n2o.y5=cat(1,n2o.y5,amundsen2015.n2o_nmolkg);
n2o.y6=cat(1,laurier2016.n2o_nmolkg,jois2016.n2o_nmolkg);
n2o.y6=cat(1,n2o.y6,amundsen2016.n2o_nmolkg);
n2o.y7=cat(1,laurier2017.n2o_nmolkg,jois2017.n2o_nmolkg);
n2o.y7=cat(1,n2o.y7,amundsen2017.n2o_nmolkg);
n2o.y8=cat(1,laurier2018.n2o_nmolkg,jois2018.n2o_nmolkg);


lon.y5=cat(1,laurier2015.lon,jois2015.lon,amundsen2015.lon);
lon.y6=cat(1,laurier2016.lon,jois2016.lon,amundsen2016.lon);
lon.y7=cat(1,laurier2017.lon,jois2017.lon,amundsen2017.lon);
lon.y8=cat(1,laurier2018.lon,jois2018.lon);

lon.all = cat(1,lon.y5,lon.y6,lon.y7,lon.y8);
pt.all = cat(1,pt.y5,pt.y6,pt.y7,pt.y8);
sal.all = cat(1,sal.y5,sal.y6,sal.y7,sal.y8);
ch4.all = cat(1,ch4.y5,ch4.y6,ch4.y7,ch4.y8);
n2o.all = cat(1,n2o.y5,n2o.y6,n2o.y7,n2o.y8);
%%

% set boundaries for east, center, west based on longitude
e = find(lon.all>=-112);
c = find(lon.all<-112&lon.all>=-158);
w = find(lon.all<-158);

lon.west = lon.all(w);
pt.west = pt.all(w);
sal.west = sal.all(w);
n2o.west = n2o.all(w);
ch4.west = ch4.all(w);

lon.cent = lon.all(c);
pt.cent = pt.all(c);
sal.cent = sal.all(c);
n2o.cent = n2o.all(c);
ch4.cent = ch4.all(c);

lon.east = lon.all(e);
pt.east = pt.all(e);
sal.east = sal.all(e);
n2o.east = n2o.all(e);
ch4.east = ch4.all(e);

% set plot properties
yt = [0 4 8]; % yaxis ticks
xt = [25 30 35]; %x axis ticks
yl = [-2.5 10]; %yaxls limit
xl = [24 35.5]; %x axis limit
%cach4 = prctile(ch4.all,[3 97]); % colorbar range for ch4
cach4 = [1 20]; % colorbar range for ch4
%can2o = prctile(n2o.all,[3 97]); %colorbar range for
can2o = [12 19]; % colorbar range for n2o
rp1 = [30 -3 3.5 3]; % position
rp2 = [30 0 3.5 7]; % position
fs=16; % font size
ms=30; % marker size
lwr=2; % line width

nr = 2; % number of rows
nc = 3;  % number of columns
fig=figure(1)
clf;
sp=tight_subplot(nr,nc,[.05 .03],[.15 .1],[.08 .04]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 6]);


% set(gcf,'color','w');
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperPosition', [0 0 12 8]);
set(gca,'linewidth',2);
subplot(sp(1))

hold on; box on;
set(gca,'linewidth',2);

[c,h]=contour(si,thetai,dens,'--','color',[0.5 0.5 0.5]);
clabel(c,h,'LabelSpacing',1000);

    [~,a]=sort(ch4.west);    
    scatter(sal.west(a),pt.west(a),ms,ch4.west(a),'filled') 
    
    rectangle('Position',rp1,'linewidth',lwr);
%    rectangle('Position',rp2,'linewidth',1.5);   
        
    ylim(yl);
    xlim(xl);
    caxis(cach4);
    set(gca,'ytick',yt);
    set(gca,'xtick',xt);
    set(gca,'yticklabel',yt);
   % xlabel(['Salinity (PSS)'])
    ylabel(['Temperature [^oC]'])
    set(gca,'fontsize',fs); 
    set(gca,'tickdir','out')
    set(gca,'ticklength',[0.02 0.02]);
    hold on

    
subplot(sp(2))
hold on; box on;
[c,h]=contour(si,thetai,dens,'--','color',[0.5 0.5 0.5]);
set(gca,'linewidth',2);
%clabel(c,h,'LabelSpacing',1000);
    
    [~,a]=sort(ch4.cent);    
    scatter(sal.cent(a),pt.cent(a),ms,ch4.cent(a),'filled') 
    rectangle('Position',rp1,'linewidth',lwr);
%    rectangle('Position',rp2,'linewidth',1.5);      
    ylim(yl)
    xlim(xl)
    caxis(cach4)
    set(gca,'ytick',yt);
    set(gca,'xtick',xt);
    set(gca,'xticklabel',[]);
    set(gca,'yticklabel',[]);
 %   xlabel(['Salinity ', char(8240)])
   % ylabel(['Temperature ',char(176), 'C'])
    set(gca,'fontsize',fs); 
    set(gca,'tickdir','out')
    set(gca,'ticklength',[0.02 0.02]);
    hold on    

subplot(sp(3))
hold on; box on;
set(gca,'linewidth',2);
[c,h]=contour(si,thetai,dens,'--','color',[0.5 0.5 0.5]);
%clabel(c,h,'LabelSpacing',1000);
    
    [~,a]=sort(ch4.east);    
    scatter(sal.east(a),pt.east(a),ms,ch4.east(a),'filled')     
    rectangle('Position',rp1,'linewidth',lwr);
 %   rectangle('Position',rp2,'linewidth',1.5);      
    ylim(yl)
    xlim(xl)
    caxis(cach4)
    set(gca,'ytick',yt);
    set(gca,'xtick',xt);
    %xlabel(['Salinity ', char(8240)])
 %   ylabel(['Temperature ',char(176), 'C'])
    set(gca,'fontsize',fs); 
    set(gca,'tickdir','out')
    set(gca,'ticklength',[0.02 0.02]);
        set(gca,'xticklabel',[]);
    set(gca,'yticklabel',[]);
    hold on

  % OPTIONAL code for adding colorbar; it makes the right plot smaller
  %  hcb=colorbar; % colorbar
  %  hcb.Label.String = 'CH_4 (nmol/kg)'; % colorbar label
    
subplot(sp(4))
hold on; box on;
set(gca,'linewidth',2);
[c,h]=contour(si,thetai,dens,'--','color',[0.5 0.5 0.5]);
clabel(c,h,'LabelSpacing',1000);

    [~,a]=sort(n2o.west);    
    scatter(sal.west(a),pt.west(a),ms,n2o.west(a),'filled') 
    rectangle('Position',rp1,'linewidth',lwr);
%    rectangle('Position',rp2,'linewidth',1.5);          
    ylim(yl)
    xlim(xl)
    caxis(can2o)
    set(gca,'ytick',yt);
    set(gca,'xtick',xt);
    set(gca,'xticklabel',xt);
    set(gca,'yticklabel',yt);
    xlabel(['Salinity [PSS]'])
    ylabel(['Temperature [^oC]'])
    set(gca,'fontsize',fs); 
    set(gca,'tickdir','out')
    set(gca,'ticklength',[0.02 0.02]);
    hold on
    
subplot(sp(5))
hold on; box on;
set(gca,'linewidth',2);
[c,h]=contour(si,thetai,dens,'--','color',[0.5 0.5 0.5]);
%clabel(c,h,'LabelSpacing',1000);
    
    [~,a]=sort(n2o.cent);    
    scatter(sal.cent(a),pt.cent(a),ms,n2o.cent(a),'filled') 
    rectangle('Position',rp1,'linewidth',lwr);
 %   rectangle('Position',rp2,'linewidth',1.5);          
    ylim(yl)
    xlim(xl)
    caxis(can2o)
    set(gca,'ytick',yt);
    set(gca,'xtick',xt);
        set(gca,'xticklabel',xt);
    xlabel(['Salinity [PSS]'])
 %   ylabel(['Temperature ',char(176), 'C'])
    set(gca,'fontsize',fs); 
    set(gca,'tickdir','out')
    set(gca,'ticklength',[0.02 0.02]);
    hold on
    
subplot(sp(6))
hold on; box on;
set(gca,'linewidth',2);
[c,h]=contour(si,thetai,dens,'--','color',[0.5 0.5 0.5]);
%clabel(c,h,'LabelSpacing',1000);
    
    [~,a]=sort(n2o.east);    
    scatter(sal.east(a),pt.east(a),ms,n2o.east(a),'filled') 
    rectangle('Position',rp1,'linewidth',lwr);
 %   rectangle('Position',rp2,'linewidth',1.5);          
    ylim(yl)
    xlim(xl)
    caxis(can2o)
    set(gca,'ytick',yt);
    set(gca,'xtick',xt);
    xlabel(['Salinity [PSS]'])
        set(gca,'xticklabel',xt);
   % ylabel(['Temperature ',char(176), 'C'])
    set(gca,'fontsize',fs); 
    set(gca,'tickdir','out')
    set(gca,'ticklength',[0.02 0.02]);
    hold on    
  % OPTIONAL code for adding colorbar; it makes the right plot smaller
  %  hcb=colorbar; % colorbar
  %  hcb.Label.String = 'N_2O (nmol/kg)'; % colorbar label
    
    
    wysiwyg;

  % version of plots without colorbar, every subplot is the same size  
  % print -dpdf -r300 20200324_TSplot.pdf;
  % print -dpng -r300 20200324_TSplot.png;    
    
  % version of plots with colorbar, subplot on right is smaller, just using
  % to print the colorbar
   % print -dpdf -r300 20200324_TSplot_colorbar.pdf;
   % print -dpng -r300 20200324_TSplot_colorbar.png;
    
