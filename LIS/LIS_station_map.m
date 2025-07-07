c5 = [1 0 0];
c6 = [138/255 43/255 215/255];
c7 = [0 1 0];
c8 = [1 0.5 0];

c1 = [253 204 138]./255;
c2 = [252 141 89]./255;
c3 = [227 74 51]./255;
c4 = [179 0 0]./255;

bathymetry2 = [0.152941182	0.227450982	0.372549027
0.148895115	0.237306774	0.389148265
0.144849062	0.247162566	0.405747503
0.140802994	0.257018358	0.422346711
0.136756927	0.266874164	0.438945949
0.132710874	0.276729941	0.455545187
0.128664806	0.286585748	0.472144425
0.126723826	0.29281801	0.481043875
0.124782853	0.299050272	0.489943326
0.122841872	0.305282533	0.498842776
0.120900892	0.311514765	0.507742226
0.118959911	0.317747027	0.516641676
0.117018938	0.323979288	0.525541127
0.115077958	0.33021155	0.534440577
0.113136977	0.336443812	0.543340027
0.111195996	0.342676073	0.552239478
0.109255023	0.348908305	0.561138928
0.107314043	0.355140567	0.570038378
0.105373062	0.361372828	0.578937829
0.103432089	0.36760509	0.587837279
0.101491109	0.373837352	0.596736729
0.099550128	0.380069613	0.605636179
0.097609147	0.386301875	0.61453563
0.095668174	0.392534107	0.62343514
0.093727194	0.398766369	0.63233459
0.091786213	0.40499863	0.64123404
0.089845233	0.411230892	0.650133491
0.08790426	0.417463154	0.659032941
0.085963279	0.423695415	0.667932391
0.084022298	0.429927677	0.676831841
0.082081325	0.436159909	0.685731292
0.080140345	0.44239217	0.694630742
0.078199364	0.448624432	0.703530192
0.076258384	0.454856694	0.712429643
0.074317411	0.461088955	0.721329093
0.07237643	0.467321217	0.730228543
0.070435449	0.473553449	0.739127994
0.068494469	0.479785711	0.748027444
0.066553496	0.486017972	0.756926894
0.064612515	0.492250234	0.765826344
0.084135473	0.503921568	0.770944774
0.103658438	0.515592933	0.776063144
0.123181395	0.527264237	0.781181574
0.142704353	0.538935602	0.786300004
0.162227318	0.550606906	0.791418374
0.181750283	0.562278271	0.796536803
0.201273233	0.573949575	0.801655233
0.220796198	0.58562094	0.806773603
0.240319163	0.597292244	0.811892033
0.259842128	0.608963609	0.817010462
0.279365093	0.620634913	0.822128892
0.298888028	0.632306278	0.827247262
0.318410993	0.643977582	0.832365692
0.337933958	0.655648947	0.837484121
0.357456923	0.667320251	0.842602491
0.376979887	0.678991616	0.847720921
0.396502852	0.69066292	0.852839351
0.416025817	0.702334285	0.857957721
0.435548753	0.714005589	0.86307615
0.455071718	0.725676954	0.86819458
0.474594682	0.737348258	0.87331295
0.494117647	0.749019623	0.87843138
0.570588231	0.813725471	0.896078467
0.647058845	0.87843138	0.913725495
0.735294133	0.908823535	0.935294122
0.867647067	0.954411767	0.967647061];
%%
%load bathymetry_colormap;

%bathymetry_ud = flipud(bathymetry);
%save bathymetry_ud_colormap bathymetry_ud;


ms=9;
fs=16;
figure(10)
clf; hold on;
set(gcf,'color','w');
set(gca,'fontsize',fs);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 10 5]);

lonlim = [-74 -71.8];
latlim = [40.5 41.35];

m_proj('lambert','lon',lonlim,'lat',latlim);
m_grid('linestyle','none','tickdirection','out','linewidth',1);

[cs,ch]=m_etopo_Ddrive('contourf',[-80:2:-8 -4 -3 -2 -1],'color','none');
m_gshhs_f('patch',[.6 .6 .6],'edgecol','none');
%m_gshhs_f('patch',[.6 .6 .6],'edgecol','none');
%m_gshhs('lc','patch','r');  % Low resolution filled coastline
%m_gshhs('fb1');             % Full resolution national borders
%m_gshhs('fr','patch',[0.6 0.6 0.6],'edgecol','none');              % Intermediate resolution rivers
%m_gshhs('fr');




s_ll = [-73.73683,40.8710
    -73.72917,40.8830
    -73.69200,40.90533
    -73.65533,40.92333
    -73.614333,40.94000
    -73.582167,40.957167
    -73.55767,40.96167];

% s_ll = [-73.72863	40.88312333
% -73.72685	40.88502
% -73.69131	40.90469
% -73.65494	40.92345833
% -73.61414	40.940435
% -73.58215	40.95661667
% -73.55767	40.96254333];

m_line(s_ll(:,1),s_ll(:,2),'color','r');

for i = 1:7
m_line(s_ll(i,1), s_ll(i,2),'marker','d','markersize',ms-5,'color','k','markerfacecol','r'); % EXR1
end;

% m_line(-73.7369, 40.8710,'marker','.','markersize',ms-5,'color','k','markerfacecol','g'); % EXR1 Aug
% m_line(-73.7360, 40.8710,'marker','.','markersize',ms-5,'color','k','markerfacecol','g'); % EXR1 Aug
% m_line(-73.7368, 40.8720,'marker','.','markersize',ms-5,'color','k','markerfacecol','g'); % EXR1 Aug
% 
% m_line(-73.5814, 40.9558,'marker','.','markersize',ms-5,'color','k','markerfacecol','g'); % EXR1 Aug
% m_line(-73.5840, 40.9572,'marker','.','markersize',ms-5,'color','k','markerfacecol','g'); % EXR1 Aug
% m_line(-73.5822, 40.9572,'marker','.','markersize',ms-5,'color','k','markerfacecol','g'); % EXR1 Aug



% m_line(-73.72917,40.8830,'marker','d','markersize',ms-5,'color','k','markerfacecol','r'); % EXRX
% 
% m_line(-73.69200,40.90533,'marker','d','markersize',ms-5,'color','k','markerfacecol','r'); % MID3
% 
% m_line(-73.65533,40.92333,'marker','d','markersize',ms-5,'color','k','markerfacecol','r'); % MID4
% 
% m_line(-73.614333,40.94000,'marker','d','markersize',ms-5,'color','k','markerfacecol','r'); % MID5
% 
% m_line(-73.582167,40.957167,'marker','d','markersize',ms-5,'color','k','markerfacecol','r'); % WLIS
% 
% m_line(-73.55767,40.96167,'marker','d','markersize',ms-5,'color','k','markerfacecol','r'); % WLI6

%m_line(-73.5355,40.96683,'marker','o','markersize',ms-5,'color','k','markerfacecol','k'); % WLI7

m_line(-73.28683,41.01117,'marker','o','markersize',ms-5,'color','k','markerfacecol','k'); % ARTG

m_line(-72.65550,41.13833,'marker','o','markersize',ms-5,'color','k','markerfacecol','k'); % CLIS

lt = -73.9;
% m_text(lt,41.3,'stations east to west','color','k','fontsize',fs-2);
% m_text(lt,41.3-0.05,'WLI6','color','r','fontsize',fs-2);
% m_text(lt,41.3-2*0.05,'WLIS','color','r','fontsize',fs-2);
% m_text(lt,41.3-3*0.05,'MID5','color','r','fontsize',fs-2);
% m_text(lt,41.3-4*0.05,'MID4','color','r','fontsize',fs-2);
% m_text(lt,41.3-5*0.05,'MID3','color','r','fontsize',fs-2);
% m_text(lt,41.3-6*0.05,'EXRX','color','r','fontsize',fs-2);
% m_text(lt,41.3-7*0.05,'EXR1','color','r','fontsize',fs-2);
% m_text(-73.25,41.05,'ARTG','color','k','fontsize',fs-2);
% m_text(-72.65,41.18,'CLIS','color','k','fontsize',fs-2);
% %m_text(-71.5, 40.9,'depth (m)','vertical','color','k','fontsize',fs);

b=flipud(m_colmap('blues'));
colormap(b);

colorbar('location','westoutside');
wysiwyg;

print -dpng -r300 20250707_LIS_Map.png;

print -depsc -r300 20250707_LIS_Map.eps;

print(gcf,'-depsc','20250707_LIS_Map.eps');
epsclean('20250707_LIS_Map.eps','20250707_LIS_Map_epsclean.eps');


%%

ms=9;
fs=16;
figure(12)
clf; hold on;
set(gcf,'color','w');
set(gca,'fontsize',fs);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 5 5]);

lonlim = [-73.8 -73.5];
latlim = [40.75 41];

m_proj('lambert','lon',lonlim,'lat',latlim);
m_grid('linestyle','none','tickdirection','out','linewidth',1);

%[cs,ch]=m_etopo_Ddrive('contourf',[-80:2:-8 -4 -3 -2 -1],'color','none');
%m_gshhs_f('patch',[.6 .6 .6],'edgecol','none');
m_gshhs_f('patch',[.6 .6 .6],'edgecol','none');
%m_gshhs('lc','patch','r');  % Low resolution filled coastline
%m_gshhs('fb1');             % Full resolution national borders
%m_gshhs('fr','patch',[0.6 0.6 0.6],'edgecol','none');              % Intermediate resolution rivers
%m_gshhs('fr');


m_line(s_ll(:,1),s_ll(:,2),'color','r','linewidth',2);

for i = 1:7
    m_line(s_ll(i,1), s_ll(i,2),'marker','d','markersize',ms,'color','k','markerfacecol','r'); % EXR1
end;
% 
% m_line(-73.7369, 40.8710,'marker','.','markersize',ms-5,'color','k','markerfacecol','g'); % EXR1 Aug
% m_line(-73.7360, 40.8710,'marker','.','markersize',ms-5,'color','k','markerfacecol','g'); % EXR1 Aug
% m_line(-73.7368, 40.8720,'marker','.','markersize',ms-5,'color','k','markerfacecol','g'); % EXR1 Aug
% 
% m_line(-73.5814, 40.9558,'marker','.','markersize',ms-5,'color','k','markerfacecol','g'); % EXR1 Aug
% m_line(-73.5840, 40.9572,'marker','.','markersize',ms-5,'color','k','markerfacecol','g'); % EXR1 Aug
% m_line(-73.5822, 40.9572,'marker','.','markersize',ms-5,'color','k','markerfacecol','g'); % EXR1 Aug

lt = -73.9;
% m_text(lt,41.3,'stations east to west','color','k','fontsize',fs-2);
% m_text(lt,41.3-0.05,'WLI6','color','r','fontsize',fs-2);
% m_text(lt,41.3-2*0.05,'WLIS','color','r','fontsize',fs-2);
% m_text(lt,41.3-3*0.05,'MID5','color','r','fontsize',fs-2);
% m_text(lt,41.3-4*0.05,'MID4','color','r','fontsize',fs-2);
% m_text(lt,41.3-5*0.05,'MID3','color','r','fontsize',fs-2);
% m_text(lt,41.3-6*0.05,'EXRX','color','r','fontsize',fs-2);
% m_text(lt,41.3-7*0.05,'EXR1','color','r','fontsize',fs-2);
% m_text(-73.25,41.05,'ARTG','color','k','fontsize',fs-2);
% m_text(-72.65,41.18,'CLIS','color','k','fontsize',fs-2);
%m_text(-71.5, 40.9,'depth (m)','vertical','color','k','fontsize',fs);

%b=flipud(m_colmap('blues'));
%colormap(b);

%colorbar;
wysiwyg;

print -dpng -r300 20250707_LIS_Map_inset.png;

print -depsc -r300 20250707_LIS_Map_inset.eps;

print(gcf,'-depsc','20250707_LIS_Map_inset.eps');

epsclean('20250707_LIS_Map_inset.eps','20250707_LIS_Map_epsclean_inset.eps');


%%

ms=9;
fs=16;
figure(10)
clf; hold on;
set(gcf,'color','w');
set(gca,'fontsize',fs);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 10 5]);

lonlim = [-73.9 -73.5];
latlim = [40.8 41.1];

m_proj('lambert','lon',lonlim,'lat',latlim);
m_grid('linestyle','none','tickdirection','out','linewidth',1);

[cs,ch]=m_etopo_Ddrive('contourf',[-80:1:-8 -4 -3 -2 -1],'color','none');
m_gshhs_f('patch',[.6 .6 .6],'edgecol','none');
%m_gshhs_f('patch',[.6 .6 .6],'edgecol','none');
%m_gshhs('lc','patch','r');  % Low resolution filled coastline
%m_gshhs('fb1');             % Full resolution national borders
m_gshhs('ir','patch',[0.6 0.6 0.6],'edgecol','none');              % Intermediate resolution rivers
%m_gshhs('fr');

s_ll = [-73.73683,40.87200
    -73.72917,40.8830
    -73.69200,40.90533
    -73.65533,40.92333
    -73.614333,40.94000
    -73.582167,40.957167
    -73.55767,40.96167];

m_line(s_ll(:,1),s_ll(:,2),'color','r');

m_line(-73.73683,40.87200,'marker','d','markersize',ms-5,'color','k','markerfacecol','r'); % EXR1

m_line(-73.72917,40.8830,'marker','d','markersize',ms-5,'color','k','markerfacecol','r'); % EXRX

m_line(-73.69200,40.90533,'marker','d','markersize',ms-5,'color','k','markerfacecol','r'); % MID3

m_line(-73.65533,40.92333,'marker','d','markersize',ms-5,'color','k','markerfacecol','r'); % MID4

m_line(-73.614333,40.94000,'marker','d','markersize',ms-5,'color','k','markerfacecol','r'); % MID5

m_line(-73.582167,40.957167,'marker','d','markersize',ms-5,'color','k','markerfacecol','r'); % WLIS

m_line(-73.55767,40.96167,'marker','d','markersize',ms-5,'color','k','markerfacecol','r'); % WLI6

%m_line(-73.5355,40.96683,'marker','o','markersize',ms-5,'color','k','markerfacecol','k'); % WLI7

%m_line(-73.28683,41.01117,'marker','o','markersize',ms-5,'color','k','markerfacecol','k'); % ARTG

%m_line(-72.65550,41.13833,'marker','o','markersize',ms-5,'color','k','markerfacecol','k'); % CLIS



%b=flipud(bathymetry);
%colormap(b);
b=flipud(m_colmap('blues'));
colormap(b);
%caxis([-4000 000]);



%m_proj('lambert','lon',lonlim,'lat',latlim);
%m_grid('tickdirection','out','linewidth',1);

%m_etopo4('contourf',[-3000 -2250 -2000 -1750 -1500 -1250 -1000 -750 -500 -250 -100 -50 -20 -5 0],'edgecolor','none');


%m_proj('mercator','lon',[-90 -48],'lat',[60 82])
%m_grid('tickdirection','out','linewidth',3,'backgroundcolor',bathymetry(1,:))
%m_gshhs_l('patch',[.6 .6 .6],'edgecol','none')

%set(findobj('tag','m_grid_color'),'facecolor','none') 


%b=flipud(bathymetry);

%colormap(m_colmap('blues'));  
%caxis([-2500 000]);
 % m_line(AN1902_CTD.Lon(c(1)),AN1902_CTD.Lat(c(1)),'marker','d','markersize',ms,'color','r','markerfacecol','r'); 
% m_line(xx,yy,'markersize',ms-5,'color','r','markerfacecol','r','linewi',5); 

%   BB2
%%%%%  m_line(-66.9933,72.7537,'marker','d','markersize',ms-5,'color','r','markerfacecol','r'); 

 % m_line(-56.803,60.686,'marker','d','markersize',ms-5,'color','r','markerfacecol','r'); 

 
%colormap(b)
colorbar;
wysiwyg;
%%
print -dpng -r300 LIS_Map.png;

print -depsc -r300 LIS_Map.eps;

print(gcf,'-depsc','LIS_Map.eps');
epsclean('LIS_Map.eps','LIS_Map_epsclean.eps');
