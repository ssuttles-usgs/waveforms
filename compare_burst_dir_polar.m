clear all ; close all ; clc ; 
% see where negative skewness is 
%load('/media/taran/DATADRIVE2/Obs_data/matfiles/skewness_steve.mat','Su_skewness','dn')
%load('C:\Users\tkalra\Desktop\Observationaldata\matfiles\matfiles\skewness_steve.mat',....
%    'Su_skewness','dn'); 
%load('/media/taran/DATADRIVE2/Obs_data/matfiles/skewness_steve.mat','Su_skewness','dn')
load('C:\Users\ssuttles\data\FireIsland\analysis\site3\mat\skewness_steve.mat',....
    'Su_skewness','dn'); 
wh=fullfile('C:\Users\ssuttles\data\FireIsland\analysis\Taran\9921whp-cal.nc')  ; 
%wh=fullfile('C:\Users\tkalra\Desktop\Observationaldata\data_FI\9921whp-cal.nc')  ; 
%wh=fullfile('/media/taran/DATADRIVE2/Obs_data/data_netcdf/9921whp-cal.nc'); % statistics filename
% convert work horse surface data of wave energy spectra to get ubr and
% Tbr..

netcdf_load(wh)
nt1=1614; nt2=1621;

d_spec(:,:,:)=double(dspec(1,1,:,:,:));
 v_spec(:,:)=double(vspec(1,1,:,:));
 p_spec(:,:)=double(pspec(1,1,:,:)); 
 s_spec(:,:)=double(sspec(1,1,:,:)); 
band_width=0.015625  ;
f=double(frequency(:,1));
df=band_width ;

isave=0 ;

%d_spec=(d_spec*0.001).^2; % converted to sq.m/Hz from sq.mm/Hz
%d_spec=d_spec*0.001*0.001; % converted sq.m/Hz/degrees
d_spec=d_spec*1e-6; %convert from mm^2/Hz to m^2/Hz
ndirs=length(direction);
d_spec=d_spec/ndirs; %normalize to direction to get in m^2/Hz/degree
% v_spec=(v_spec*0.001).^2; 
% p_spec=(p_spec*0.001).^2 ;

count_post=1; count_neg=1 ; 

%for t=nt1:nt2;
  d_spec(d_spec<0)=0.0;
  d_spec(d_spec>1e12)=0.0;
%end
% 
%dspec_3d(:,:,:)=double(d_spec(1,1,:,:,:)); 


dspec1=d_spec(:,:,nt1);
dspec2=d_spec(:,:,nt2); 
% % nt=1166; % corresponding to valentine's day
% % nt=1424; % correspondign to highest negative Su
% % nt=1166;  % corresponding to highest positive Su

for t=1:64
   direction_2d(:,t)=direction(:);
end     
for tt=1:360
   frequency_2d(tt,:)=frequency(:);
end 


figure(1)
set(gcf,'position',[200 50 900 700])
ax0=axes('position',[0.15 0.25 0.4 0.5],'box','on');
S.dirs=direction;
S.freqs=frequency;
S.dspec=dspec1;
tstr=sprintf('Burst %d',nt1)
ca=[-2.5 -0.5];
hp0=plrplt_dspec_ax(ax0,S,tstr,ca,1)


ax1=axes('position',[0.55 0.25 0.4 0.5],'box','on');
S.dirs=direction;
S.freqs=frequency;
S.dspec=dspec2;
tstr=sprintf('Burst %d',nt2);
hp1=plrplt_dspec_ax(ax1,S,tstr,ca,0)

colormap(jet)

% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 8])
% print('-dpng','-r300','pngfiles\compare_burst_dspec.png')
% 
figure(2)
set(gcf,'position',[200 50 900 700])
ax0=axes('position',[0.15 0.25 0.4 0.5],'box','on');
S.dirs=direction;
S.freqs=frequency;
S.dspec=dspec1;
tstr=sprintf('Burst %d',nt1);
ca=[-2.5 -0.5];
hp0=plrplt_dspec_ax_freq(ax0,S,tstr,ca,1)


ax1=axes('position',[0.55 0.25 0.4 0.5],'box','on');
S.dirs=direction;
S.freqs=frequency;
S.dspec=dspec2;
tstr=sprintf('Burst %d',nt2);
hp1=plrplt_dspec_ax_freq(ax1,S,tstr,ca,0)

colormap(jet)
% 
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 8])
% print('-dpng','-r300','pngfiles\compare_dspec_burst.png')