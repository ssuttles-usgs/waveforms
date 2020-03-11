clear all ; close all ; clc; 
% Use the workhorse data with UBSPECFUN code that is using JONSWAP spectrum
% Code written by Tarandeep S. Kalra and Chris Sherwood

% Enter sediment information in meters
d50 = 0.4e-3 ; %0.2e-3;
d90 = 1.3*d50;

% Near bottom current data
% This is constant
deg2rad=pi/180.0; 
umag_curr=0.0; %.2814;% abs(0.0);
phi_curwave=0.0;% 79.92*deg2rad ;% 0.0*deg2rad;
% Zref is a reference height for computing current friction factor 
Zref=0.04 ;
% delta is the reference height at which current velocity is computed (Wave boundary layer thickness) 
delta=0.2;

nt1=1; nt2= 2044;
%  load('skewness_steve.mat','Su_skewness','ur_maj_rot','vr_min_rot',........
%             'Hrmsu','Ubr','Ur','dn','jtb_rec','ur_cube')

% LOAD the observational data from workhorse from Fire Island 
url='/media/taran/DATADRIVE2/Obs_data/data_netcdf/9921whp-cal.nc'
netcdf_load(url)
Hs(:)=squeeze(wh_4061(1,1,:));
Td(:)=squeeze(wp_peak(1,1,:));
h(:)=squeeze(hght_18(1,1,:)); % extract depth; 

 % h=depth; % Depth from ADV  
% %ntime=end 
 for i=1:length(Hs)
     if (Hs(i)>100);
        Hs(i)=0.0;
     end
     if (Td(i)>30); 
         Td(i)=0.0;
     end 
     [uhat(i),Tbav(i)]=ubspecfun(Hs(i),Td(i),h(i) ); 
 end
% 

% Question- Would need to use peak wave period and not bottom wave period
waveavgd_stress_term=1; 
surface_wave=1;  
current_timeperiod==0 ; 
% NOTE THAT I AM USING TBOT AND NOT PEAK WAVE PERIOD 
% 
%t_end=2000 ;
 for i=nt1:nt2
   [bedldx(i), bedldy(i), Ur(i), R(i), Beta(i)]=vandera_function_workhorse(i, Hs(i), Td(i), h(i), d50, d90, .....
                                                  umag_curr, phi_curwave, uhat(i), .....
                                                  Zref, delta, waveavgd_stress_term, surface_wave);
 end 
 save('vandera_bedld_workhorse_ubspecfun.mat','bedldx','R','Beta','Ur')
%  figure(1)
% plot(uhat(1:t_end),bedldx(1:t_end),'*--')
% xlabel('uhat')
% ylabel('bedld in x')
% print('-dpng','-r100','fig1.png')
% 
%  figure(2)
%  plot(Hs(1:t_end),Ur(1:t_end),'.')
% xlabel('Hs')
% ylabel('Ursell number') 
% print('-dpng','-r100','fig2.png')
% 
%  figure(3)
%  plot(Hs(1:t_end),R(1:t_end),'.')
% ylabel('R')
% xlabel('Hs') 
% print('-dpng','-r100','fig3.png')
% 
% % 
% figure(4)
% plot(Tbot(1:t_end),'r.')
% hold on
% plot(Td(1:t_end),'.')
% ylabel('Tbot')
% xlabel('Td')
% print('-dpng','-r100','fig5.png')

% figure(4)
% plot(Hs(1:t_end),bedldy(1:t_end),'.')
% xlabel('Hs')
% ylabel('bedld in y')
