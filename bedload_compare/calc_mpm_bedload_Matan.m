clear all ; close all ; clc; 

nt1=1; nt2= 1893; 
% 
% 
load('C:\Users\ssuttles\data\Matanzas\proc\mat\puv_proc_Matan_iwaves','UBS','PUV','depth')
%%       
gamma=8; 
% make sure same d50 used for the vandera formula
d50 = 0.23e-3 ; %0.2e-3;
d90 = 1.3*d50;
rho0 = 1025.0;
g = 9.81 ;            % m/s2
rhos= 2650.0;

s=rhos/rho0; 
smgd=(s-1.0)*g*d50;
osmgd=1.0/smgd; 
%theta_cr=0.0 ; 

smgd3=(s-1)*g*d50.^3; 
theta_cr=theta_cr_calc(d50,rhos)  ; 

% calculate for each burst 
nsamps=length(UBS(1).ur) ; 

%jt = time+time2/(3600*24*1000);
%dn=j2dn(time,time2);
% theta_cr =0.   ; % WRONG ASSUMPTION


 for t=1:nt2
    disp('.')
    if rem(t,10)==0
        disp(t)
    end
    Su_skewness(t)=mean(UBS(t).ur.^3)/(std(UBS(t).ur)).^3;
   [qb_measured(t)]=func_calc_mpm (UBS(t).ur,PUV(t).omegar,......
                        PUV(t).ubr,d50,nsamps,osmgd,theta_cr,.......
                        gamma,smgd3) ; 
 end
  
% 
% plot(qb_measured)
%  ylim([-0.05 0.15])
%  xlim([dn(nt1) dn(nt1+300)]);
%  datetick('x',2) % '
% % hold on
% plot(bedldx(1,1:301)*3000)
clf
dt=3600 ; 

bedload_measured=cumtrapz(qb_measured) 

save('C:\Users\ssuttles\data\Matanzas\proc\mat\mpm_only_ss_Matan.mat','Su_skewness','qb_measured','bedload_measured')
