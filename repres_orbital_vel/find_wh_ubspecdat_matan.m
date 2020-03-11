% written by Tarandeep S Kalra with help from Steve Suttles
% it calls the output of ADV based Hrmsu, ubr, Tr and compares
% with the workhorse data to get the work horse based wave orbital velocity and representative
% time period 

clear all ; close all ; clc; 
 wh=fullfile('C:\Users\ssuttles\data\Matanzas\proc\released\11101whVp-cal.nc'); % statistics filename
% convert work horse surface data of wave energy spectra to get ubr and
% Tbr..

netcdf_load(wh)
nt1=1; nt2=length(burst); 

 jtb = double(time)+double(time2)/(3600*24*1000); 
      dnsb = datestr(datenum(gregorian(jtb))); 

Hs(:)=squeeze(wh_4061(1,1,:));
h_1d(:)=double(hght_18(1,1,:));
s(:,:)=double(vspec(1,1,:,:));  


f=double(frequency(:,1));
df=diff(f);
df=df(1);

isave=0 ; 

s=(s*0.001).^2; % converted to sq.m/Hz from sq.mm/Hz 

for t=nt1:nt2; 
 % for ft=1:length(f)
  %if(s<-1e8)
     s(s<0)=0.0;
     s(s>1e12)=0.0;
  %end 
    [ubr(t),Tbr(t)]=ubspecdat(squeeze(h_1d(t)),s(:,t)',f(:,1)',df);
    Hs_wh(t)=Hs(t);   
    Hs_wh(Hs_wh>100)=0.0;
    %if (Hs_wh(i)>100);
    %    Hs(i)=0.0;
    %end
    
end 
uhat_wh=ubr;
Tr_wh=Tbr;


save('C:\Users\ssuttles\data\Matanzas\proc\mat\vspec_uhat_tr_matan.mat','uhat_wh','Tr_wh')

