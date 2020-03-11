%find corrected depth for Matanzas site A (mooring# 1110), using corrected pressure and brange from downlooking Signature1000, and MicroCat TC logger at 1109 (surface) and 1110 (near-bed).
clear;
%files needed
sigfn='C:\Users\ssuttles\data\Matanzas\proc\released\11102sigs_echo1-cal.nc'; %corrected pressure (P_1ac) and brange
mcfn1='C:\Users\ssuttles\data\Matanzas\proc\released\11091mc-a.nc'; %seacat logger with temp and salinity data
mcfn2='C:\Users\ssuttles\data\Matanzas\proc\released\11105mc-a.nc'; %seacat logger with temp and salinity data
    
   
%% find where to trim raw data
 %use start & stop time from stats processed file
time=ncread(sigfn,'time');
time2=ncread(sigfn,'time2') ; 

jt = time+time2/(3600*24*1000);
dn = j2dn(time,time2);

p1ac=ncread(sigfn,'P_1ac');
p1ac=squeeze(p1ac);

%clear dn jt time time2
%% Find time series of pressure sensor elevation above seabed
    
    %use brange from processed sig file and instrument offsets to get time series
    %of pressure sensor elevation above seabed
    
    brange=squeeze(ncread(sigfn,'brange'));%units are in m
        
    %find all of the offsets 
    zoff = 0;
    
    z_init = ncreadatt(sigfn,'P_1','initial_sensor_height');
    zp = brange; % elevation of pressure measurements [m] accounting for variable brange
    
    %zpmf=medfilt(zp,15); %smooth using median filter
%% need to convert corrected pressure to depth using temp & salinity data from deployed TC data

%surface T/C
mc(1).cond=squeeze(ncread(mcfn1,'C_51'));
mc(1).temp=squeeze(ncread(mcfn1,'T_28'));
mc(1).dn = nctime2dn(mcfn1)';
mc(1).depth=ncreadatt(mcfn1,'C_51','sensor_depth');

%nearbed T/C
mc(2).cond=squeeze(ncread(mcfn2,'C_51'));
mc(2).temp=squeeze(ncread(mcfn2,'T_28'));
mc(2).dn = nctime2dn(mcfn2)';
mc(2).depth=ncreadatt(mcfn2,'C_51','sensor_depth');

% NOTE- - Note near-bed Conductivity became fouled after just a couple
%of weeks so just use surface. there is little evidence of significant
%stratification in the overlapping record.

%find salinity from T&C record
c3515=sw_c3515*100/1000; %convert units to S/m
cond=mc(1).cond;
cndr=cond./c3515;
sal=sw_salt(cndr,mc(1).temp,mc(1).depth);


rhosw = sw_dens0(sal,mc(1).temp);
rhoswi=interp1(mc(1).dn,rhosw,dn);

depac = p1ac.*10^4./rhoswi./9.81; %depth of sensor
    
depth_corrected = zp+depac; % time series of depth [meters]
depth_corrected=depth_corrected';


save C:\Users\ssuttles\data\Matanzas\proc\mat\11102sig_depth_corrected depth_corrected p1ac zp z_init dn