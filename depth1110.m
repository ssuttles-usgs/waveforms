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

sig.jt = double(time)+double(time2)/(3600*24*1000);
sig.dn = j2dn(double(time),double(time2));

sig.p1ac=ncread(sigfn,'P_1ac');
sig.p1ac=squeeze(sig.p1ac);

clear jt time time2
%% Find time series of pressure sensor elevation above seabed
    
    %use brange from processed sig file and instrument offsets to get time series
    %of pressure sensor elevation above seabed
    
    sig.brange=squeeze(ncread(sigfn,'brange'));%units are in m
        
    %find all of the offsets 
    sig.zoff = 0;
    
    sig.pz_init = ncreadatt(sigfn,'P_1','initial_sensor_height');
    sig.zp = sig.brange; % elevation of pressure measurements [m] accounting for variable brange
    
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
sig.rhoswi=interp1(mc(1).dn,rhosw,sig.dn);

sig.depac = sig.p1ac.*10^4./sig.rhoswi./9.81; %depth of sensor
    
sig.depth_corrected = sig.zp+sig.depac; % time series of depth [meters]
%depth_corrected=depth_corrected';


clear time time2
save C:\Users\ssuttles\data\Matanzas\proc\mat\11102sig_depc sig

%% Find corrected depth for 11109vec
vecbfn = fullfile('C:\Users\ssuttles\data\Matanzas\proc\released\11109vecb-cal.nc'); % burstfile name
vecsfn = fullfile('C:\Users\ssuttles\data\Matanzas\proc\released\11109vecs-a.nc'); % statistics filename

time=ncread(vecsfn,'time');
time2=ncread(vecsfn,'time2');

vec.jt = double(time)+double(time2)/(3600*24*1000);
vec.dn = j2dn(double(time),double(time2));

vec.p1=squeeze(ncread(vecsfn,'P_4023'));

brangei=interp1(sig.dn,sig.brange,vec.dn); %intrepolated brange from Sig1K
p1aci=interp1(sig.dn,sig.p1ac,vec.dn);%use interpolated pressure from Sig1K

  %find all of the offsets  
    vec.zoff = ncreadatt(vecsfn,'/','ADVsampleVolumeOffset')*10; % need to multiply by 10 because metadata says 0.0149 sand should be appprox 15cm or 0.15 m
    vec.z_init = ncreadatt(vecsfn,'u_1205','initial_instrument_height');
    vec.pz_init = ncreadatt(vecsfn,'P_4023','initial_instrument_height');
    vec.pz_off=vec.pz_init-sig.pz_init; %offset between Sig1K and vector pressure
    vec.p1ac=p1aci-vec.pz_off;    
    vec.pdelz = vec.pz_init-(vec.z_init-vec.zoff); % elevation diff. between velocity and pressure- NEED to include zoff to get distance between velocity and pressure measuremens
    vec.zr = brangei-(vec.pdelz-vec.pz_off); % time series of Vector sample locations
    vec.zp = vec.zr+vec.pdelz; % elevation of pressure measurements [m] accounting for variable brange
    vec.rhoswi=interp1(mc(1).dn,rhosw,vec.dn);
    vec.depac = vec.p1ac.*10^4./vec.rhoswi./9.81; %depth of sensor
    vec.depth_corrected=vec.zp+vec.depac;
    
clear time time2 brangei p1aci
 save C:\Users\ssuttles\data\Matanzas\proc\mat\11109vec_depc vec
 
 %% Find corrected depth for SeaGauge 1110
  sgwfn='C:\Users\ssuttles\data\Matanzas\proc\released\111010sgw-cal.nc';
 sgrfn='C:\Users\ssuttles\data\Matanzas\proc\released\111010sgr-cal.nc';
 
time=ncread(sgwfn,'time');
time2=ncread(sgwfn,'time2');

sg.jt = double(time)+double(time2)/(3600*24*1000);
sg.dn = j2dn(double(time),double(time2));

p1ac=squeeze(ncread(sgrfn,'P_1ac'));

sg.p1ac=mean(p1ac)';

brangei=interp1(sig.dn,sig.brange,sg.dn); %intrepolated brange from Sig1K


  %find all of the offsets  
    sg.zoff = 0;
    sg.pz_init = ncreadatt(sgrfn,'P_1ac','initial_instrument_height');
    sg.pz_off=sg.pz_init-sig.pz_init; %offset between Sig1K and vector pressure
    sg.zp = brangei+sg.pz_off; % elevation of pressure measurements [m] accounting for variable brange
    sg.rhoswi=interp1(mc(1).dn,rhosw,sg.dn);
    sg.depac = sg.p1ac.*10^4./sg.rhoswi./9.81; %depth of sensor
    sg.depth_corrected=sg.zp+sg.depac;
    
clear time time2 brangei
 save C:\Users\ssuttles\data\Matanzas\proc\mat\111010sg_depc sg
 
 %% 
 
 whvfn='C:\Users\ssuttles\data\Matanzas\proc\released\11101whVp-cal.nc';
 whafn='C:\Users\ssuttles\data\Matanzas\proc\released\11101wh-a.nc';

% 
time=ncread(whvfn,'time');
time2=ncread(whvfn,'time2');

whv.jt = double(time)+double(time2)/(3600*24*1000);
whv.dn = j2dn(double(time),double(time2));

brangei=interp1(sig.dn,sig.brange,whv.dn); %intrepolated brange from Sig1K
p1aci=interp1(sig.dn,sig.p1ac,whv.dn);%use interpolated corrected pressure from Sig1K

  %find all of the offsets  
    whv.pz_init = ncreadatt(whafn,'P_1','initial_sensor_height');
    whv.pz_off=whv.pz_init-sig.pz_init; %offset between Sig1K and vector pressure
    whv.p1ac=p1aci-whv.pz_off;    
    whv.zp = brangei+whv.pz_off; % elevation of pressure measurements [m] accounting for variable brange
    whv.rhoswi=interp1(mc(1).dn,rhosw,whv.dn);
    whv.depac = whv.p1ac.*10^4./whv.rhoswi./9.81; %depth of sensor
    whv.depth_corrected=whv.zp+whv.depac;
    
clear time time2 brangei p1aci
 save C:\Users\ssuttles\data\Matanzas\proc\mat\11101whv_depc whv

 


