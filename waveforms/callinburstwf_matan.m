clear
%load raw puv data from near-bed adv

projstr='Matanzas 2018'
%set burst number to process
bn=56
load C:\Users\ssuttles\data\Matanzas\proc\mat\puv_proc_Matan_iwaves raw depth fs zr zp dn

%% % LOAD the observational data from workhorse from Fire Island 
fn_wh='C:\Users\ssuttles\data\Matanzas\proc\released\11101whVp-cal.nc'
dn_Hs=nctime2dn(fn_wh);
Hs(:)=squeeze(ncread(fn_wh,'wh_4061'));
vspec=squeeze(ncread(fn_wh,'vspec'));
freq=squeeze(ncread(fn_wh,'frequency'));
%Td(:)=squeeze(wp_peak(1,1,:));
instht=ncreadatt(fn_wh,'hght_18','initial_sensor_height')
%h(:)=squeeze(hght_18)+instht; % extract depth; 

%import corrected depth from ADV 9917 extrenal pressure sensor
%depc=load('..\mat\9917advs_depth_corrected.mat');
h(:)=depth(:,1);

% load vspec data of uhat and Tr 
load('C:\Users\ssuttles\data\Matanzas\proc\mat\workhorse_matanzas_inlet.mat','uhat_emp','Tbr_combine')

u=raw(bn).u;
v=raw(bn).v;


 %u check for and interpret NaN
      iu=find(~isnan(u));
      ui=interp1(iu,u(iu),[1:length(u)],'linear','extrap')';
      
      %v check for and interpret NaN
      iv=find(~isnan(v));
      vi=interp1(iv,v(iv),[1:length(v)],'linear','extrap')';
 
u=ui;
v=vi;

%make input data to m-fcn
datin.bn=bn;
datin.dn=dn(bn);
datin.depth=depth(bn);
datin.fs=fs;
datin.zr=zr(bn)
datin.zp=zp(bn)
datin.u=u;
datin.v=v;
datin.p=raw(bn).p;
datin.t=[0:1/datin.fs:(length(datin.u)-1)./datin.fs]';


whin.dn=dn_Hs(bn);
whin.Hs=Hs(bn);
whin.vspec=vspec(:,bn);
whin.freq=freq;
whin.instht=instht;
whin.h=h(bn);
whin.uhat=uhat_emp(bn);
whin.Tr=Tbr_combine(bn);


%set filttype
filttype=2; %remove ig

nwf=0; %no wf sub-sampling
%nwf=20; %sub-sample wf burst

mkinburstwf(projstr,bn,datin,whin,filttype,nwf)

save2pdf(sprintf('pdfs\\burstwaveformview_matan_b%d',bn),gcf,300)