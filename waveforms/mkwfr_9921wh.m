clear;

%% % LOAD the observational data from workhorse from Fire Island 
fn_wh='C:\Users\ssuttles\data\FireIsland\analysis\Taran\9921whp-cal.nc'
dn_Hs=nctime2dn(fn_wh);
Hs(:)=squeeze(ncread(fn_wh,'wh_4061'));
vspec=squeeze(ncread(fn_wh,'vspec'));
freq=squeeze(ncread(fn_wh,'frequency'));
%Td(:)=squeeze(wp_peak(1,1,:));
instht=ncreadatt(fn_wh,'hght_18','initial_sensor_height')
%h(:)=squeeze(hght_18)+instht; % extract depth; 

%import corrected depth from ADV 9917 extrenal pressure sensor
depc=load('..\mat\9917advs_depth_corrected.mat');
h(:)=depc.depth_corrected;

% load vspec data of uhat and Tr 
load('..\mat\vspec_uhat_tr_depc.mat','uhat_wh','Tr_wh')


for ii=1:length(h)
    wfr(ii)=findwfr_wh(Hs(ii),uhat_wh(ii),Tr_wh(ii),h(ii));
end

save mat\9921wh_wfr wfr