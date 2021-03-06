%skewness FI14 plots

%load data

%uplookng ADCP-WH parameterized waveformsfrom Ruessink etal and Abreu etal eqns.
wh=load('C:\Users\ssuttles\data\FireIsland\analysis\site3\mat\9921wh_wfr_Su.mat'); 

%near-bed adv bandpass filtered and rotated into along wave direction
advbp=load('C:\Users\ssuttles\data\FireIsland\analysis\site3\mat\puv_proc_FI_iwaves2.mat')

%uplooking wh-V 5-beam ADCP processed wave data
f.vr='C:\Users\ssuttles\data\FireIsland\analysis\site3\data\9922whr-cal404-trm.cdf'

whv.dn=nctime2dn(f.vr);
whv.p1=squeeze(ncread(f.vr,'press'));
whv.vel=squeeze(ncread(f.vr,'vel'));
whv.strk=squeeze(ncread(f.vr,'strk'));


%plot 1- burst avg'd skewness whole deployment param v iwaves
figure(1)
set(gcf,'position',[100 100 1100 600])
sk=[[advbp.UBS.ur_sk]' [wh.wfr.u_sk]'];
hp=plot(advbp.dn,sk);
set(hp(2),'color','k');
set(hp(2),'linewidth',1.5);
datetick('x')
xlim=get(gca,'xlim');
line([xlim(1) xlim(2)],[0 0],'color',[0.5 0.5 0.5]);

legend({'adv: ur' 'wh: wfr_{param}' '0'},'box','off','fontsize',12)
ylabel('S_u','fontsize',14)
xlabel('2014');
set(gca,'fontsize',14)


saveas(gcf,'..\pngfiles\plt1_skew.png','png')

%add pressure skewness- find from raw signal then filter
nt1=1; nt2=2044 %all burst
T_long= 20; %longest wave period in bandpass filter (s)
T_short = 4; %shortest wave period in bandpass filter (s)
for ii=nt1:nt2
    p=advbp.raw(ii).p;
    p=detrend(p);
    p=iwavesbp(p,advbp.fs,T_long,T_short);
    p_sk(ii)=skewness(p);
end

figure(2)
clf
set(gcf,'position',[100 100 1100 600])
sk=[[advbp.UBS.ur_sk]' p_sk' [wh.wfr.u_sk]'];
hp=plot(advbp.dn,sk);
set(hp(3),'color','k');
set(hp(3),'linewidth',1.5);
set(hp(2),'color',[0.1 0.5 0.1])
set(hp(2),'linestyle','-')
datetick('x')
xlim=get(gca,'xlim');
line([xlim(1) xlim(2)],[0 0],'color',[0.5 0.5 0.5]);

legend({'adv: ur' 'adv: p' 'wh: wfr_{param}'  '0'},'box','off','fontsize',12)
ylabel('S_u','fontsize',14)
xlabel('2014');
set(gca,'fontsize',14)

saveas(gcf,'..\pngfiles\plt1_skew_p.png','png')

if 0 %turn on-off surface track plot
    %add skewness of surface tracking from WH-V at site 3
    whv.dn=nctime2dn(f.vr);
    whv.p1=squeeze(ncread(f.vr,'press'));
    whv.vel=squeeze(ncread(f.vr,'vel'));
    whv.strk=squeeze(ncread(f.vr,'strk')); %units are millimeters
    whv.samplerate=ncreadatt(f.vr,'/','WavesMonCfg.SampleRate');

    strk=squeeze(whv.strk(1,:,:))./1000; %convert to meters

    for ii=1:length(strk)
        s(:,ii)=detrend(strk(:,ii));
        s(:,ii)=iwavesbp(s(:,ii),whv.samplerate,T_long,T_short);
        s_sk(ii)=skewness(s(:,ii));
    end

    %despike s_sk using mdeian filter threshold method
    s_ski=despike(s_sk,0.5,2);

     hold on
     hp4=plot(whv.dn(1:2044),s_ski(1:2044),'color','m')
     legend({'adv: ur' 'adv: p' 'wh: wfr_{param}'  '0' 'whv: strk'},'box','off','fontsize',12)
     set(gca,'ylim',[-0.5 0.5])
end


%% No filtering

%near-bed adv unfiltered and rotated into along wave direction
advnf=load('C:\Users\ssuttles\data\FireIsland\analysis\site3\mat\puv_proc_FI_nofilt.mat')

%plot 1- burst avg'd skewness iwave v unfiltered v wh_param
figure(3)
set(gcf,'position',[100 100 1100 600])
sk=[[advbp.UBS.ur_sk]' [advnf.UBS.ur_sk]' [wh.wfr.u_sk]'];
hp=plot(advbp.dn,sk);
set(hp(3),'color','k');
set(hp(3),'linewidth',1.5);
datetick('x')
xlim=get(gca,'xlim');
line([xlim(1) xlim(2)],[0 0],'color',[0.5 0.5 0.5]);

legend({'adv: ur_{bp}' 'adv: ur_{nofilt}' 'wh: wfr_{param}' '0'},'box','off','fontsize',12)
ylabel('S_u','fontsize',14)
xlabel('2014');
set(gca,'fontsize',14)


saveas(gcf,'..\pngfiles\plt1_skew_nofilt.png','png')


%% waveforms - all direct

load('C:\Users\ssuttles\data\FireIsland\analysis\site3\mat\waveforms\mat\9917adv_waveforms_burstwfcalcs_bp')
for ii=1:length(b)
dn=[b(ii).wf.dn];

tt=(dn-dn(1)).*86400;
[t,I]=unique(tt);
ub=[b(ii).wf.ub];
ub=ub(I);
ubi=interp1(t,ub,[t(1):1/8:t(end)]);
wf(ii).ub_sk=skewness(ubi);
%Su(ii)=mean(ubi.^3)./std(ubi).^3;
end

%plot 1- burst avg'd skewness whole deployment param v iwaves
figure(4)
clf
set(gcf,'position',[100 100 1100 600])
sk=[[advbp.UBS.ur_sk]' [wf.ub_sk]' [wh.wfr.u_sk]'];
hp=plot(advbp.dn,sk);
set(hp(3),'color','k');
set(hp(3),'linewidth',1.5);
datetick('x')
xlim=get(gca,'xlim');
line([xlim(1) xlim(2)],[0 0],'color',[0.5 0.5 0.5]);

legend({'adv: ur_{bp}' 'adv: wf_{alldirect-bp}' 'wh: wfr_{param}' '0'},'box','off','fontsize',12)
ylabel('S_u','fontsize',14)
xlabel('2014');
set(gca,'fontsize',14)


saveas(gcf,'..\pngfiles\plt1_skew_wfall.png','png')


