function [UBS,PUV,wf,wh,wfr,wfrwt]=mkinburstwf(projstr,bn,datin,whin,filttype,nwf)


%make input data to m-fcn
fs=datin.fs;
[wf,UBS,PUV]=burstwfcalcs(datin,filttype);

dnsb=datestr(datin.dn);
if 1
     figure; clf
     plot(UBS.ur,UBS.vr,'.')
     title(sprintf('Rotated Burst Velocities burst = %d , %s',bn,dnsb))
     xlabel('ur, [m ^. s^{-1}]')
     ylabel('vr, [m ^. s^{-1}]')
     set(gca,'xlim',[-0.5 0.5]);
     set(gca,'ylim',[-0.5 0.5])
     pause(0.1)
end

UBSnofilt = ubstatsr(detrend(datin.u), detrend(datin.v), fs );

%% % LOAD the observational data from workhorse from Fire Island 

%find waveform from surface wave data
[r, phi, Ur] = skewness_params(whin.Hs, whin.Tr, whin.h);
% P=sqrt(1-r^2);
% b=r/(1+P);
% b1=b*cos(phi)/abs(cos(phi));

[T_c, T_t, T_cu, T_tu, umax, umin, RR, beta] = ...
    abreu_points(r, phi, whin.uhat, whin.Tr); %changes sign of phi in m-fcn

 
t=0:pi/64:2*pi;
wh.T=whin.Tr;
wh.t=t*wh.T/(2*pi);
[wh.u,wh.a,wh.b,wh.phi,wh.R,wh.beta,wh.alpha]=abreu_fstream(t,r,-phi,whin.uhat,whin.Tr,1,2,1); %change sign of phi to correct known issus with skewness_params m-fcn
%[u,a,b,phi,R,beta,alpha]=abreu_fstream(t,inp1,inp2,U1s,T,iop1,iop2,iapp)
tz=findzs(wh.u,wh.t); %find zero-crossing
tz=tz(tz>0 & tz<wh.T); %make sure it is crest to trough crossing
tz=tz(1)
wh.umax=max(wh.u);
wh.umin=min(wh.u);
wh.Tc=tz;
wh.Tt=wh.T-tz;
wh.Tcu=wh.t(wh.u==wh.umax);
wh.Ttu=wh.t(wh.u==wh.umin)-wh.Tc;


%%


tall=([wf.dn]-[wf(1).dn(1)])*86400;
uball=[wf.ub];

len=length(wf);

if nwf==0
    istart=2
else
    istart=1
end

for jj=istart:2
    
    if jj==1
        figure(1)
        clf
        scrsz=get(0,'screensize');
        %set(gcf,'position',[scrsz(3)+50 -50 scrsz(3)*0.8 scrsz(4)*0.8]);
        set(gcf,'position',[50 50 scrsz(3)*0.8 scrsz(4)*0.8]);
        ax0=axes('position',[0.1 0.725 0.8 0.2]);
        ax1=axes('position',[0.1 0.475 0.8 0.2]);
        ax2=axes('position',[0.1 0.1 0.2 0.3])
        ax3=axes('position',[0.4 0.1 0.2 0.3])
        ax4=axes('position',[0.7 0.1 0.2 0.3])

    else
        nwf=len;
        figure(2)
        clf
        scrsz=get(0,'screensize');
        %set(gcf,'position',[scrsz(3)+50 -50 scrsz(3)*0.8 scrsz(4)*0.8]);
        set(gcf,'position',[50 50 scrsz(3)*0.8 scrsz(4)*0.8]);
        ax0=axes('position',[0.1 0.725 0.8 0.2]);
        ax1=axes('position',[0.1 0.475 0.8 0.2]);
        ax2=axes('position',[0.1 0.1 0.2 0.3])
        ax3=axes('position',[0.4 0.1 0.2 0.3])
        ax4=axes('position',[0.7 0.1 0.2 0.3])
    end



for ii=1:len-nwf+1

if exist('ax0')
    cla(ax0,'reset')
end
    
if exist('ax1')
    cla(ax1,'reset')
end
if exist('ax2')
    cla(ax2,'reset')
end
if exist('ax3')
    cla(ax3,'reset')
end

if exist('ax4')
    cla(ax4,'reset')
end

axes(ax0)
plot(datin.t,datin.u)
hold on
plot(datin.t,datin.v,'r')
set(ax0,'ylim',[-1.5 1.5]);
ylabel('[m ^. s^{-1}]')
legend({'u_{direct}' 'v_{direct}'},'Location','NorthEast','box','off')
text(ax0,25,1.2,'a)')

axes(ax1)    
plot(tall,uball,'b')
hold on
set(ax1,'ylim',[-1.5 1.5]);
ylabel('[m ^. s^{-1}]')
xlabel('burst time, [s]')
text(ax1,25,1.2,'b)')

swf=ii;
cnt=0;
for i=swf:swf+nwf-1
    cnt=cnt+1;
axes(ax2)
Uw(cnt)=(wf(i).umax-wf(i).umin)./2;
uhat(cnt)=sqrt((2/wf(i).T).*trapz(wf(i).t,wf(i).ub.^2));
plot(wf(i).t./wf(i).T,wf(i).ub./Uw(cnt))
hold on
end

set(ax2,'box','on')
line([0 1],[0 0],'color',[0.75 0.75 0.75])
set(ax2,'ylim',[-2 2]);
xlabel(ax2,'t/T')
ylabel(ax2,'u/U_w')
text(0.02,1.8,'c)')

wfr=findwfr_abreufit(wf(swf:swf+nwf-1));

%find significant waveforms using highest 1/3 of wave orbital amps
wfi=wf(swf:swf+nwf-1);

[Uw_sort,I]=sort(Uw,'descend');
nsig=round(length(I)/3);
isig=I(1:nsig);

wfrs=findwfr_abreufit(wfi(isig));
%cla(ax2)
for jj=1:length(isig)
plot(wfi(isig(jj)).t./wfi(isig(jj)).T,wfi(isig(jj)).ub./Uw(isig(jj)),'--k')
end

line([0 1],[0 0],'color',[0.75 0.75 0.75])
set(ax2,'ylim',[-2 2]);
xlabel(ax2,'t/T')
ylabel(ax2,'u/U_w')
text(0.02,1.8,'c)')
set(ax2,'box','on')


wfrwt=findwfr_abreufit_wt(wf(swf:swf+nwf-1));



% wfr=findwfr(wf(swf:swf+nwf-1));
% tic
% wfr.wf7=findwf7(wfr);
% toc

axes(ax3)
line([0 14],[0 0],'color',[0.75 0.75 0.75])
set(ax3,'xlim',[0 14]);
set(ax3,'ylim',[-1 1]);
box on
hold on
wfrpts=[0 0; wfr.Tcu wfr.umax; wfr.Tc 0; wfr.Tc+wfr.Ttu wfr.umin; wfr.T 0];
plot(wfrpts(:,1),wfrpts(:,2),'ok','MarkerFaceColor','r','MarkerSize',4)
%plot(wfr.wf7.t,wfr.wf7.ub,'-k')
plot(wfr.t,wfr.u,'-k');
plot(wh.t,wh.u,'r');
plot(wfrs.t,wfrs.u,'--k');
plot(wfrwt.t,wfrwt.u,'m')
legend({'0' 'wfr pts' 'wfr_{direct}' 'wh_{param}' 'wfr_{direct_{sig}}' 'wfr_{direct_{wt}}'},'box','off')
xlabel('t,[s]')
ylabel('u, [m ^. s^{-1}]')
text(ax3,0.25,0.8,'d)')


axes(ax1)
iplt=swf:swf+nwf-1;
tplt=([wf(iplt).dn]-[wf(1).dn(1)])*86400;
ubplt=[wf(iplt).ub];

Su=mean(ubplt.^3)./std(ubplt).^3;


plot(tplt,ubplt,'r')

%find t & ub significant continuous waveforms
t0=0;
tsig=[];
ubsig=[];
for kk=1:length(isig);
    tsig=[tsig wfi(isig(kk)).t+t0];
    ubsig=[ubsig wfi(isig(kk)).ub];
    t0=tsig(end);
end
   

legend({'wf_{direct} all' 'wf_{direct} selected'},'box','off')



%plot power spectra of near-bed adv data segment and wh-vspec
axes(ax4)
%interpolate ub data to be continuous time-series
% dt=1/4; % 4 Hz
% fsi=1/dt;
% 
% if 1 %all waveforms in section
% ti=tplt(1):dt:tplt(end);
% [tt,ind]=unique(tplt);
% ui=interp1(tplt(ind),ubplt(ind),ti);
% end
% 
% if 0 %significant waveforms only
% ti=tsig(1):dt:tsig(end);
% [tt,ind]=unique(tsig);
% ui=interp1(tsig(ind),ubsig(ind),ti);
% end

% plot spectra- raw 
ui=UBSnofilt.ur;
nsegs=16;
nsamps=length(ui);
pow2=nextpow2(nsamps/nsegs);
nfft=2^pow2
nx = length(ui);
na = nx/nfft;
w = hanning(floor(nx/na));
[Uxx,F]=pwelch(ui,w,0,nfft,fs);
wns=wavenumber(2*pi*F,datin.depth*ones(size(F)));
%Kz=velx(wns,datin.zr,datin.depth);
cutoff=0.1;
z=datin.zr;
KZ=cosh(datin.zr*wns)./sinh(datin.depth*wns);
FLKZ=(2*pi*F).*KZ;
%FLKZ=KZ;
%find the cutoff and replace the values below it
cutindex=(FLKZ < cutoff);
FLKZ(cutindex)=cutoff;

%Cut off spectrum and apply linear transformation

Nxx=Uxx./FLKZ.^2;

plot(F(~cutindex),Nxx(~cutindex));

hold on

%spectra filtered
ui=UBS.ur;
nsegs=16;
nsamps=length(ui);
pow2=nextpow2(nsamps/nsegs);
nfft=2^pow2
nx = length(ui);
na = nx/nfft;
w = hanning(floor(nx/na));
[Uxx,F]=pwelch(ui,w,0,nfft,fs);
wns=wavenumber(2*pi*F,datin.depth*ones(size(F)));
%Kz=velx(wns,datin.zr,datin.depth);
cutoff=0.1;
z=datin.zr;
KZ=cosh(datin.zr*wns)./sinh(datin.depth*wns);
FLKZ=(2*pi*F).*KZ;
%FLKZ=KZ;
%find the cutoff and replace the values below it
cutindex=(FLKZ < cutoff);
FLKZ(cutindex)=cutoff;

%Cut off spectrum and apply linear transformation

Nxx=Uxx./FLKZ.^2;

plot(F(~cutindex),Nxx(~cutindex),'--k');


plot(whin.freq,(whin.vspec./1000).^2,'r')

set(ax4,'xscale','log')
set(ax4,'yscale','log')
set(ax4,'ylim',[1e-6 1e2])

xlabel('Frequency, [Hz]')
ylabel('Spectral Power, [m^2 ^. Hz^{-1}]')

legend({'direct-unfiltered' 'direct-bandpass' 'param'},'box','off','Location','SouthEast')

text(ax4,1e-2,20,'e)')
title(ax0,sprintf('%s, site 3 Burst Waveforms from direct and parameterized methods- burst %d , %s ',projstr,bn,datestr(datin.dn)),'fontsize',12)

box on

figure(gcf)
pause(1)

clear Uw uhat
end

end
 
%%
if 0
        figure(3)
        clf
        scrsz=get(0,'screensize');
        %set(gcf,'position',[scrsz(3)+50 -50 scrsz(3)*0.8 scrsz(4)*0.8]);
        set(gcf,'position',[100 100 scrsz(3)*0.7 scrsz(4)*0.5]);
        ax5=axes('position',[0.075 0.15 0.4 0.7]);
        ax6=axes('position',[0.575 0.15 0.4 0.7]);
        axttl=axes('position',[0.1 0.9 0.9 0.05],'visible','off','color','none')
        
        axes(ax5)
        line([0 14],[0 0],'color',[0.75 0.75 0.75])
        set(ax5,'xlim',[0 14]);
        set(ax5,'ylim',[-1 1]);
        box on
        hold on
        wfrpts=[0 0; wfr.Tcu wfr.umax; wfr.Tc 0; wfr.Tc+wfr.Ttu wfr.umin; wfr.T 0];
        plot(wfrpts(:,1),wfrpts(:,2),'ok','MarkerFaceColor','r','MarkerSize',4)
        %plot(wfr.wf7.t,wfr.wf7.ub,'-k')
        plot(wfr.t,wfr.u,'-k');
        plot(wh.t,wh.u,'r');
        legend({'0' 'wfr pts' 'wfr_{direct}' 'wh_{param}'},'box','off')
        xlabel('t, [s]')
        ylabel('u, [m ^. s^{-1}]')
        text(0.2,0.8,'a)')
        
        axes(ax6)
        plot(F(~cutindex),Nxx(~cutindex));
        hold on
        plot(freq,(vspec(:,bn)./1000).^2,'r')

        set(ax6,'xscale','log')
        set(ax6,'yscale','log')
        set(ax6,'ylim',[1e-6 1e2])

        xlabel('Frequency, [Hz]')
        ylabel('Spectral Power, [m^2 ^. Hz^{-1}]')
        text(9e-3,15,'b)')
        legend({'wf_{direct}' 'wh_{vspec}'},'box','off')
        
        %text(axttl,0.5,0.5,sprintf('Fire Island 2014, site 3 Burst Waveforms from direct and parameterized methods- burst %d , %s ',...
         %   bn,datestr(datin.dn)),'fontsize',12,'horizontalalignment','center')
            
        saveas(gcf,sprintf('pngfiles\\waveform_compare_b%d_%s.png',bn,datestr(datin.dn,30)),'png')
        save(sprintf('mat\\waveform_compare_b%d_%s.mat',bn,datestr(datin.dn,30)),'wh','wfr','datin')
end
        
