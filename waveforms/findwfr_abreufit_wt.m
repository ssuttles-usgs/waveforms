function wfr=findwfr_abreufit(wf)

% for i=1:length(wf)
% ubh(i)=wf(i).umax-wf(i).umin;
% T(i)=wf(i).T;
% end
% Uw=mean(ubh)/2;
% T=mean(T);

for i=1:length(wf)
Uw(i)=(wf(i).umax-wf(i).umin)/2;

% Use Eqn 8 from van der A for representative near-bed orbital amplitude
try
uhat(i)=sqrt((2/wf(i).T).*trapz(wf(i).t,wf(i).ub.^2));
catch
uhat(i)=NaN;
end

T(i)=wf(i).T;
R(i)=wf(i).R;
%Su(i)=wf(i).Su;
t=wf(i).t*2*pi/T(i);
u=wf(i).ub./Uw(i);
[bb1]=abreu_fit(R(i),t',u');
b1(i)=bb1;
end

t=0:pi/64:2*pi;


R_wt=sum(R.*Uw.^2.*T)./sum(Uw.^2.*T);
b1_wt=sum(b1.*Uw.^2.*T)./sum(Uw.^2.*T)
Uw_wt=sum(Uw.*Uw.^2.*T)./sum(Uw.^2.*T);
T_wt=sum(T.*Uw.^2.*T)./sum(Uw.^2.*T);
%Su_wt=sum(Su.*Uw.^2.*T)./sum(Uw.^2.*T);

inp1=mean(R_wt);
inp2=mean(b1_wt);
iop1=3;%inp1=R and inp2=b1;
iop2=2;%velocity starts @ 0
iapp=1;%beta is approximate

wfr.dn=wf(1).dn(1);
wfr.Uw=Uw_wt;
wfr.T=T_wt;
wfr.t=t*T_wt./(2*pi);

[wfr.u,wfr.a,wfr.b,wfr.phi,wfr.R,wfr.beta,wfr.alpha]=abreu_fstream(t,inp1,inp2,wfr.Uw,wfr.T,iop1,iop2,iapp);
%wfr.Su=Su_wt;
tz=findzs(wfr.u,wfr.t); %find zero-crossing
tz=tz(tz>0 & tz<wfr.T); %make sure it is crest to trough crossing
tz=tz(1)
wfr.umax=max(wfr.u);
wfr.umin=min(wfr.u);
wfr.Tc=tz;
wfr.Tt=wfr.T-tz;
wfr.Tcu=wfr.t(wfr.u==wfr.umax);
wfr.Ttu=wfr.t(wfr.u==wfr.umin)-wfr.Tc;
% wfr.Ac=trapz([wfr.t(wfr.t<tz) tz],[wfr.u(wfr.t<tz) 0]);
% wfr.At=trapz([tz wfr.t(wfr.t>tz)],[0 wfr.u(wfr.t>tz)]);


