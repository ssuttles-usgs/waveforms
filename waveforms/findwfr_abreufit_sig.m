function wfr=findwfr_abreufit_sig(wf)

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
t=wf(i).t*2*pi/T(i);
u=wf(i).ub./Uw(i);
[bb1]=abreu_fit(R(i),t',u');
b1(i)=bb1;
end

t=0:pi/64:2*pi;

[Uw_sort,I]=sort(Uw,'descend');

nsig=round(length(I)/3);
isig=I(1:nsig);

inp1=mean(R(isig));
inp2=mean(b1(isig));
iop1=3;%inp1=R and inp2=b1;
iop2=2;%velocity starts @ 0
iapp=1;%beta is approximate

wfr.Uw=mean(Uw(isig));
wfr.T=mean(T(isig));
wfr.t=t*mean(T(isig))./(2*pi);

[wfr.u,wfr.a,wfr.b,wfr.phi,wfr.R,wfr.beta,wfr.alpha]=abreu_fstream(t,inp1,inp2,wfr.Uw,wfr.T,iop1,iop2,iapp);
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


