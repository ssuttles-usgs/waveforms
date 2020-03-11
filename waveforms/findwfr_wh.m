function wfr=findwfr_wh(Hs,uhat,Tr,h)
%find waveform from surface wave data
%INPUT:
%
%Hs- significant wave height from surface [m]
%uhat - orbital velocity amplitutde [m/s]
%Tr - representative wave period [s]
%h - total water depth [m]

[r, phi, Ur] = skewness_params(Hs,Tr,h);
% P=sqrt(1-r^2);
% b=r/(1+P);
% b1=b*cos(phi)/abs(cos(phi));

[T_c, T_t, T_cu, T_tu, umax, umin, RR, beta] = ...
    abreu_points(r, phi, uhat, Tr); %changes sign of phi in m-fcn

 
t=0:pi/64:2*pi;
wfr.T=Tr;
wfr.t=t*wfr.T/(2*pi);
[wfr.u,wfr.a,wfr.b,wfr.phi,wfr.R,wfr.beta,wfr.alpha]=abreu_fstream(t,r,-phi,uhat,Tr,1,2,1); %change sign of phi to correct known issus with skewness_params m-fcn
%[u,a,b,phi,R,beta,alpha]=abreu_fstream(t,inp1,inp2,U1s,T,iop1,iop2,iapp)
wfr.Uw=(max(wfr.u)-min(wfr.u))/2
wfr.umax=max(wfr.u);
wfr.umin=min(wfr.u);
wfr.Ur=Ur;

try
    tz=findzs(wfr.u,wfr.t); %find zero-crossing
    tz=tz(tz>0 & tz<wfr.T); %make sure it is crest to trough crossing
    tz=tz(1)
    wfr.Tc=tz;
    wfr.Tt=wfr.T-tz;
    wfr.Tcu=wfr.t(wfr.u==wfr.umax);
    wfr.Ttu=wfr.t(wfr.u==wfr.umin)-wfr.Tc;

catch
    
    wfr.Tc=NaN;
    wfr.Tt=NaN;
    wfr.Tcu=NaN;
    wfr.Ttu=NaN;
end
    