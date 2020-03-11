%% abreu_points
function [DTc, DTt, DTcu, DTtu, umax, umin, RR, beta] = ...
    abreu_points(r, phi, Uw, T );
%
%  Calculate umax, umin, and phases of asymmetrical wave orbital velocity
%
%  Use the asymmetry parameters from Ruessink et al, 2012
%  to get the umax, umin and phases of asymettrical wave
%  orbital velocity to be used by Van Der A.
%  T_c is duration of crest
%  T_cu Duration of accerating flow within crest half cycle
%
omega=2.0*pi/T;
%
phi_new=-phi;

% Malarkey and Davies (Under equation 16b)
P=sqrt(1.0-r*r);
%
% Malarkey and Davies (Under equation 16b)
%
b=r/(1.0+P);
%
% Appendix E of Malarkey and Davies
%
c=b*sin(phi_new);
%
cff1=4.0*c*(b*b-c*c)+(1.0-b*b)*(1.0+b*b-2.0*c*c) ;
cff2=(1.0+b*b)^2.0-4.0*c*c ;
ratio=cff1/cff2; 
if(ratio>1.0) 
  ratio=1.0  ;
end
if(ratio<-1.0) 
  ratio=-1.0 ;
end
tmc=asin(ratio) ;

%
cff1=4.0*c*(b*b-c*c)-(1.0-b*b)*(1.0+b*b-2.0*c*c);
cff2=(1.0+b*b)^2.0-4.0*c*c;
ratio=cff1/cff2; 
if(ratio>1.0) 
  ratio=1.0  ;
end
if(ratio<-1.0) 
  ratio=-1.0 ;
end
tmt=asin(ratio) ;


if(tmt<0.0)
    tmt=tmt+2.0*pi;
end
if(tmc<0.0)
    tmc=tmc+2.0*pi;
end
%
% Non dimensional umax and umin, under E5 in Malarkey and Davies
%
umax=1.0+c;
umin=umax-2.0;
%
%       Dimensionalize
%
umax=umax*Uw;
umin=umin*Uw;
%
% phase of zero upcrossing and downcrossing (radians)
%
 tzu=asin(b*sin(phi_new)) ;
tzd=2.0*acos(c)+tzu;
%
% MD, equation 17
%
RR=0.5*(1.0+b*sin(phi_new));
%
% MD, under equation 18
%
if(r<=0.5)
    F0=1.0-0.27*(2.0*r)^(2.1);
else
    F0=0.59+0.14*(2.0*r)^(-6.2);
end
%
% MD, Equation 15a,b
%
if(r >= 0.0 && r<0.5)
    betar_0=0.5*(1.0+r);
elseif(r>0.5 && r<1.0)
    cff1=4.0*r*(1.0+r);
    cff2=cff1+1.0;
    betar_0=cff1/cff2;
end
    betar_0=0.5 ; 

disp('Warning that beta is hardwired to be 0.5')
%
% MD, Equation 18
%
cff=sin((0.5*pi-abs(phi_new))*F0)/sin(0.5*pi*F0);
beta=0.5+(betar_0-0.5)*cff;
%
% MD, Table 1, get asymmetry parameterization
% using GSSO (10a,b)
%
cff=sqrt(2.0*(1.0+b*b)^3.0);
Sk=3.0*sin(phi_new)/cff;
Ak=-3.0*cos(phi_new)/cff;
%
% These are the dimensional fractions of wave periods needed by Van der A eqn.
% TSK - Check source of these equations
%
w=1.0/omega;
DTc=(tzd-tzu)*w;
DTt=T-DTc;
DTcu=(tmc-tzu)*w;
DTtu=(tmt-tzd)*w;
%
T_tu=tzd*w;
T_cu=tzu*w;
T_c=tmc*w;
T_t=tmt*w;
%
%fprintf(fid,'R, beta %f, %f\n', RR, beta);
return
end % function abreu_points
