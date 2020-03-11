function [u,a,b,phi,R,beta,alpha]=abreu_fstream(t,inp1,inp2,U1s,T,iop1,iop2,iapp)
%
% [u,a,b,phi,R,beta,alpha]=abreu_fstream(t,inp1,inp2,U1s,T,iop1,iop2,iapp)
%
% Function to calculate the dimensional or non-dimensional free-stream velocity, u,
% and acceration, a at phases specified by t of Abreu et al. (2010).  The function 
% allows the skewness and asymmetry to be specified through (r,phi), (b,phi) or
% (R,b1). The function also returns b, phi, R, beta and alpha. beta is either
% calculated approximately if iapp=1 or exactly if iapp=0. All other measures of
% skewness and asymmetry can be readily calculated from b and phi (see Table 1, 
% Malarkey and Davies, 2012).
%
% INPUTS
% t     - phase (0 <= t < 2*pi)                                        [vector]
% inp1  - r (-), b (-)  or R (-)  (see iop1)                           [value]
% inp2  - phi (rad) or b1 (-)   (see iop1)                             [value]
% U1s   - dimensional velocity amplitude (for non-dim. U1s=1)          [value]
% T     - dimensional wave period (for non-dim. T=2*pi)                [value]
% iop1  - iop1=1, for r=inp1, phi=inp2; iop1=2, for b=inp1, phi=inp2;  [integer]
%         iop1=3, for R=inp1, b1 =inp2
% iop2  - iop2=2 velocity starts from zero up-crossing                 [integer]
% iapp  - iapp=1 beta is approximated else beta is calculated exactly  [integer]
%
% OUTPUTS
% u     - velocity at t                                                [vector]
% a     - acceleration at t                                            [vector]
% b     - skewness/asymmetry parameter r=2b/(1+b^2)                    [value]
% phi   - skewness/asymmetry parameter                                 [value]
% R     - umax/(umax-umin)                                             [value]
% beta  - amax/(amax-amin)                                             [value]
% alpha - tmax/pi                                                      [value]
%
% Input parameter ranges
% 0 <= r < 1; -pi < phi <= pi; 0 <= b < 1; 0 <= R < 1; -1 < b1 < 1.
% phi is defined with the opposite sign to Abreu et al. (2010), see Malarkey and
% Davies (2012).  When R and b1 are used as inputs, then function forces 
% b1=sign(b1)*|c| when |b1| < |c|, where c=2R-1.
%
% References
% Abreu, T., P.A. Silva, F. Sancho and A. Temperville, 2010. Analytical approximate
% wave form for asymmetric waves. Coastal Engineering, 57, 656-667.
% Malarkey, J. and Davies, A.G. 2012. Free-stream velocity descriptions under waves 
% with skewness and asymmetry. Coastal Engineering, 68, 78-95.
%
% Jonathan Malarkey 2012

global b r phi c

% Calculation of initial parameters from inputs
A1s=2*pi*U1s/T;
if iop1<=2
  if iop1==1; r=inp1; b=r/(1+sqrt(1-r^2)); end
  if iop1==2; b=inp1; r=2*b/(1+b^2); end  
  phi=inp2;   c=b*sin(phi); R=(1+c)/2;
  b1=b; if abs(phi)>pi/2; b1=-b; end
else
  R=inp1; b1=inp2; c=2*R-1;
  if abs(b1)<abs(c)
    ['warning |b1| must be less than or equal to |c| = ',num2str(abs(c)),...
     ', b1 is being set to sign(b1)|c|']
    b1=sign(b1)*abs(c)
  end
  b=abs(b1); r=2*b/(1+b^2);
  phi=0; if b~=0; phi=atan2(c,b1*sqrt(1-c^2/b^2)); end
end

u=U1s*sin(t); a=A1s*cos(t); beta=0.5; alpha=0.5;
if b~=0
% Calculation of u and a
  gam=asin(c); t1=t; if iop2==2; t1=t+gam; end; t2=t1-phi;
  u=U1s*(1-b^2).*(sin(t1)-c)./(1+b^2-2*b*cos(t2));
  a=A1s*calca(t1);
% Calculation of alpha
  alpha1=pi/2;
  if abs(phi)~=pi/2
    alpha1=asin((4*c*(b^2-c^2)+(1-b^2)*(1+b^2-2*c^2))/((1+b^2)^2-4*c^2));
  end
  alpha=(alpha1-gam)/pi; if abs(phi)>pi/2; alpha=(pi-alpha1-gam)/pi; end
% Calculation of beta
  [beta]=calcbeta(iapp);
end

function [a]=calca(t1)
global b r phi c
t2=t1-phi;
a=(1-b^2)*((1+b^2)*cos(t1)-2*b*cos(phi)+2*b*c*sin(t2))./(1+b^2-2*b*cos(t2)).^2;

function [beta]=calcbeta(iapp)
global b r phi c
if iapp==1 | phi==0 | abs(phi)==pi | abs(phi)==pi/2
  % Calculate beta approximately (approximation is exact if phi = 0, +/-pi/2, +/-pi)
  beta0=(1+r)/2; F0=1-0.27*(2*r)^2.1;
  if r>0.5; beta0=4*r.*(1+r)./(4*r.*(1+r)+1); F0=0.59+0.14*(2*r)^(-6.2); end
  beta=0.5+(beta0-0.5)*sin(F0*(pi/2-abs(phi)))/sin(0.5*pi*F0);
else
  % Calculate beta exactly
  P=sqrt(1-r^2); sp=sin(phi); cp=cos(phi);
  B4=r^2*(1-(r*sp)^2);
  B3=2*r*cp*(1-2*r^2+r^2*(1-2*P)*sp^2);
  B2=6*P*(r*sp)^2*(1-(1-P)*sp^2)+4*P^4-3*P^2;
  B1=-2*r*cp*(1-2*(r*cp)^2+(2-3*r^2-2*P^3)*sp^4);
  B0=(1-P)*sp^2*(-1-3*P+4*r^2+4*P*r^2-(1-P+2*r^2)*sp^2+r^2*(1-P)*sp^4)+(2*r*P)^2-1;
  xx=roots([B4 B3 B2 B1 B0]); ii=find(imag(xx)==0 & abs(real(xx))<=1);
  tm=acos(xx(ii)); if phi<0; tm=2*pi-acos(xx(ii)); end
  am=calca(tm); beta=max(am)/(max(am)-min(am));
end
