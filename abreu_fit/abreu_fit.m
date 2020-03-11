function [b1]=abreu_fit(R,t,u)
%
% [b1]=abreu_fit(R,t,u)
%
% Function to find b1 for free-stream velocity description of Abreu et al. (2010) 
% given t, u and R, and u is non-dimensional (umax-umin = 2). Depending on the length
% of t and u, the function either does an exact fit to tmax or tmin if the length is
% <= 2 (it is assumed that tmax < tmin) or a least-squares-fit if the length > 2, as 
% described in Malarkey and Davies (2012). This function calls abreu_fstream.m.
%
% INPUTS
% R  - umax/(umax-umin)                [value]
% t  - phase (0 <= t < 2*pi)           [column vector]
% u  - non-dimensional velocity at t   [column vector]
%
% OUTPUTS
% b1    b1=b*cos(phi)/abs(cos(phi))    [vector when length of t = 2 otherwise value]
%
% where b=r/(1+P), P=sqrt(1-r^2) and phi has the opposite sign to Abreu et al. (2010).
%
% References
% Abreu, T., P.A. Silva, F. Sancho and A. Temperville, 2010. Analytical approximate
% wave form for asymmetric waves. Coastal Engineering, 57, 656-667.
% Malarkey, J. and Davies, A.G. 2012. Free-stream velocity descriptions under waves 
% with skewness and asymmetry. Coastal Engineering, 68, 78-95.
%
% Jonathan Malarkey 2012

Nt=length(t);
if Nt>2
  b1=lsf(R,t,u);
else
  b1=[];
  for i=1:Nt
    ifit=1; if u(i)<0; ifit=2; end
    b1=[b1; findb1(R,t(i),ifit)];         
  end 
end

function [b1]=lsf(R,t,u)
Nt=length(t); c=2*R-1; Nb2=101; Nb=2*Nb2; bmax=0.99;
b=[-bmax+(bmax-abs(c))*[0:Nb2-1]/(Nb2-1) abs(c)+(bmax-abs(c))*[0:Nb2-1]/(Nb2-1)];
if c==0; b=-bmax+2*bmax*[0:Nb-1]/(Nb-1); end
uu=u*ones(1,Nb);
uuf=zeros(Nt,Nb);
for i=1:Nb
  [uuf(:,i),ad,bd,phid,Rd,betad,alphad]=abreu_fstream(t,R,b(i),1,2*pi,3,2,1);
end
rms=sqrt(mean((uu-uuf).^2,1)); [rmsmin,imin]=min(rms);
b1=b(imin);

function [b1]=findb1(R,t1,ifit)
zer=2e-16;
c=2*R-1;    t12=t1+asin(c); B1=(cos(t12))^2;
B2=2-B1+2*c*(c-2*sin(t12)); B3=B1+4*c^2*(c-sin(t12))^2;
ss=sign(pi/2-t12); if ifit==2; ss=sign(t12-3*pi/2); end
if B1<=zer; ss=1; b1=abs(c); else; b1=ss*sqrt((B2-sqrt(B2^2-B1*B3))./B1); end
