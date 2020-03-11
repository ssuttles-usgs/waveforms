function [U1oUA,R,tmax,t0,tmin,U0,u,a]=elfrinknd(t,H,L,xi)
%
% [U1oUA,R,tmax,t0,tmin,U0,u,a]=elfrinknd(t,H,L,xi)
%
% Function to calculate the non-dimensional free stream velocity (u) and acceleration
% (a) at t shallow water based on the characterisation of Elfrink et al. (2006). The 
% velocity and acceleration are non-dimensionalsed by U1s, such that umax-umin=2. The
% procedure is as described in Malarkey and Davies (2012).
%
% INPUTS
% t     - phase 0 <= t <= 2*pi                                  [value]
% H     - Hs/hs wave height to water depth ratio                [value]
% L     - Ls/hs wavelength to water depth ratio                 [value]
% xi    - surf similarity parameter (=tan(theta)/sqrt(Hs/L0s))  [value]
%
% OUTPUTS
% R     - umax/(umax-umin)                                      [value]
% U1oUA - U1s/UAirys                                            [value] 
% U0    - Velocity correction                                   [value]
% t0    - down zero-crossing phase                              [value]
% tmax  - phase of umax                                         [value]
% tmin  - phase of umin                                         [value]
% u     - non-dimensional velocity                              [vector]
% a     - non-dimensional acceleration                          [vector]
%
% where Hs and Ls are the dimensional wave height and wave length tan(theta) is the 
% beach slope, L0s is the deep-water wavelength, U1s = (umax-umin)/2 and UAirys is 
% the velocity amplitude from linear theory.
%
% References
% Elfrink, B., D.M. Hanes and B.G. Ruessink. 2006. Parameterization and simulation of
% near bed orbital velocities under irregular waves in shallow water. Coastal 
% Engineering, 53, 915-927.
% Malarkey, J. and Davies, A.G. 2012. Free-stream velocity descriptions under waves 
% with skewness and asymmetry. Coastal Engineering, 68, 78-95.
%
% Jonathan Malarkey 2012                                 % Elfink et al. (2006)
                                                         % parameter/constant names
     
Ur=H*L^2;                                                % Ur (Ursell number)

% Calculation of U1oUA                                       
A1=.4758*1.16-0.0145; B1=1.2001*1.16;                    % 1.2001*b2; 0.4758*b2+a2
N1=2*xi*sqrt(Ur/L)+1;                                    % sqrt(d3*Ur/L)
P1=N1*((L^0.5-tanh(abs(3*xi+2*L/Ur)))/(Ur*N1^2+1))^0.5;  % d5
U1oUA=B1*P1+A1;                                          % U2 (= Us/2*UAiry)

% Calculation of R
A2=0.38989;  B2=0.5366;                                  % a1; b1;
N2=xi*(11-L)-abs(L-10-H+abs(xi));                        % c3-c2
P2=(1-sqrt(H)*sqrt(abs(xi)+tanh(tanh(abs(N2)/Ur))));     % P1/sqrt(H)
R=B2*sqrt(H)*P2+A2;                                      % U1 (= UC/(UC+UT))

% Calculation of tmax
A3=-0.0005; B3=-0.2615; C3=-9.3852;                      % a3; b3;
N3=tanh(9.8496^2*xi^3*L*H^3)+tanh(H*L*xi)+L-1;           % e3
tmax=2*pi*(B3*tanh(C3/N3)+A3);                           % 2*pi*T1

% Calculation of t0
A4=0.5028; B4=0.0958;                                    % a4; b4
C4=0.02899*0.0113; D4=3.5667e-4; E4=0.1206;              %
N4=D4*H*L^3/E4;                                          % f4 
if xi~=0;                                                %
  N4=H*tanh(D4*xi*L^4)/tanh(E4*L*tanh(tanh(xi)));        %
end                                                      %
P4=H*tanh(C4*xi*L^3)-tanh(N4);                           % P4
t0=2*pi*(B4*P4+A4);                                      % 2*pi*T0

% Calculation of tmin
A5=0.9209; B5=-0.5623*4.1958;                            % a5; 4.1958*b5
N5=sqrt(abs(L-sqrt(Ur)+sqrt(2.5185/L)-4.6505)/H);        % g3
N6=abs(xi)+xi+L+0.9206;                                  % g1+g8
N7=abs(1.195*xi+0.195*abs(L+sqrt(3.0176/H)-5.2868+H));   % g7
N8=abs(abs(L+xi)+abs(xi)+xi-9.8976);                     % g5
tmin=2*pi*(B5/(N5+N6+N7+N8)+A5);                         % 2*pi*T2

% Calculation of U0 and corrections for large U0
U0=min(0.5*R,2*(t0-2*pi*(1-R))/(t0-tmax));               % U0
a1=max(0.99,(U0*tmax-4*pi*(1-R))/(t0*(U0-2)));           % a1

% Correction of t0, U1oUA, R and U0
t0=a1*t0;
if a1<0.99999
  U1oUA=U1oUA*(2*pi-0.25*(t0-tmax))*R/(2*pi-t0);
  R=(2*pi-t0)/(2*pi-0.25*(t0-tmax));
  U0=0.5*R;
end

% Calculation of u and a
umax=2*R;  umin=2*(1-R);                                 % UC/Us; UT/Us
u=zeros(size(t)); a=zeros(size(t));                      % U/Us

A1 =0.5*pi/tmax;        [ii1,tp] =find(t<=tmax);        tt1 =t(ii1)/tmax;
A01=0.5*pi/(t0-tmax);   [ii01,tp]=find(t>tmax & t<=t0); tt01=(t(ii01)-tmax)/(t0-tmax);
A02=0.5*pi/(tmin-t0);   [ii02,tp]=find(t>t0 & t<=tmin); tt02=(t(ii02)-t0)/(tmin-t0);
A2 =0.5*pi/(2*pi-tmin); [ii2,tp] =find(t>tmin);         tt2 =(t(ii2)-tmin)/(2*pi-tmin);

u(ii1) =umax*sin(0.5*pi*tt1);                    u(ii02)=-umin*sin(0.5*pi*tt02);
a(ii1) =A1*umax*cos(0.5*pi*tt1);                 a(ii02)=-A02*umin*cos(0.5*pi*tt02);
u(ii01)=umax*cos(0.5*pi*tt01)-U0*sin(pi*tt01);   u(ii2) =-umin*cos(0.5*pi*tt2);
a(ii01)=umax*sin(0.5*pi*tt01)+2*U0*cos(pi*tt01); a(ii01)=-A01*a(ii01);
a(ii2) =A2*umin*sin(0.5*pi*tt2);
