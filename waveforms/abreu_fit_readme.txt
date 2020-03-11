abreu_fit_readme.txt
Example using elfrinknd.m, abreu_fit.m and abreu_fstream.m MATLAB scripts as
described in Malarkey and Davies (2012). 

Consider the a wave with 100 points in the wave cycle where H = 0.6, 
L = 15, slope = 1/40

(1) Calculating U1oUA, R, tmax, t0, tmin and U0 and u and a using Elfrink et al.
(2006) characterisation

>> t=[0:99]'*pi/50; H=0.6; L=15; slope=1/40;
>> xi=slope*sqrt(L/(H*tanh(2*pi/L)));
>> [U1oUA,R,tmax,t0,tmin,U0,u,a]=elfrinknd(t,H,L,xi);

(2) Fitting to the ABR description to u on a least-squares-fit basis to find b1

>> [b1ls]=abreu_fit(R,t,u);
>> [ufls,afls,bls,phils,R,betals,alphals]=abreu_fstream(t,R,b1ls,1,2*pi,3,2,2);

(3) Fitting to the ABR description on the basis of tmax and tmin (b1(1) and b1(2)
respectively).

>> [b1]=abreu_fit(R,[tmax; tmin],[2*R; 2*(R-1)]);
>> [uf11,af11,b11,phi11,R,beta11,alpha11]=abreu_fstream(t,R,b1(1),1,2*pi,3,2,2);
>> [uf12,af12,b12,phi12,R,beta12,alpha12]=abreu_fstream(t,R,b1(2),1,2*pi,3,2,2);

The abreu_fit.m function will work on any continuous, non-dimensional, free-stream 
velocity and time provided that 0 <= t < 2pi and umax-umin = 2 (umax = 2R and umin 
= 2(R-1)) and umax occurs before umin and t = 0 corresponds to the zero up-crossing.
Further information for each function is available in the function header.

Malarkey, J. and Davies, A.G. 2012. Free-stream velocity descriptions under waves
with skewness and asymmetry. Coastal Engineering, 68, 78-95.


Jonathan Malarkey 2012
