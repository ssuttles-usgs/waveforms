clear;
load('mat\9917adv_waveforms_burstwfcalcs_bp.mat')
len=length(b(1).UBS.ur);
fs=8;
nseg=16;
pow2=nextpow2(8400/nseg);
pow2=nextpow2(8400/nseg)
nfft=2^pow2
nx = len;
na = nx/nfft;
w=hanning(nfft)
for ii=1:2044
 disp(ii)
u=b(ii).UBS.ur;
[Uxx,F]=pwelch(u,w,0,nfft,8);
m0=trapz(F,Uxx);
Umo(ii)=4*sqrt(m0)/2;

%find Tp
 iTp=find(Uxx == max(Uxx));
 FF=F;
 
try
 Tp(ii)=1./FF(iTp(1));
catch
    Tp(ii)=NaN;
end
end
save mat\ur_spec_Um0_Tp Umo Tp