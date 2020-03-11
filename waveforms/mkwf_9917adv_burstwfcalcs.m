clear
load ..\mat\puv_proc_FI_iwaves2 raw depth fs zr zp dn
%load('..\mat\puv_proc_FI_iwaves_depc.mat')
for ii=1:length(depth)
  %make input data to m-fcn
  ii
datin.bn=ii;
datin.dn=dn(ii);
datin.depth=depth(ii);
datin.fs=fs;
datin.zr=zr(ii)
datin.zp=zp(ii)
datin.u=raw(ii).u;
datin.v=raw(ii).v;
datin.p=raw(ii).p;
datin.t=[0:1/datin.fs:(length(datin.u)-1)./datin.fs]';

%set filttype
filttype=1; %remove ig

[b(ii).wf,b(ii).UBS,b(ii).PUV]=burstwfcalcs(datin,filttype);  
%b(ii).wf=findurwaveform(UBS(ii).ur',fs,dn(ii))
end

save mat\9917adv_waveforms_burstwfcalcs b