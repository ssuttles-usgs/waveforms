clear
load C:\Users\ssuttles\data\Matanzas\proc\mat\puv_proc_Matan_iwaves raw depth fs zr zp dn
%load('..\mat\puv_proc_FI_iwaves_depc.mat')
for ii=1:length(depth)
  %make input data to m-fcn
  ii
  
  %u check for and interpret NaN
      iu=find(~isnan(raw(ii).u));
      ui=interp1(iu,raw(ii).u(iu),[1:length(raw(ii).u)],'linear','extrap')';
      
      %v check for and interpret NaN
      iv=find(~isnan(raw(ii).v));
      vi=interp1(iv,raw(ii).v(iv),[1:length(raw(ii).v)],'linear','extrap')';
 
datin.bn=ii;
datin.dn=dn(ii);
datin.depth=depth(ii);
datin.fs=fs;
datin.zr=zr(ii)
datin.zp=zp(ii)
datin.u=ui;
datin.v=vi;
datin.p=raw(ii).p;
datin.t=[0:1/datin.fs:(length(datin.u)-1)./datin.fs]';



%set filttype
filttype=2; %remove ig using bandpass (butter 2nd order)

[b(ii).wf,b(ii).UBS,b(ii).PUV]=burstwfcalcs(datin,filttype);  
%b(ii).wf=findurwaveform(UBS(ii).ur',fs,dn(ii))
end

save ..\..\..\data\Matanzas\proc\mat\11109vec_waveforms_burstwfcalcs_bp b