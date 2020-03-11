%clear
load('..\..\..\data\Matanzas\proc\mat\11109vec_waveforms_burstwfcalcs_bp')
for ii=1:length(b)
    disp('.')
    if rem(ii,50)==0
        disp(ii)
    end
    wf=b(ii).wf;
    wfr(ii)=findwfr_abreufit_wt(wf);
    
end

save ..\..\..\data\Matanzas\proc\mat\11109vec_wfr_abreufit_burstwfcalcs_wt wfr