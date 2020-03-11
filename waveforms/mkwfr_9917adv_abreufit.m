clear
load('mat\9917adv_waveforms_burstwfcalcs_bp')
for ii=1:length(b)
    disp('.')
    if rem(ii,50)==0
        disp(ii)
    end
    wf=b(ii).wf;
    wfr(ii)=findwfr_abreufit(wf);
    
end

save mat\9917adv_wfr_abreufit_burstwfcalcs wfr