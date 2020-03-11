wh=load('mat\9921wh_wfr.mat')
ave=load('mat\9917adv_wfr_abreufit.mat')
wt=load('mat\9917adv_wfr_abreufit_burstwfcalcs_wt.mat')
sig=load('mat\9917adv_wfr_abreufit_burstwfcalcs_sig.mat')
load('mat\Umo.mat')

Uwh=[wh.wfr.Uw];
Usig=[sig.wfr.Uw];
Uwt=[wt.wfr.Uw];
Uave=[ave.wfr.Uw];
Uplt=[Uwh' Uave' Uwt' Usig' Umo']

figure
scrsz=get(0,'screensize')
set(gcf,'position',[100 100 scrsz(3)*0.8 scrsz(4)*0.8])
plot(Uplt)
set(gca,'xlim',[100 250])
legend({'U_{wh}' 'U_{ave}' 'U_{wt}' 'U_{sig}' 'U_{m0}'},'box','off','fontsize',12)
title('FI14- Site3 Wave Orbital Velocity Amplitude Comparison- during mid-Feb event','fontsize',12)
xlabel('burst')
ylabel('wave orbital amplitude, [m ^. s^{-1}]')
save2pdf('Uall_compare',gcf,300)