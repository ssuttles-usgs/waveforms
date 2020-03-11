clear
wh=load('C:\Users\ssuttles\data\Matanzas\proc\mat\workhorse_matanzas_inlet.mat')
ave=load('C:\Users\ssuttles\data\Matanzas\proc\mat\11109vec_wfr_abreufit_burstwfcalcs.mat')
wt=load('C:\Users\ssuttles\data\Matanzas\proc\mat\11109vec_wfr_abreufit_burstwfcalcs_wt.mat')
%sig=load('mat\9917adv_wfr_abreufit_burstwfcalcs_sig.mat')
%load('mat\ur_spec_Umo_Tp.mat')

Uwh=[wh.uhat_emp(1:1893)];
%Usig=[sig.wfr.Uw];
Uwt=[wt.wfr.Uw];
Uave=[ave.wfr.Uw];
Uplt=[Uwh' Uave' Uwt']

figure
scrsz=get(0,'screensize')
set(gcf,'position',[100 100 scrsz(3)*0.8 scrsz(4)*0.8])
plot(Uplt)
%set(gca,'xlim',[100 250])
legend({'U_{wh}' 'U_{ave}' 'U_{wt}'},'box','off','fontsize',12)
title('Matanzas- SiteA Wave Orbital Velocity Amplitude Comparison','fontsize',12)
xlabel('burst')
ylabel('wave orbital amplitude, [m ^. s^{-1}]')
save2pdf('pdfs\Uall_compare_matan',gcf,300)