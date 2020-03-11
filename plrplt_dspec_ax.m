function hp=plrplt_dspec_ax(ax,S,tstr,ca,cb)
%function hp=plrplt_dspec_ax(ax,S,tstr,ca,cb)
%INPUT
%S- dspec + wavee stats  (struct):
%   S.dspec-  directional spectra [m^2/s^/Hz}
%   S.freqs- frequencies [Hz]
%   S.dirs- directionas [degT] FROM convention
%   tstr- title of plot (string)
%   ca- axis limits in log10 units; i.e. 2 = 100, 2-element vector [lim1 lim2]

freqs=S.freqs;
dirs=S.dirs;
Dspec=S.dspec;
ineg=find(Dspec <= 0);
Dspecplt=Dspec;
Dspecplt(ineg)=1e-9;

if exist('ca')
    clim1=ca(1);
    clim2=ca(2);
else
clim1=-8;
clim2=3;
end

% hf1=figure
% fpos=get(hf1,'position')
% set(hf1,'position',[fpos(1)*0.75 fpos(2)*0.5 fpos(3)*1.33 fpos(4)*1.33])

hp=polar(pi,20)
set(get(hp,'Parent'),'fontweight','bold')
set(get(hp,'Parent'),'fontsize',12)
set(gca,'view',[90 270])
set(gca,'linewidth',1.5)
get(gcf,'currentaxes')
hpos=get(ax,'position')
set(ax,'position',[hpos(1) hpos(2) hpos(3)*0.9 hpos(4)*0.9])

upperlim=ceil(log10(max(max(Dspecplt))))
rem((log10(max(max(Dspecplt)))),1)
%lowlim=upperlim-2
lowlim=clim1;


fhot=hot;
fhot=fhot(64:-1:1,:);
%colormap('default')

caxis([clim1 clim2]);

%caxis('manual');
cmap=colormap;
%bclr=get(gcf,'color');

if exist('cb') & cb==1
h=colorbar('WestOutside');
set(h,'ylim',[clim1 clim2]);
set(h,'ytick',[clim1:1:clim2]);
%set(h,'yticklabel',[10.^get(h,'ytick')])
%set(gcf,'currentaxes',h);
hpos=get(h,'position')
set(h,'position',[hpos(1)-0.1 hpos(2) hpos(3)/2 hpos(4)])
hxl=xlabel('log10([m^2 ^. s ^{-1} ^. deg^{-1}])','fontsize',10,'fontweight','bold');
hxl_pos=get(hxl,'position')
set(hxl,'position',[hxl_pos(1) hxl_pos(2) hxl_pos(3)])
end

%axes(ax)
%title('Directional Wave Spectra Polar Plot')


if isnan(Dspec)
    
    x0=7*cos(315*pi/180);
    y0=7*sin(315*pi/180);
    text(x0,y0,'No Directional Data')

else
    
[r,c]=size(Dspec);
%x=fliplr(dirs+S.xaxisdir)+180;
% idx=find(x > 360);
% x(idx)=x(idx)-360;
x=dirs;
y=freqs;
[X,Y]=meshgrid(x,y');

hold on;
    
[px,py]=pol2cart(X*pi/180,1./Y);
contourf(px',py',log10(Dspecplt),[lowlim:0.05:upperlim]);
%pcolor(px,py,log10(Dspecplt))
shading interp

end


caxis([clim1 clim2]);
tstr=sprintf('Directional Wave Spectra: %s',tstr)
text(25,0,tstr,'fontsize',10,'fontweight','bold','horizontalalignment','center')

if 0
axtxt=axes('units','pixels','position',[fpos(3)*0.75 fpos(4)*0.05 fpos(3)*0.2 fpos(4)*0.2],'color','none','visible','off')
set(axtxt,'ydir','reverse')
text(1,0,sprintf('%s',datestr(S.dn,0)),'fontsize',8,'fontweight','bold','horizontalalignment','right')
text(1,0.2,sprintf('Burst # = %d',S.nb),'fontsize',8,'fontweight','bold','horizontalalignment','right')
text(1,0.4,sprintf('Hsig = %5.2f m',S.Hs),'fontsize',8,'fontweight','bold','horizontalalignment','right')
text(1,0.6,sprintf('Tp = %6.2f s',S.Tp),'fontsize',8,'fontweight','bold','horizontalalignment','right')
text(1,0.8,sprintf('DTp = %d ^oT',S.DTp),'fontsize',8,'fontweight','bold','horizontalalignment','right')
%text(1,1,sprintf('Dp = %d ^oT',Dp),'fontsize',8,'fontweight','bold','horizontalalignment','right')
%disp('You are Here')
%pause
end

if 0
axwnd=axes('units','pixels','position',[fpos(3)*0.75 fpos(4)*0.65 fpos(3)*0.25 fpos(4)*0.25],'color','none','visible','off')
hp=quiver(0,0,wnde,wndn)
set(hp,'linewidth',1.5)
set(hp,'MarkerSize',5)
set(hp,'MaxHeadSize',5)
set(gca,'xlim',[-15 15])
set(gca,'ylim',[-15 15])
set(gca,'visible','off')
set(hp,'Autoscale','off')
set(hp,'linewidth',1.5)
set(hp,'MarkerSize',10)
set(hp,'MaxHeadSize',5)
line([-5 5],[-15 -15],'color','k','linewidth',1)
line([-5 -5],[-15 -14],'color','k','linewidth',1)
line([5 5],[-15 -14],'color','k','linewidth',1)
text(0,-16.5,'Scale= 10 [m ^. s^{-1}]','horizontalalignment','center','fontsize',8,'fontweight','bold')
text(0,15,wstr,'horizontalalignment','center','fontsize',8,'fontweight','bold')
end
