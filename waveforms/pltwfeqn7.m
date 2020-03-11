function wf7=mkwfeqn7(Uw,T,r,phi)

t=-1:T/128:T+1;
om=2*pi/T
for ii=1:length(t)
f=sqrt(1-r^2);    
ut(ii)=Uw*f*((sin(om.*t(ii))+(r*sin(phi)/(1+f)))/(1-r*cos(om.*t(ii)+phi)));
end
wf7=calcwfparams(ut,t);
ut=wf7.ub;
t=wf7.t;

if 1
    figure
    scrsz=get(0,'screensize')
    set(gcf,'position',[scrsz(3)+scrsz(3)*0.1 -50 scrsz(3)*0.8 scrsz(4)*0.8])
    set(gca,'xlim',[0 16])
    xlim=get(gca,'xlim');
    line(xlim,[0 0],'color',[0.75 0.75 0.75])
    hold on
    plot(t,ut)
    pltpts=[0 0;wf7.Tcu wf7.umax;wf7.Tc 0;wf7.Tc+wf7.Ttu wf7.umin;T 0];
    plot(pltpts(:,1),pltpts(:,2),'or','MarkerFaceColor','r')
end
%%
function wf=calcwfparams(ut,t)
%function wf=calcwfparams(ut,t)
%script to find waveform parameters from velocity time-series of a single
%wave form
%
%INPUT:
%ut = nearbed velocity time series(m/s)
%t = time (s)
%
%OUTPUT:
%wf - sturct of individual waveforms in burst
%   .t = wave form time (s)
%   .ub = wave nearbed velocity
%   .Tc = time crest (s)
%   .Tt = time trough (s)
%   .umax = amplitude crest (m/s)
%   .umin = amplitude trough (m/s)
%   .Ac = area crest (m/s * s)
%   .At = area trough (m/s *s)


try
    %find zero crossing times
    tz=findzs(ut,t);
   
    %loop to find crest/trough info
    for ii=1:length(tz)-1;
        %find half wave signal
        idx=find(t>tz(ii) & t<tz(ii+1));
        s(ii)={[0 ut(idx) 0]};
        ts(ii)={[tz(ii) t(idx) tz(ii+1)]};
        A(ii)=trapz(ts{ii},s{ii});
        S_max(ii)=max(s{ii});
        S_min(ii)=min(s{ii});
        T(ii)=tz(ii+1)-tz(ii);
    end

    cnt=0;
    xlim=[0 15];
    if sign(A(1))>=0
        for i=1:2:length(ts)-1;
        cnt=cnt+1;
        wf(cnt).t=[ts{i}-ts{i}(1) ts{i+1}(2:end)-ts{i}(1)];
        wf(cnt).ub=[s{i} s{i+1}(2:end)];
        wf(cnt).T=T(i)+T(i+1);
        wf(cnt).umax=S_max(i);
        wf(cnt).umin=S_min(i+1);
        wf(cnt).Tc=T(i);
        wf(cnt).Tt=T(i+1);
        wf(cnt).Tcu=wf(cnt).t(wf(cnt).ub==wf(cnt).umax);
        wf(cnt).Ttu=wf(cnt).t(wf(cnt).ub==wf(cnt).umin)-wf(cnt).Tc;
        wf(cnt).Ac=A(i);
        wf(cnt).At=A(i+1);
        wf(cnt).R=wf(cnt).umax/(wf(cnt).umax-wf(cnt).umin);
        wf(cnt).alpha=2*wf(cnt).Tcu/wf(cnt).T;
        wf(cnt).Su=mean(wf(cnt).ub.^3)/(std(wf(cnt).ub)).^3;

    %     plot(wf(cnt).t,wf(cnt).s)
    %     set(gca,'xlim',xlim)
    %     line(xlim, [0 0],'color',[0.5 0.5 0.5])
    %     set(gca,'ylim',[-1.5 1.5])
    %     pause(0.5)
    %     cla
        end

    else
        for i=2:2:length(ts)-2
        cnt=cnt+1;
        wf(cnt).t=[ts{i}-ts{i}(1) ts{i+1}(2:end)-ts{i}(1)];
        wf(cnt).ub=[s{i} s{i+1}(2:end)];
        wf(cnt).T=T(i)+T(i+1);
        wf(cnt).umax=S_max(i);
        wf(cnt).umin=S_min(i+1);
        wf(cnt).Tc=T(i);
        wf(cnt).Tt=T(i+1);
        wf(cnt).Tcu=wf(cnt).t(wf(cnt).ub==wf(cnt).umax);
        wf(cnt).Ttu=wf(cnt).t(wf(cnt).ub==wf(cnt).umin)-wf(cnt).Tc;
        wf(cnt).Ac=A(i);
        wf(cnt).At=A(i+1);
        wf(cnt).R=wf(cnt).umax/(wf(cnt).umax-wf(cnt).umin);
        wf(cnt).alpha=2*wf(cnt).Tcu/wf(cnt).T;
        wf(cnt).Su=mean(wf(cnt).ub.^3)/(std(wf(cnt).ub)).^3;
    %     plot(wf(cnt).t,wf(cnt).s)
    %     set(gca,'xlim',xlim)
    %     line(xlim, [0 0],'color',[0.5 0.5 0.5])
    %     set(gca,'ylim',[-1.5 1.5])
    %     pause(0.5)
    %     cla
        end
    end
catch
        wf.t=NaN;
        wf.ub=NaN;
        wf.T=NaN;
        wf.umax=NaN;
        wf.umin=NaN;
        wf.Tc=NaN;
        wf.Tt=NaN;
        wf.Tcu=NaN;
        wf.Ttu=NaN;       
        wf.Ac=NaN;
        wf.At=NaN;
        wf.R=NaN;
        wf.alpha=NaN;
        wf.Su=NaN;
end

return


%%
function tz=findzs(s,t)
%function tz=findzs(s,t)
%
%m-fcn to find zero crossing times of a signal s
%
%INPUT:
% s = signal
% t = time
len=length(s);

idx=1;
cnt=0;
c=sign(s(idx));
while idx < len
cnt=cnt+1;
c=sign(s(idx));
while c==sign(s(idx)) & idx < len
idx=idx+1;
end
tz(cnt)=interp1([s(idx-1) s(idx)],[t(idx-1) t(idx)],0);
end
return