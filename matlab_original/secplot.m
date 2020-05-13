function [s1,s2,s3] = secplot(nc,ec,zc,figno,newlst,tstart,title1,title2,title3, mode, col)
%  [s1,s2,s3] = secplot(nc,ec,zc,figno,newlst,tstart,title1,title2,title3, mode, color)
% 

if ~exist('col')
    col='k';
end

ndt=0.025;
aflag=-1;
if mode==3
    figure(figno);clf;
    ns=size(newlst,1);
    % Plot sections
    s1=subplot(1,3,1);
    h=section(nc,tstart,ndt,aflag,newlst); set(h,'Color',col)
    title(title1)
    ylabel([' ']);
    ax1=gca;
    set(ax1,'YTickLabel',[newlst,repmat('.N',ns,1)],'YTick',[1:ns],'FontSize',14)
    xlim([tstart tstart+ndt*length(nc)])
    s2=subplot(1,3,2);
    h=section(ec,tstart,ndt,aflag,newlst); set(h,'Color',col)
    title(title2)
    ylabel([' ']);
    ax1=gca;
    set(ax1,'YTickLabel',[newlst,repmat('.E',ns,1)],'YTick',[1:ns],'FontSize',14)
    xlim([tstart tstart+ndt*length(nc)])
    s3=subplot(1,3,3);
    h=section(zc,tstart,ndt,aflag,newlst); set(h,'Color',col)
    ylabel([' ']);
    ax1=gca;
    set(ax1,'YTickLabel',[newlst,repmat('.Z',ns,1)],'YTick',[1:ns],'FontSize',14)
    title(title3)
    xlim([tstart tstart+ndt*length(nc)])
    
elseif mode==2
    figure(figno);clf;
    ns=size(newlst,1);
    % Plot sections
    s1=subplot(1,2,1);
    h=section(nc,tstart,ndt,aflag,newlst); set(h,'Color',col)
    title(title1)
    ylabel([' ']);
    ax1=gca;
    set(ax1,'YTickLabel',[newlst,repmat('.N',ns,1)],'YTick',[1:ns],'FontSize',14)
    xlim([tstart tstart+ndt*length(nc)])
    s2=subplot(1,2,2);
    h=section(ec,tstart,ndt,aflag,newlst); set(h,'Color',col)
    title(title2)
    ylabel([' ']);
    ax1=gca;
    set(ax1,'YTickLabel',[newlst,repmat('.E',ns,1)],'YTick',[1:ns],'FontSize',14)
    xlim([tstart tstart+ndt*length(nc)])
    s3=[];
end

return
