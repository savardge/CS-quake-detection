function [] = secplot_psd(nc,ec,zc,figno,newlst,title1,title2,title3, mode)

ndt=0.025;
aflag=-1;
if mode==3
    figure(figno);
    ns=size(newlst,1);
    % Plot sections
    subplot(1,3,1)
    section_psd(nc,ndt,newlst);
    title(title1)
    ylabel([' ']);
    ax1=gca;
    set(ax1,'YTickLabel',newlst,'YTick',[1:ns],'FontSize',14)
    xlim auto
    subplot(1,3,2)
    section_psd(ec,ndt,newlst);
    title(title2)
    ylabel([' ']);
    ax1=gca;
    set(ax1,'YTickLabel',newlst,'YTick',[1:ns],'FontSize',14)
    xlim auto
    subplot(1,3,3)
    section_psd(zc,ndt,newlst);
    ylabel([' ']);
    ax1=gca;
    set(ax1,'YTickLabel',newlst,'YTick',[1:ns],'FontSize',14)
    title(title3)
    xlim auto
    
elseif mode==2
    figure(figno);
    ns=size(newlst,1);
    % Plot sections
    subplot(1,2,1)
    section_psd(nc,ndt,newlst);
    title(title1)
    ylabel([' ']);
    ax1=gca;
    set(ax1,'YTickLabel',newlst,'YTick',[1:ns],'FontSize',14)
    xlim auto
    subplot(1,2,2)
    section_psd(ec,ndt,newlst);
    title(title2)
    ylabel([' ']);
    ax1=gca;
    set(ax1,'YTickLabel',newlst,'YTick',[1:ns],'FontSize',14)
    xlim auto
end

return
