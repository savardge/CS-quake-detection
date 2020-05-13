function [h] = secplot_pretty(t, nc,ec,zc,figno,ordlst, tsN,tsE, tp, flagS,flagP)
addpath('/home/savardge/TOOLBOXES/Stransform/') 
addpath('/home/savardge/TOOLBOXES/') 
% addpath('/home/savardge/TOOLBOXES/samexaxis/')

h=figure(figno); clf; set(gcf, 'Color', [1,1,1]);
ns=size(nc,1);


for is=1:ns
    tsE(is)=tsE(is)-tsN(is);
    tp(is)=tp(is)-tsN(is);
    t(is,:)=t(is,:)-tsN(is);
    tsN(is)=tsN(is)-tsN(is);
end

for ii=1:ns
    
    subtightplot(ns,3,1+3*(ii-1)); plot(t(ii,:),nc(ii,:),'k'); axis tight; 
    hold on; h=vline(tsN(ii),'m','S'); if flagS(ii)==1; set(h,'LineWidth',8);uistack(h, 'bottom'); end; 
    if ii==ns
        set(gca,'YTick',0,'YTickLabel',[ordlst(ii,:)])
    else
        set(gca,'YTick',0,'XTickLabel','','YTickLabel',[ordlst(ii,:)])
    end
    if ii==1
        title('N')
    end
    
    subtightplot(ns,3,2+3*(ii-1)); plot(t(ii,:),ec(ii,:),'k'); axis tight; 
    hold on;  h=vline(tsE(ii),'m','S'); if flagS(ii)==1; set(h,'LineWidth',8); uistack(h, 'bottom');end; 
    if ii==ns
        set(gca,'YTick',0,'YTickLabel','')
    else
        set(gca,'YTick',0,'XTickLabel','','YTickLabel','')
    end
    if ii==1
        title('E')
    end
    
    subtightplot(ns,3,3+3*(ii-1)); plot(t(ii,:),zc(ii,:),'k'); axis tight; 
    hold on; vline(tsN(ii),'m--','S');
    hold on; h=vline(tp(ii),'m','P'); if flagP(ii)==1; set(h,'LineWidth',8); uistack(h, 'bottom'); end; 
    if ii==ns
        set(gca,'YTick',0,'YTickLabel','')
    else
        set(gca,'YTick',0,'XTickLabel','','YTickLabel','')
    end
    if ii==1
        title('V')
    end
    
    
end

  
return
