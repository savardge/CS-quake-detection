function [stnproj,steproj,stzproj] = secplot_st(nc,ec,zc,figno,ordlst, mode)
addpath('/home/savardge/TOOLBOXES/Stransform/')
addpath('/home/savardge/TOOLBOXES/samexaxis/')

figure(figno); clf
ns=size(nc,1);
maxfreq=140;
if mode==3
    for ii=1:ns
        [stn,t,f] = st(nc(ii,:),0,maxfreq,0.025,1); stn=abs(stn);
        [ste,t,f] = st(ec(ii,:),0,maxfreq,0.025,1); ste=abs(ste);
        [stz,t,f] = st(zc(ii,:),0,maxfreq,0.025,1); stz=abs(stz);
        stnproj(ii,:) = sum(abs(stn),1)./max([max(stn(:)),max(ste(:)),max(stz(:))]);
        steproj(ii,:) = sum(abs(ste),1)./max([max(stn(:)),max(ste(:)),max(stz(:))]);
        stzproj(ii,:) = sum(abs(stz),1)./max([max(stn(:)),max(ste(:)),max(stz(:))]);
        subplot(ns,3,1+3*(ii-1)); imagesc(t,f,abs(stn));axis xy;colorbar;
        subplot(ns,3,2+3*(ii-1)); imagesc(t,f,abs(ste));axis xy;colorbar;
        subplot(ns,3,3+3*(ii-1)); imagesc(t,f,abs(stz));axis xy;colorbar;
        
    end
    samexaxis('abc','xmt','on','ytac','join','yld',1)
    secplot(stnproj,steproj,stzproj,figno,ordlst,0,'N','E','Z', 3);
elseif mode==2
    for ii=1:ns
        [stn,t,f] = st(nc(ii,:),0,maxfreq,0.025,1); stn=abs(stn);
        [ste,t,f] = st(ec(ii,:),0,maxfreq,0.025,1); ste=abs(ste);
        stnproj(ii,:) = sum(abs(stn),1)./max([max(stn(:)),max(ste(:))]);
        steproj(ii,:) = sum(abs(ste),1)./max([max(stn(:)),max(ste(:))]);
        subtightplot(ns,2,1+2*(ii-1)); imagesc(t,f,abs(stn));axis xy;colorbar;
        subtightplot(ns,2,2+2*(ii-1)); imagesc(t,f,abs(ste));axis xy;colorbar;
    end
%     samexaxis('abc','xmt','on','ytac','join','yld',1)
    secplot(stnproj,steproj,[],figno+1,ordlst,0,'N','E','Z', 2);
end



return
