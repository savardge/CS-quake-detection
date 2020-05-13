function [] = plot_wf_powspec(ncomp,ecomp,zcomp,stalst,sdfm,col)

% Standardize input
ncomp = zscore(ncomp')';
ecomp = zscore(ecomp')';
zcomp = zscore(zcomp')';

ndt=0.025;
ns=size(ncomp,1);
nfftne=2^nextpow2(size(ncomp,2));
nfftz=2^nextpow2(size(zcomp,2));
[p,fne] = periodogram(ncomp(1,:),hamming(size(ncomp,2)),nfftne,1/ndt,'power');
psn = zeros(ns,length(fne)); pse=psn;
[p,fz] = periodogram(zcomp(1,:),hamming(size(zcomp,2)),nfftz,1/ndt,'power'); 
psz=zeros(ns,length(fz));

for is=1:ns
    [psn(is,:),f] = periodogram(ncomp(is,:),hamming(size(ncomp,2)),nfftne,1/ndt,'power');
    [pse(is,:),f] = periodogram(ecomp(is,:),hamming(size(ncomp,2)),nfftne,1/ndt,'power');
    [psz(is,:),f] = periodogram(zcomp(is,:),hamming(size(zcomp,2)),nfftz,1/ndt,'power');
end
ix=find(fne>10,1,'first');
fne=fne(1:ix);
psn=psn(:,1:ix);
pse=pse(:,1:ix);
ix=find(fz>10,1,'first');
fz=fz(1:ix);
psz=psz(:,1:ix);

subplot(1,3,1)
if col == 'r'; hold on; else; hold off; end
h=section(psn,fne(1),fne(2)-fne(1),-1,stalst);
set(h,'Color',col);
ylabel([' ']);
xlabel('Frequency Hz')
xlim([1 10])
ylim([0,ns+1])
set(gca,'YTickLabel',stalst,'YTick',[1:ns])
subplot(1,3,2)
if col == 'r'; hold on; else; hold off; end
h=section(pse,fne(1),fne(2)-fne(1),-1,stalst);
set(h,'Color',col);
title(['Time: ',sdfm])
ylabel([' ']);
xlabel('Frequency Hz')
xlim([1 10])
ylim([0,ns+1])
set(gca,'YTickLabel',stalst,'YTick',[1:ns])
subplot(1,3,3)
if col == 'r'; hold on; else; hold off; end
h=section(psz,fz(1),fz(2)-fz(1),-1,stalst);
set(h,'Color',col);
ylabel([' ']);
xlabel('Frequency Hz')
xlim([1 10])
ylim([0,ns+1])
set(gca,'YTickLabel',stalst,'YTick',[1:ns])

figure(30)
subplot(1,3,1); plot(fne,sum(psn,1)./ns);title('N component stack of PSD')
subplot(1,3,2); plot(fne,sum(pse,1)./ns);title('E component stack of PSD')
subplot(1,3,3); plot(fz,sum(psz,1)./ns);title('Z component stack of PSD')

return
