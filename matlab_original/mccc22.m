function [tdel,rmean,sigr,r,tcc] = mccc22(seis,twin)
% Function MCCC determines optimum relative delay times for a
% set of seismograms based on the VanDecar & Crosson multi-channel
% cross-correlation algorithm. SEIS is the set of seismograms. It
% is assumed that this set includes the window of interest and
% nothing more since we calculate the correlation functions in the
% Fourier domain. DT is the sample interval and TWIN is the window
% about zero in which the maximum search is performed (if TWIN is
% not specified, the search is performed over the entire correlation
% interval).
% If parameter shw=1, the cross correlation between waveforms will be
% displayed.
ndt=0.025;

seis=seis-repmat(mean(seis,2),1,size(seis,2));

% Set nt to twice length of seismogram section to avoid
% spectral contamination/overlap. Note we assume that
% columns enumerate time samples, and rows enumerate stations.
nt = 2^nextpow2(2*size(seis,2)-1);
ns=size(seis,1);
tcc=zeros(ns,ns);

% Temporary fix to make Genevieve's version compatible with original.
if size(twin) == [1,1]
    twin=ones(ns,ns).*twin;
end

% Determine relative delay times between all pairs of
% traces.
ftseis=fft(seis,nt,2);

%% First compute autocorrelations.
acf=real(ifft(ftseis.*conj(ftseis),nt,2)); %./(nt/2);
% acf=seis.*seis;
sigt=sqrt(max(acf,[],2));

%% Cross-correlation
r=zeros(ns,ns);

% [pair(:,1),pair(:,2)]=ind2sub(ns,find(triu(ones(ns,ns),1)));

%%%%%%%%%

for is=1:ns-1
    L=numel(is:ns);
    ccfu=real(ifft((conj(ftseis(ones(1,L).*is,:)).*ftseis(is:ns,:)),nt,2)); %./(nt/2);
    
    for js=1:(L-1) %is+1:ns
        itw=min(nt/2,round(twin(is,is+js)./(2*ndt)));
        mask=zeros(1,nt); mask(:,[(1:itw+1)';(nt-itw+1:nt)'])=1;
        [cmax,tcc(is,is+js)]=max(ccfu(js+1,:).*mask,[],2);
        
        figure(30);
        subplot(3,1,1); plot(seis(is,:))
        subplot(3,1,2); plot(seis(js+1,:))
        subplot(3,1,3); plot(fftshift(ccfu(js+1,:))./(sigt(is).*sigt(is+js)));ylim([0 0.5])
        sigt(is)
        sigt(js)
        pause

        % Compute estimate of cross correlation coefficient.
        r(is,is+js)=cmax./(sigt(is).*sigt(is+js));
        
    end
    tcc(is,is)=0;
    r(is,is)=0;
end



%%%%%%
%
% for k=1:size(pair,1)
%         ccfu=real(ifft((conj(ftseis(pair(k,1),:)).*ftseis(pair(k,2),:)),nt));
% %         mask=zeros(1,nt);
% %         mask([(1:itw(pair(k,1),pair(k,2))+1)';(nt-itw(pair(k,1),pair(k,2))+1:nt)'])=1;
% %         [cmax,tcc(pair(k,1),pair(k,2))]=max(ccfu.*mask);
%      [cmax,tcc(pair(k,1),pair(k,2))]=max(ccfu([(1:itw(pair(k,1),pair(k,2))+1)';(nt-itw(pair(k,1),pair(k,2))+1:nt)']));
%
%     if tcc(pair(k,1),pair(k,2))>(itw(pair(k,1),pair(k,2))+1)
%         tcc(pair(k,1),pair(k,2))=tcc(pair(k,1),pair(k,2))-(itw(pair(k,1),pair(k,2))+1);
%     end
%
%         r(pair(k,1),pair(k,2))=cmax/(sigt(pair(k,1))*sigt(pair(k,2)));
%
% end

%%%%%%%%%%%%%

% Fisher's transform of cross-correlation coefficients to produce
% normally distributed quantity on which Gaussian statistics
% may be computed and then inverse transformed.
z=0.5*log((1+r)./(1-r));
zmean=(sum(z,2)'+sum(z,1))/(ns-1);
rmean=(exp(2*zmean)-1)./(exp(2*zmean)+1); % inverse fischer transform

%Correct negative delays.
ix=find(tcc>nt/2);
tcc(ix)=tcc(ix)-nt;

% Subtract 1 to account for sample 1 at 0 lag.
tcc=tcc-1;

% Multiply by sample rate.
tcc=tcc*ndt;

% Use sum rule to assemble optimal delay times with zero mean.
% tdel=(-sum(tcc,2)'+sum(tcc,1))/(ns);
tdel=zeros(1,ns);
for is=1:ns
    tdel(is)=(-sum(tcc(1:is-1,is))+sum(tcc(is,is+1:ns)))/ns; % eq. 6 from paper
end

% Compute associated residuals (unweighted).
res=zeros(ns,ns);
% for is=1:ns-1
%     for js=is+1:ns
%         res(is,js)=tcc(is,js)-(tdel(is)-tdel(js)); % eq. 7 from paper
%     end
% end

for is=1:ns-1
    res(is,(is+1):ns)=tcc(is,(is+1):ns)'-(tdel(is)-tdel((is+1):ns)');
end

% Standard deviation of the residual
sigr=zeros(1,ns);
for is=1:ns
    sigr(is)=sqrt((sum(res(is,is+1:ns).^2)+sum(res(1:is-1,is).^2))/(ns-2));
end

return