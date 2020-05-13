function [tdel,rmean,sigr,r,tcc,tdelw,resw,sigrw] = mccc_target(seis,twinmin,twinmax,weight)
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
nt=size(seis,2)*2;
ns=size(seis,1);
tcc=zeros(ns,ns);

% Determine relative delay times between all pairs of
% traces.

%% First compute autocorrelations.
acf=real(ifft(fft(seis,nt,2).*conj(fft(seis,nt,2)),nt,2));
sigt=sqrt(max(acf,[],2));

%% Cross-correlation
r=zeros(ns,ns);

for is=1:ns-1
    %   ffis=conj(fft(seis(is,:),nt));
    for js=is+1:ns
        % Mask (Set width of window around 0 time to search for maximum.)
        % Allow for itw longer than length of correlation.
        
        itwmin=min(nt/2,ceil(abs(twinmin(is,js))./(2*ndt)));
        itwmax=min(nt/2,ceil(abs(twinmax(is,js))./(2*ndt)));
        mask=zeros(1,nt);
        if (twinmin(is,js)<0 && twinmax(is,js)<0)
            % both negative lags. twinmin > twinmax because negative
            mask(nt-itwmin+1:nt-itwmax)=1.0;
        elseif (twinmin(is,js)>0 && twinmax(is,js)>0)
            % both positive lags. twinmin < twinmax because positive
            mask(itwmin:itwmax+1)=1.0; 
        elseif (twinmin(is,js)<0 && twinmax(is,js)>0)
            % min lag is negative, max lag is positive
            mask(1:itwmax+1)=1.0; % positive lags
            mask(nt-itwmin+1:end)=1.0; % negative lags
        end
      
        % CC
        % ffjs=fft(seis(js,:),nt);
        ccfu=real(ifft((conj(fft(seis(is,:),nt)).*fft(seis(js,:),nt)),nt));
        ccf=ccfu.*mask;        
        [cmax,tcc(is,js)]=max(ccf);

        % Compute estimate of cross correlation coefficient.
        r(is,js)=cmax/(sigt(is)*sigt(js));

        %% Display
        %     if shw==1
%                 sig=std(ccf([1:itw,nt-itw:nt]));
%                 figure(5); clf;
%                 hold on;
%                 plot((-(nt-1)/2:(nt-1)/2).*ndt,fftshift(ccfu));
%                 plot((-(nt-1)/2:(nt-1)/2).*ndt,fftshift(ccf),'r');axis tight %xlim(gca,[-2*itw*ndt,2*itw*ndt]);
%                 if tcc(is,js)<(nt-1)/2
%                     plot(tcc(is,js)*ndt,cmax,'ro');
%                 else
%                     plot(-(nt-tcc(is,js))*ndt,cmax,'mo');
%                 end
%                 hold off;
%                 figure(6);clf;
%                 subplot(2,1,1);plot((0:(nt-1)/2).*ndt,seis(is,:));
%                 subplot(2,1,2);plot((0:(nt-1)/2).*ndt,seis(js,:));
%                 pause;
%         %     end
        %
        %%
        
    end
end

% Fisher's transform of cross-correlation coefficients to produce
% normally distributed quantity on which Gaussian statistics
% may be computed and then inverse transformed.
z=0.5*log((1+r)./(1-r));
zmean=(sum(z,2)'+sum(z,1))/(ns-1);
rmean=(exp(2*zmean)-1)./(exp(2*zmean)+1); % inverse fischer transform

% Correct negative delays.
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
for is=1:ns-1
    for js=is+1:ns
        res(is,js)=tcc(is,js)-(tdel(is)-tdel(js)); % eq. 7 from paper
    end
end

% Standard deviation of the residual
sigr=zeros(1,ns);
for is=1:ns
   sigr(is)=sqrt((sum(res(is,is+1:ns).^2)+sum(res(1:is-1,is).^2))/(ns-2));
end

%% Weighted delay times

if weight==1
    % Coefficient matrix
    k=0;
    A=zeros(ns,ns);
    for is=1:ns-1
        for js=is+1:ns
            A(is+k,is)=1;
            A(is+k,js)=-1;
            k=k+1;
        end
    end
    
    A(all(A==0,2),:)=[]; % clean rows of zeros
    A= [A; ones(1,size(A,2))]; % add row of ones
    
    % Weight matrix
    v=zeros((ns*ns-ns)/2,1);
    k=1;
    for i=1:ns-1
        for j=i+1:ns
            v(k)=r(i,j).^(2); %sqrt(sigr(i).^2 + sigr(j).^2);
            k=k+1;
        end
    end
    v(~isfinite(v))=0;
    v=v./max(v); % Normalize
    v=[v;1];
    W=diag(v,0);
   
   
    
    
    % Vectorize tcc
    u=zeros((ns*ns-ns)/2,1);
    k=1;
    for i=1:ns-1
        for j=i+1:ns
            u(k)=tcc(i,j);
            k=k+1;
        end
    end
    tccv=[u;0];
    
    % Calculate weighted delay times
    tdelw=((A'*W*A)\A'*W)*tccv;
    
    % Compute associated residuals (weighted)
    resw=zeros(ns,ns);
    for is=1:ns-1
        for js=is+1:ns
            resw(is,js)=tcc(is,js)-(tdelw(is)-tdelw(js)); % eq. 7 from paper
        end
    end
    
    % Standard deviation of the weighted residual
    sigrw=zeros(1,ns);
    for is=1:ns
        sigrw(is)=sqrt((sum(resw(is,is+1:ns).^2)+sum(resw(1:is-1,is).^2))/(ns-2));
    end
end

end
