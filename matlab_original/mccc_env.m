function [tdel,rmean,sigr,r,tcc,tdelw,resw,sigrw] = mccc_env(seis,lenw, twinmin,twinmax,weight)
% Function MCCC determines optimum relative delay times for a
% set of seismograms based on the VanDecar & Crosson multi-channel
% cross-correlation algorithm. 

ndt=0.025;

% Zero-mean
% seis=seis-repmat(mean(seis,2),1,size(seis,2));

% Set nt to twice length of seismogram section to avoid
% spectral contamination/overlap. Note we assume that
% columns enumerate time samples, and rows enumerate stations.
% nt=size(seis,2)*2;
nt = size(seis,2);
ns=size(seis,1);
tcc=zeros(ns,ns);

% Window on which actual correlation is performed (exclude "buffing" before
% and after
imid = (nt-1)/2 +1; % middle sample if nt is odd
iwin = (imid-lenw/ndt) : (imid+lenw/ndt); % +/- lenw seconds around the middle sample

%% Cross-correlation
r=zeros(ns,ns);

for is=1:ns-1
   
    for js=is+1:ns
        
        xa=seis(is,iwin); %-mean(seis(is,iwin));
        na=length(xa); % MUST BE ODD
        xb=seis(js,:); % remove mean later

        % Set mask (Allow for itw longer than length of correlation.)
   
        if (twinmin(is,js)<0 && twinmax(is,js)<0)
            % both negative lags. twinmin > twinmax because negative
            itwmin=min((nt-1)/2-(na-1)/2,ceil(abs(twinmin(is,js))./ndt));
            itwmax=min((nt-1)/2-(na-1)/2,ceil(abs(twinmax(is,js))./ndt));
            nminlag = imid - itwmin;
            nmaxlag = imid - itwmax;
        elseif (twinmin(is,js)>0 && twinmax(is,js)>0)
            % both positive lags. twinmin < twinmax because positive         
            itwmin=min((nt-1)/2-(na-1)/2 - na,ceil(abs(twinmin(is,js))./ndt));
            itwmax=min((nt-1)/2-(na-1)/2 - na,ceil(abs(twinmax(is,js))./ndt));
            nminlag = imid + itwmin;
            nmaxlag = imid + itwmax;
        elseif (twinmin(is,js)<0 && twinmax(is,js)>0)
            % min lag is negative, max lag is positive
            itwmin=min((nt-1)/2-(na-1)/2, ceil(abs(twinmin(is,js))./ndt));
            itwmax=min((nt-1)/2-(na-1)/2, ceil(abs(twinmax(is,js))./ndt));
            nminlag = iwin(1) - itwmin;
            nmaxlag = iwin(1) + itwmax;
            lags=nminlag:nmaxlag;
            lagslab=-(itwmin*ndt):ndt:(itwmax*ndt);
        end
        
        % CC
        aa=sqrt(xa*xa');
        ccf=zeros(1,length(lags));
        
        for ik=1:length(lags)
            ia=lags(ik);
            xb2=xb(ia:ia+na-1); %-mean(xb(ia:ia+na-1));
            bb=sqrt(xb2*xb2');
            ccf(ik)=xa*xb2'/(aa*bb);
        end
        iz=~isfinite(ccf);
        ccf(iz)=0;

        % Compute estimate of cross correlation coefficient.
        [cmax,ix]=max(ccf);
        tcc(is,js)=lagslab(ix);
        r(is,js)=cmax;

        %% Display
        %     if shw==1
%         figure(3);clf
%         plot(seis(is,:));hold on;plot(seis(js,:),'r');
%                 figure(4);clf
% %                 plot(1:nt,padarray(xa,[0 iwin(1)-1],'both'),1:nt,xb)
%                 plot(xa);hold on;plot(xb(lags(ix):lags(ix)+na-1),'r')
%                 figure(5); clf;
%                 plot(lagslab,ccf);
%                 hold on
%                 vline(tcc(is,js),'r')
%                 pause;
% %         %     end
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
