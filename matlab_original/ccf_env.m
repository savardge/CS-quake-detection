function [ccf] = ccf_env(ncompa,ncompb,twinmin,twinmax)

% Simple function to compute running correlation coefficient
% of short window (master envelope = ncompa) time series across long window 
% (other envelopes = ncompb) time series. 
%   twinmin: lower bound for negative lags
%   twinmax: upper bound for positive lags
% Could be made more efficient.

ncompa=ncompa-mean(ncompa);
na=length(ncompa); % MUST BE ODD
nb=length(ncompb); % MUST BE ODD
if rem(na,2)==0 || rem(nb,2)==0
    error('Waveform lengths must be odd!')
end

nminlag=(twinmin/ndt); 
nmaxlag=(twinmax/ndt);
nzero=(na-1)/2 +1;

lags=(nzero-nminlag):(nzero+nmaxlag);
% lagslab=-nminlag:nmaxlag;
aa=sqrt(ncompa*ncompa');
ccf=zeros(1,length(lags));

for ik=1:length(lags)
   ia=lags(ik);
   ncompb2=ncompb(ia:ia+na-1)-mean(ncompb(ia:ia+na-1));
   bb=sqrt(ncompb2*ncompb2');
   ccf(ik)=ncompa*ncompb2'/(aa*bb);
end
iz=~isfinite(ccf);
ccf(iz)=0;

return
