function [ccf] = ccfnew(ncompa,ncompb)

% Simple function to compute running correlation coefficient
% of short window (S-LFE = ncompa) time series across long window 
% (P-LFE = ncompb) time series. Could be made more efficient.

ncompa=ncompa-mean(ncompa);
na=length(ncompa);
nb=length(ncompb);
aa=sqrt(ncompa*ncompa');
ccf=zeros(1,nb-na+1);

for ia=1:nb-na+1
   ncompb2=ncompb(ia:ia+na-1)-mean(ncompb(ia:ia+na-1));
   bb=sqrt(ncompb2*ncompb2');
   ccf(ia)=ncompa*ncompb2'/(aa*bb);
end
iz=~isfinite(ccf);
ccf(iz)=0;

return
