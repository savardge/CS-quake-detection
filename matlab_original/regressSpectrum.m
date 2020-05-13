function [ PApprox, f, integ ] = regressSpectrum( uk )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
addpath('/home/savardge/TOOLBOXES/logfit/')
ns=size(uk,1); nn=size(uk,2); nfft=2^nextpow2(nn);
stack=[];
integ=zeros(ns,ns);
for is=1:ns-1
    for js=is+1:ns

        seis1 = uk(is,:);
        seis2 = uk(js,:);

        [Cxy, F ] = mscohere( seis1, seis2, hanning(nfft/2) , fix(0.5*nn) , nfft, 40 );
        if isempty(stack)
            stack=zeros(size(Cxy));
        end
        [Pxx] = pwelch( seis1, hanning(nfft/2) , fix(0.5*nn) , F, 40 );
        [Pyy] = pwelch( seis2, hanning(nfft/2) , fix(0.5*nn) , F, 40 );
        Pave = (Pxx./max(Pxx) + Pyy./max(Pyy))/2;
        stack=stack+Cxy.*Pave;
        integ(is,js)=trapz(F,Cxy.*Pave);
    end
end
ind=find(F>1,1,'first'):find(F<8,1,'last');
f=F(ind);
[slope, intercept] = logfit(f,stack(ind),'logy');
PApprox=(10^intercept)*(10^slope).^f;

ind2=find(f>1,1,'first'):find(f<3,1,'last');
ind3=find(f>3,1,'first'):find(f<5,1,'last'); 
if length(ind3)>length(ind2)
    ind3=ind3(1:length(ind2));
elseif length(ind3)<length(ind2)
    ind3=[ind3 +ind3(end)+(1:(length(ind2)-length(ind3)))];
end
PApprox(ind2)=flipud(PApprox(ind3));

end

