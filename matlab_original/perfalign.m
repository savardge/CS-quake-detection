function [tdel,sigr,cc]= perfalign( seis, twin )
%PERFALIGN: Apply mccc recursively until there is no further improvement in alignment

% global weight ndt
weight=0;
ndt=0.025;

if isstruct(twin)==1 %then assymmetric windows
    twinmin = twin.minim ;
    twinmax = twin.maxim ;
end

% initialize
prevmax=-5;
curmax=0;
tdel=zeros(1,size(seis,1));


while (curmax-prevmax)>1.0 % while we can improve the fit (0-lag autocorr.)
    
    prevmax=curmax;
    
    if isstruct(twin)==1 %then assymmetric windows
        [curtdel,~,sigr,cc]=mccc_target(seis,twinmin,twinmax,weight);
    else
        [curtdel,~,sigr,cc]=mccc2(seis,twin,weight);
    end
    
    % shift
    seis=shift(seis,ndt,curtdel);
    % stack and 0-lag autocorrelation
    
    curmax=sum(sum(seis,1),2);
    
%     hist=[hist; curmax];
%     curmax-prevmax
    
    % update delay times
    tdel=tdel+curtdel;
    
end

% figure;plot(hist)

end

