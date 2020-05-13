function [ix,tdel,rmean,sigr,cc,cmax] = ccsearch(seis,twin,smax)
% funtion [ix,tdel,rmean,sigr,cc,cmax] = ccsearch(seis,ndt,twin,smax)
% FUNCTION [IX,TDEL,RMEAN,SIGR,CC,CMAX] = CCSEARCH(SEIS,NDT,TWIN,SMAX)
% OUTPUTS - 
% IX : indices of retained stations.
% TDEL : time delays for retained stations.
% RMEAN : mean corrrelation coefficients for retained stations.
% SIGR : rms timing residuals for retained stations.
% CC : cross correlation coefficient matrix.
% INPUTS - 
% SEIS : seismogram section
% NDT : sample interval
% TWIN : window about zero to search in correlation window.
% SMAX : maximum allowable residual

%  twin
% global weight 
weight=0;

dc=0.01;
pflag=1;
ix=[];
tdel=[];

% figure(10);section(seis,0,ndt);
% pause;

% Compute CC matrix once.
[~,rmean,sigr,cc]=mccc2(seis,twin,weight);
% [~,rmean,sigr,cc]=mccc22(seis,twin);
cmax=max(max(cc));

% Loop through, gradually decreasing cross correlation 
% threshold to allow more stations until we start to see inconsistent times.
while pflag && cmax > 0
    [ir,ic]=find(cc>cmax);
    ix2=unique([ir;ic]);
    if ~isempty(ix2) && length(ix2) > 2
%         col=[cc ; cc']; col=reshape(col(col>0),ns-1,ns);
        if weight==0
            [tdel2,rmean2,sigr2]=mccc2(seis(ix2,:),twin,weight);
%             [tdel2,rmean2,sigr2]=mccc22(seis(ix2,:),twin);
            if max(sigr2) > smax  || (size(ix2,1) == size(seis,1))
                pflag=0;
            else
                ix=ix2;
                tdel=tdel2;
                rmean=rmean2;
                sigr=sigr2;
            end
        elseif weight==1
            [tdel2,rmean2,sigr2,cc2,tcc,tdel2w,res2w,sigr2w]=mccc2(seis(ix2,:),twin,weight);
            if max(sigr2w) > smax   || (size(ix2,1) == size(seis,1))
                pflag=0;
            else
                ix=ix2;
                tdel=tdel2w;
                rmean=rmean2;
                sigr=sigr2w;
            end
        end
        
    end
    cmax=cmax-dc;
end

return
