function cc = crosscortwin(seis,b,twin)
% cross-correlate vector a with b, with a shift of 0:n points
ns=size(seis,1);

cc=zeros(ns,ns);

 % Mask (Set width of window around 0 time to search for maximum.)
itw=round(twin./(2*ndt));
% Allow for itw longer than length of correlation.
itw(itw>nt/2)=nt/2;


[pair(:,1),pair(:,2)]=ind2sub(ns,find(triu(ones(ns,ns),1)));

l=length(a);
b1=[b,b]; % wrap b around instead of padding with zeros

for i = 0:n
  cc(i+1) = a*b1(1+i:l+i)'; % cc=a dot b, with b shifted by i
end

return



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