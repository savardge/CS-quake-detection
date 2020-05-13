function [ maxCC ] = updateMax( maxCC, new, strt )
% This function updates maximum CC coefficients and corresponding window
% start for each lag time.
val=squeeze(maxCC(1,:,1));
ind=new>val;
maxCC(1,ind,1)=new(ind);
maxCC(1,ind,2)=strt;

end

