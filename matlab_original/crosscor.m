function cc = crosscor(a,b,n)
% cross-correlate vector a with b, with a shift of 0:n points

l=length(a);
b1=[b,b]; % wrap b around instead of padding with zeros
cc=zeros(n,1);
for i = 0:n
  cc(i+1) = a*b1(1+i:l+i)'; % cc=a dot b, with b shifted by i
end

return