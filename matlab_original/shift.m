function [sfn] = shift(fn,dt,t0)
% FUNCTION [SFN] = SHIFT(FN,DT,T0)
% Function SHIFT produces shifted time series SFN from
% input time series FN where T0 is the desired time
% shift, DT is the sampling interval. If array of 
% traces is passed then T0 should be an array too.
% Uses FFT (ie sinc(x)) interpolation.

% Preliminaries.
i=sqrt(-1);
n1=max(size(fn));
m1=min(size(fn));
omega=(0:n1-1)*2*pi/(n1*dt);

% Loop through traces. Check for even/odd length since it affects
% definition of Nyquist frequency.
if(mod(n1,2)==0)
  for i1=1:m1
    ffn=fft(fn(i1,:),n1);
    ffn(1:n1/2)=ffn(1:n1/2).*exp(-i*omega(1:n1/2)*t0(i1));
    ffn(n1/2+1)=ffn(n1/2+1)*cos(pi*t0(i1)/dt);

% Negative frequencies (faster just to use complex conjugate relation for
% real signal).
%  ffn(n1/2+2:n1)=ffn(n1/2+2:n1).*exp(-i*(omega(n1/2+2:n1)*t0(i1)-2*pi*t0(i1)/dt));

    ffn(n1/2+2:n1)=conj(ffn(n1/2:-1:2));
    dum=real(ifft(ffn,n1));
    sfn(i1,:)=dum(1:n1);
  end
else
  for i1=1:m1
    ffn=fft(fn(i1,:),n1);
    ffn(1:(n1+1)/2)=ffn(1:(n1+1)/2).*exp(-i*omega(1:(n1+1)/2)*t0(i1));
    ffn((n1+1)/2+1:n1)=conj(ffn((n1+1)/2:-1:2));
    dum=real(ifft(ffn,n1));
    sfn(i1,:)=dum(1:n1);
  end
end
return
