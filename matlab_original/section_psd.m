function [h]=section_psd(seis,dt,stalst)

% SECTION(SEIS,DT,STALST) plots a Power spectral density section of 
% seismograms in 2D array SEIS. DT is sample interval 
ns=size(seis,1);
nfft=2^nextpow2(size(seis,2));
[p,f] = periodogram(seis(1,:),hamming(size(seis,2)),nfft,1/dt,'power');
p = zeros(ns,length(f)); 
aflag=1;
idbadf=zeros(1,ns);
areapk=zeros(1,ns);
for is=1:ns
    [p(is,:),f] = periodogram(seis(is,:),hamming(size(seis,2)),nfft,1/dt,'power');
    S=p(is,:);
    [m I]=max(S); % first peak
   if (I-5<1) % to avoid negative subscript indices
       areapk(is) =(sum(S(1:I+5))/sum(S)*100);
   elseif (I+5>length(S))
       areapk(is) =(sum(S(I-5:length(S)))/sum(S)*100);
   else
       areapk(is) =(sum(S(I-5:I+5))/sum(S)*100);
   end
   
   % Apply the thresholds
   fid=fopen('FreqThresh.dat','a');
   fprintf(fid,'%d\t%5.2f\t%5.2f\t%5.2f\n',1,(sum(S(f<2))/sum(S)*100),(sum(S(f>5))/sum(S)*100),areapk(is));
   fclose(fid);
   
    if (sum(S(f<2))/sum(S)*100)>70 || (sum(S(f>5))/sum(S)*100)>70 || areapk(is)>70
        idbadf(is)=1; % then we remove this station
    end
end


ix=find(f>10,1,'first');
f=f(1:ix);
p=p(:,1:ix);

ny=size(p,1);

if nargin < 4
  aflag=-1;
end

ic=fix(rand(1)*7)+1;
col=['y','m','c','r','g','b','k'];

yaxe=[1:ny];
if aflag < 0
  for iy=1:ny
    xymat(iy,:)=iy-0.7*p(iy,:)/max(abs(p(iy,:))+0.0000001);
  end
elseif aflag == 0
  for iy=1:ny
    xymat(iy,:)=iy-8.0*p(iy,:)/max(max(abs(p))+0.0000001);
  end
else
  for iy=1:ny
    xymat(iy,:)=iy-1.0*p(iy,:)/aflag;
  end
end

h=plot(f,xymat,['k-']);
axis('ij');
xlab=xlabel('Frequency [Hz]');
ylab=ylabel([' ']);
ax1=gca;
if ~isempty(stalst)
    set(ax1,'YTickLabel',stalst,'YTick',[1:ny],'FontSize',14)
end
ylim([0,ny+1]);
set(gca,'FontName','Helvetica','FontSize',14,'Clipping','off','layer','top');
set(xlab,'FontName','Helvetica','FontSize',14);
set(ylab,'FontName','Helvetica','FontSize',14);

% figure(31);
% loglog(f,sum(p(1:end-1,:),1));title('stack of PSD')
% xlim([1 8])

end