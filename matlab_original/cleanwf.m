function [ newuk,newlst ] = cleanwf( stalst, uk )
%CLEANWF : remove stations with empty waveforms or low energy (noisy)



% PARAMETERS -----------------------
Fs=1/0.025;
T=70; %70 threshold for % energy in the low or high frequencies
Tfhigh=6; %6; %  high frequency to look over for threshold T
Tflow=2; %2 low frequency to look below for threshold T
Tfmax=70; %70; % threshold for max of energy in the first peak
Tlow=60;

% Extract N & E channels 
ns=size(uk,1)/2;
nc=uk(1:ns,:);
ec=uk(ns+1:2*ns,:);

% Remove stations with empty window NaN
idnan=all(isnan(nc),2);
idkeep=find(~idnan);
nc=nc(idkeep,:);
ec=ec(idkeep,:);
stalst=stalst(idkeep,:);


% Keep waveforms not empty
idkeep=find((abs(sum(nc,2))>0)+(abs(sum(ec,2))>0));
nc=nc(idkeep,:);
ec=ec(idkeep,:);
stalst=stalst(idkeep,:);
ns=size(stalst,1);

% If there is no station left, exit
if size(stalst,1)==0
    newuk=[];
    newlst=[];
    return
end

% Envelopes.
ne=abs(hilbert(nc.').');
ee=abs(hilbert(ec.').');

% Remove low SNR signals (60 as threshold)
idlown=((sum(ne,2)./max(ne,[],2))<Tlow);
idlowe=((sum(ee,2)./max(ee,[],2))<Tlow);
idkeep=find(~(idlown+idlowe));
nc=nc(idkeep,:);
ec=ec(idkeep,:);
stalst=stalst(idkeep,:);
ns=size(stalst,1);

idbadf=zeros(ns,1);
% Normalize displacements to unit maximum amplitude and check power
% spectrum
for is=1:ns
    % Detrend
    nc(is,:)=detrend(nc(is,:));
    ec(is,:)=detrend(ec(is,:));
    % remove mean
    nc(is,:)=nc(is,:)-mean(nc(is,:));
    ec(is,:)=ec(is,:)-mean(ec(is,:));
    % normalize
    amax=max(max(abs([nc(is,:);ec(is,:)])));
    nc(is,:)=nc(is,:)./amax;
    ec(is,:)=ec(is,:)./amax;
   
    
    % Check spectrum for N
    N = length(nc(is,:));
    xdft = fft(nc(is,:));
    xdft = xdft(1:N/2+1);
    S = (1/(Fs*N)).*abs(xdft).^2;
    S(2:end-1) = 2*S(2:end-1);
    
    f = 0:Fs/N:Fs/2;
%     [S,f]=periodogram(nc(is,:),[],'onesided',size(nc,2),Fs);

    
   [~, I]=max(S); % first peak
   if (I-5<1) % to avoid negative subscript indices
       areapk =(sum(S(1:I+5))/sum(S)*100);
   elseif (I+5>length(S))
       areapk =(sum(S(I-5:length(S)))/sum(S)*100);
   else
       areapk =(sum(S(I-5:I+5))/sum(S)*100);
   end
   
   % Apply the thresholds
%    disp(['% over ',num2str(Tfhigh),' Hz for station ',stalst(is,:),'.N : ',num2str((sum(S(f>Tfhigh).*f(f>Tfhigh))/sum(f.*S)*100))])
%    if (sum(S(f<Tflow).*f(f<Tflow))/sum(f.*S)*100)>T || (sum(S(f>Tfhigh).*f(f>Tfhigh))/sum(f.*S)*100)>T || areapk > Tfmax
%        figure(20);subplot(2,1,1); plot(nc(is,:))
%        subplot(2,1,2); plot(f,S); title([stalst(is,:),'.N'])
%        pause(1);
%    end
%    disp(['% under ',num2str(Tflow),' Hz for station ',stalst(is,:),'.N : ',num2str((sum(S(f<Tflow).*f(f<Tflow))/sum(f.*S)*100))])
%    (sum(S(f<Tflow))/sum(S)*100)>T
%    disp(['% in dominant peak for station ',stalst(is,:),'.N : ',num2str(areapk)])
%    areapk > Tfmax
   if (sum(S(f<Tflow).*f(f<Tflow))/sum(f.*S)*100)>T || (sum(S(f>Tfhigh).*f(f>Tfhigh))/sum(f.*S)*100)>T || areapk > Tfmax
        idbadf(is)=1; % then we remove this station
%         disp(['Remove station ',stalst(is,:),' .N'])
        
%         conf=input('confirm? (1=yes) ');
%         if conf==1
%             fid=fopen('FreqThresh.dat','a');
%             fprintf(fid,'%d\t%5.2f\t%5.2f\t%5.2f\n',0,(sum(S(f<2))/sum(S)*100),(sum(S(f>5))/sum(S)*100),areapk);
%             fclose(fid);
%         end
    end
    
    % Check spectrum for E
    xdft = fft(ec(is,:));
    xdft = xdft(1:N/2+1);
    S = (1/(Fs*N)).*abs(xdft).^2;
    S(2:end-1) = 2*S(2:end-1);
    f = 0:Fs/N:Fs/2;
%     [S,f]=periodogram(ec(is,:),[],'onesided',size(ec,2),1/ndt);

   [~, I]=max(S); % first peak
   if (I-5<1)
       areapk =(sum(S(1:I+5))/sum(S)*100);
   elseif (I+5>length(S))
       areapk =(sum(S(I-5:length(S)))/sum(S)*100);
   else
       areapk =(sum(S(I-5:I+5))/sum(S)*100);
   end
   
   % Apply the thresholds
%    disp(['% over ',num2str(Tfhigh),' Hz for station ',stalst(is,:),'.E : ',num2str((sum(S(f>Tfhigh).*f(f>Tfhigh))/sum(f.*S)*100))])
%    if (sum(S(f<Tflow).*f(f<Tflow))/sum(f.*S)*100)>T || (sum(S(f>Tfhigh).*f(f>Tfhigh))/sum(f.*S)*100)>T || areapk > Tfmax
%        figure(20);subplot(2,1,1); plot(ec(is,:))
%        subplot(2,1,2); plot(f,S); title([stalst(is,:),'.N'])
%        pause(1);
%    end
%    disp(['% under ',num2str(Tflow),' Hz for station ',stalst(is,:),'.E : ',num2str((sum(S(f<Tflow).*f(f<Tflow))/sum(f.*S)*100))])
%    (sum(S(f<Tflow).*f(f<Tflow))/sum(f.*S)*100)>T
%    disp(['% in dominant peak for station ',stalst(is,:),'.E : ',num2str(areapk)])
%    areapk > Tfmax
   if (sum(S(f<Tflow).*f(f<Tflow))/sum(f.*S)*100)>T || (sum(S(f>Tfhigh).*f(f>Tfhigh))/sum(f.*S)*100)>T || areapk > Tfmax
        idbadf(is)=1; % then we remove this station
%         disp(['Remove station ',stalst(is,:),' .E'])
       
%         conf=input('confirm? (1=yes) ');
%         if conf==1
%             fid=fopen('FreqThresh.dat','a');
%             fprintf(fid,'%d\t%5.2f\t%5.2f\t%5.2f\n',0,(sum(S(f<2))/sum(S)*100),(sum(S(f>5))/sum(S)*100),areapk);
%             fclose(fid);
%         end
    end
    
end

idkeep=find(~idbadf);
nc=nc(idkeep,:);
ec=ec(idkeep,:);
stalst=stalst(idkeep,:);
ns=size(stalst,1);

newuk=[nc;ec];
newlst=stalst;
if ns < 4
    newuk=[];
    newlst=[];
end

return

