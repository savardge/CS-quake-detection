function flag=plotspec( seis,show )
%PLOTSPEC : plot power spectrum for each waveform in seis
T=45; % treshold for % of variance in band
Tflow=2;
Tfhigh=5;

global ndt

sz=size(seis,1);
figure(11); clf;
count=0;
for is=1:sz
    [S,f]=periodogram(seis(is,:),[],'onesided',size(seis,2),1/ndt);
    [m I]=max(S);
    
    if show==1
       
        subplot(sz,1,is)
        plot(f,S)
        hold on
        plot(f(I),m,'r*')
%         title(['% of variance in 1st peak: '
%         num2str((m*f(I))/sum(f.*S)*100)]) 
        title(['% <' num2str(Tflow) 'Hz: ' num2str(sum(S(f<Tflow).*f(f<Tflow))/sum(f.*S)*100) '    % >' num2str(Tfhigh) 'Hz: ' num2str((sum(S(f>Tfhigh).*f(f>Tfhigh))/sum(f.*S)*100)), '    % in 1st peak: ' num2str(sum(S(I-5:I+5))/sum(S)*100)])
        xlabel('Frequency (Hz)')
        ylabel('Power per freq.')
        xlim([0,10])
    end
    
    if (sum(S(f<Tflow).*f(f<Tflow))/sum(f.*S)*100)>T || (sum(S(f>Tfhigh).*f(f>Tfhigh))/sum(f.*S)*100)>T
        count=count+1;
    end
end
if count>2
    flag=1; % Low freq. noise: Expulse detection
else
    flag=0;
end
end

