for is=1:ns
    figure(1);[coefs,sgram]=cwt(ec(is,:),1:16,'mexh','scal');
    figure(2);plot(sum(sgram,1))
    pause
end

%% Morlet Fourier factor for future conversion scale to frequency

% w0=5;
% MorletFourierFactor = 4*pi/(w0+sqrt(2+w0^2));
% freq = 1./(scales.*MorletFourierFactor);
    