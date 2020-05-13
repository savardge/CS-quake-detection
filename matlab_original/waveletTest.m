load wf_048.mat

%% Get channels from uk=[ncomp;ecomp;zcomp]

ta=0;
tb=60;
ns=size(uk,1)/3;

ncomp=uk(1:ns,:);
ecomp=uk(ns+1:2*ns,:);
zcomp=uk(2*ns+1:end,:);

L=size(uk,2);

t=0:ndt:ndt*(L-1);

%% Display
figure;
subplot(1,3,1);
section(ncomp,0,ndt);
ylabel([' ']);
% ylim([0,ns+1])
% xlim([ta,tb]);
ax1=gca;
set(ax1,'YTickLabel',ordlst,'YTick',[1:ns],'FontSize',14)
subplot(1,3,2);
section(ecomp,0,ndt);
ylabel([' ']);
% ylim([0,ns+1])
% xlim([ta,tb]);
ax1=gca;
set(ax1,'YTickLabel',ordlst,'YTick',[1:ns],'FontSize',14)
subplot(1,3,3);
section(zcomp,0,ndt);
ax1=gca;
set(ax1,'YTickLabel',ordlst,'YTick',[1:ns],'FontSize',14)
ylabel([' ']);

%% Cross-correlate


% TIME THE PROCESS
starttime = now;
tic;
t0 = cputime;

% llag=800;
% lwindow=120;

ind=1441;
% trace1=ecomp(9,ind:end);
% trace2=ecomp(10,ind:end);
% trace3=ecomp(11,ind:end);
% trace4=ecomp(12,ind:end);
% trace5=ecomp(13,ind:end);
% trace6=ecomp(14,ind:end);

trace1=ncomp(9,ind:end);
trace2=ncomp(10,ind:end);
trace3=ncomp(11,ind:end);
trace4=ncomp(12,ind:end);
trace5=ncomp(13,ind:end);
trace6=ncomp(14,ind:end);

tt=t(ind:end);

figure;
subplot(6,1,1);plot(tt,trace1)
subplot(6,1,2);plot(tt,trace2)
subplot(6,1,3);plot(tt,trace3)
subplot(6,1,4);plot(tt,trace4)
subplot(6,1,5);plot(tt,trace5)
subplot(6,1,6);plot(tt,trace6)

len=60;
strt=461;
[ wndw ] = extractWindow( trace1, strt, len );
subplot(6,1,1);hold on; plot(tt(strt:strt+len-1),wndw,'r-')

MAXLAG=length(trace1)-1;

figure;

[C1, lags1]=CCGM(trace1,wndw,MAXLAG);
subplot(6,1,1);plot((lags1-strt).*ndt,C1);

[C2, lags2]=CCGM(trace2,wndw,MAXLAG);
subplot(6,1,2);plot((lags2-strt).*ndt,C2);

[C3, lags3]=CCGM(trace3,wndw,MAXLAG);
subplot(6,1,3);plot((lags3-strt).*ndt,C3);

[C4, lags4]=CCGM(trace4,wndw,MAXLAG);
subplot(6,1,4);plot((lags4-strt).*ndt,C4);

[C5, lags5]=CCGM(trace5,wndw,MAXLAG);
subplot(6,1,5);plot((lags5-strt).*ndt,C5);

[C6, lags6]=CCGM(trace6,wndw,MAXLAG);
subplot(6,1,6);plot((lags6-strt).*ndt,C6);

% DISPLAY RUN TIMES
tclock = toc;
tcpu = cputime-t0;
if tclock>20
    disp(['Clock time to complete correlations: ' num2str(tclock,'%5.1f') ' s']);
    disp(['CPU time to complete correlations:   ' num2str(tcpu,'%5.1f') ' s']);
end;


%% Fourier transform

sig=ncomp(11,:);
[f Ft]= fourier( sig, ndt );

plot(f,2*abs(Ft));
xlabel('Frequency (Hz)')
ylabel('Amplitude')

%% GISMO

station=ordlst(9,:);
channel='';
samplerate=1/ndt;
starttime=datenum([2005,3,13,10,0,0]);
data=ncomp(9,:);
w=waveform(station,channel,samplerate,starttime,data)

C=correlation(w)

%% Short windowed transform

% figure;
% specgm = WindowFT(sig);

%% Wavelet transform

% figure
% cwtnPGC=cwt(ncomp(9,:),1:32,'mexh','3Dabslvl')
% figure
% cwtnKELB=cwt(ncomp(10,:),1:32,'mexh','3Dabslvl')
% figure
% cwtnMGCB=cwt(ncomp(11,:),1:32,'mexh','3Dabslvl')