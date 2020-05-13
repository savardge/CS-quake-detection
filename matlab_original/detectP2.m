clear all
format compact
addpath('/home/savardge/TOOLBOXES/grep/')
addpath('/home/savardge/TOOLBOXES/subtightplot/')
addpath('/home/savardge/TOOLBOXES/Stransform/')
global ndt
ndt=0.025;
cd('/home/savardge/hypfiles/VI/')
events=dir('tt_05*.arc');
ne=size(events,1);
cd('/home/savardge/')
masterstationlist='    ';
stalst=['GOWB';'PFB ';'TWGB'; 'TWBB'; 'TSJB'; 'LZB '; 'TWKB'; 'MGCB'; 'KELB'; 'KLNB'; 'PGC ';'SHDB';'SSIB'; 'SILB'; 'SNB ';'KHVB';'SHVB';'VGZ ';'LOP3';'SOK4';'GLBC';'LCBC';'JRBC'];
% stalst=['GOWB';'PFB '; 'TWBB'; 'TSJB'; 'LZB '; 'TWKB'; 'MGCB';  'KLNB'; 'PGC ';'SHDB';'SSIB'; 'SILB'; 'SNB ';'VGZ ';'GLBC';];
%         stalst=['GOWB';'PFB ';'TWGB'; 'TWBB'; 'TSJB'; 'LZB '; 'TWKB'; 'MGCB'; 'KELB'; 'PGC ';'SHDB';'SSIB'; 'SILB'; 'SNB ';'KHVB';'SHVB';'VGZ ';'LOP3';'SOK4';'GLBC';'LCBC';'JRBC'];
% stalst=['GLBC';'JRBC';'KLNB';'LCBC';'LZB ';'MGCB';'PFB ';'PGC ';'SHDB';'SHVB';'SILB';'SNB ';'SSIB';'TSJB';'TWBB';'TWKB'; 'VGZ '];
ns=size(stalst,1);
ColorSet=varycolor(size(stalst,1));
figure(4);set(gca,'ColorOrder',ColorSet);hold all

for ie=1:ne
    
    event=events(ie).name(4:9)
    fname=['tt_',event,'.arc'];
    dir1='/home/savardge/hypfiles/VI/';
    % dir1='/home/savardge/hypfiles/SW/'
    dir2='/mnt/data4/data/bostock/CASC/Data/Stations/';
    % dir2='/mnt/data4/data/bostock/CAFE/Data/Stations/'
    
    [xs,ys]=read_crd(stalst);
    
    ndt=0.025;
    [Aout,fb]=grep('-n ', 'GSC', fullfile(dir1,fname));
    flag=0;
    is=0;
    fid=fopen(fullfile(dir1,fname),'r');
    % Open Hypoinverse pick file.
    fidl=fopen(['tt',event,'_P'],'w');
    
    for k=1:fb.mlines % For each detection do
        tline = fgetl(fid);
        
        
        if ismember(k,fb.line) && strncmp(tline,'200',3)
            %%%%%%%%%%%%%%%%%% LOCATION HEADER %%%%%%%%%%%%%%%%%%
            flag=1;
            hr=str2double(tline(9:10));
            mn=str2double(tline(11:12));
            hrmin=hr*3600+mn*60;
            sec=str2double(tline(13:16))/100;
            Otime=hr*3600+mn*60+sec;
            sprintf('%s @ %d:%d:%.2f',tline(1:8),hr,mn,sec)
            lat=str2double(regexprep(tline(17:18), ' ', '0'))+(str2double(regexprep(tline(20:21), ' ', '0'))+str2double(regexprep(tline(22:23), ' ', '0'))/100)/60; % verify min/sec
            lon=str2double(regexprep(tline(24:26), ' ', '0'))++(str2double(regexprep(tline(28:29), ' ', '0'))+str2double(regexprep(tline(30:31), ' ', '0'))/100)/60; % verify min/sec
            depth=str2double(tline(32:36))/100;
            HorizErr=str2double(tline(86:89))/100;
            sprintf('lat: %.2f, lon: %.2f, depth: %.2f.  Horiz. error: %.2f',lat,lon,depth,HorizErr)
            
            
        elseif flag==1 && ~isempty(tline) && isletter(tline(1))  % Station line
            %%%%%%%%%%%%%%%%%% STATION LINES %%%%%%%%%%%%%%%%%%
            if (strcmp(tline(14:16),'IPU'))
                Ptime=str2double(tline(30:34))/100;
            else
                is=is+1;
                sprintf('%s',tline)
                ordlst(is,:)=tline(1:4);
                Stime(is)=str2double(tline(42:46))/100+hrmin;
                Epidistance(is)=str2double(tline(75:78));
                Azimuth(is)=str2double(tline(92:94));
                EmergenceAngle(is)=str2double(tline(79:81));
                %                     sprintf('Epidistance: %d, Azimuth: %d, EmergenceAngle: %d',Epidistance(is),Azimuth(is),EmergenceAngle(is))
                istrt=fix((Stime(is)-6)/ndt)+1;
                iend=fix((Stime(is)+4)/ndt)+1;
                if istrt > (24*3600)/ndt % origin time the day before
                    istrt=1; iend=401;
                end
                if iend-istrt < 400
                    iend=iend+(400-(iend-istrt));
                elseif iend-istrt > 400
                    iend=iend-((iend-istrt)-400);
                end
                [ncomp(is,:), ecomp(is,:),zcomp(is,:)]=loc2win2(event,istrt,iend,ordlst(is,:),dir2);
            end
            
            
        elseif isempty(tline)
            disp('All info read in')
            %%%%%%%%%%%%%%%%%% PLOT WAVEFORMS %%%%%%%%%%%%%%%%%%
            flag=0;
            disp([num2str(is),' stations read'])
            
            masterstationlist = unique([masterstationlist; ordlst],'rows');
            % Plot waveforms
            %         figure(1);clf;
            %         section([ncomp;ecomp],0,ndt,-1,ordlst);
            if HorizErr<3
                
                % Plot waveforms used in detection
                secplot(ncomp,ecomp,zcomp,1,ordlst,0,'N','E','Z', 3);
                %                 secplot_psd(ncomp,ecomp,zcomp,10,ordlst,'N','E','Z', 3)
                
                
                % Plot for all stations according to predicted arrival-----------
                disp('Plotting all stations for aligned arrivals')
                [xhyp,yhyp,~] = deg2utm(lat,-lon);
                xhyp=xhyp/1000;yhyp=yhyp/1000;
                zhyp=depth;
                xep=sqrt((xs-xhyp).^2+(ys-yhyp).^2); % distance station to template hypocenter horizontaly
                phi=atan2(ys-yhyp,xs-xhyp); % azimuth station to template
                [tpp,tss,~,~,~,~,pr]=ttimes_hess(zhyp,xep,1.73,phi);
                twins=[tss'-1.0;tss'+3.5];
                twinp=[tpp'-4.0;tss'];
                iwinp=round(twinp/0.025);
                iwins=round(twins/0.025);
                npad=max([iwins(2,:)-iwins(1,:),iwinp(2,:)-iwinp(1,:)])+fix(1/ndt);
                
                % Read detections into time series.  NDCOMPP and NDCOMPS are the P and S wave windows, respectively.
                [ndcomps,edcomps,zdcomps]=read_data(stalst,event,Otime,iwins,npad);
                [ndcompp,edcompp,zdcompp]=read_data(stalst,event,Otime,iwinp,npad);
                
                % Filtering
                [b,a]=butter(4,[1.5 7]/(40/2),'bandpass');
                ndcomps=zscore(ndcomps')'; edcomps=zscore(edcomps')'; zdcompp=zscore(zdcompp')';
                for is=1:size(ndcomps,1)
                    ndcomps(is,:)=filtfilt(b,a,ndcomps(is,:));
                    edcomps(is,:)=filtfilt(b,a,edcomps(is,:));
                    zdcompp(is,:)=filtfilt(b,a,zdcompp(is,:));
                end
                
                % Combine N and E S-waves in one array.
                seis=[ncomp;ecomp];
                cwinx=ndt; % CC window for mccc_coeff, not really important since best waveforms already aligned.
                
                % Align
                [tdel,rmean,sigr,cc] = mccc2(seis,1,0);
                nx=size(seis,1);
                seis=shift(seis,ndt,tdel);
                [lia,locb]=ismember(cellstr(ordlst),cellstr(stalst));
                iwinx=[iwins(:,locb),iwins(:,locb)];
            
                % NLS and NLP contain the length of the extraction windows. These are used in mccc_coeff to
                % measure the cross correlation coefficient between the P trace and the S traces.
                nls=iwinx(2,:)-iwinx(1,:)+1;
                nlp=iwinp(2,:)-iwinp(1,:)+1;
                
                
                % Cycle through P-waves-----------------------------------
                
                k=0; % Counter for number of P detections.
                idp=zeros(1,ns);
                tdp=zeros(1,ns);
                % Initiatilize polarization flag and window lengths.
                ipol=zeros(ns,1);
                nt=length(zdcompp);
               
                
                for is=1:ns
                    
                    nlw=[nls,nlp(is)]';
                    % Test first for positive polarities
                    seispU=[seis(ix,:);zdcompp(is,:)];
                    %                 seispU=zscore([seis(ix,:);zdcompp(is,:)]')';
                    
                    % This version of mccc deals with the final (P) trace using a 1-side correlation
                    % with the shorter S window scanning along the longer P-window (as above this
                    % typically leads to approximately 5 s worth of correlation coefficient)..
                    [tdelpU,rmeanpU,sigrpU,ccpU,ccmaxpU]=mccc_coeff(seispU,ndt,cwinx,nlw,0);
                    
                    flagup=0;
                    
                    %                 if (max(sigrpU) <= smax) && rmeanpU(end)>ccmin %(mean(ccpU(1:end-1,end))>0.4)
                    if (max(sigrpU) <= smax) && max(ccpU(:,end))>ccmin %(mean(ccpU(1:end-1,end))>0.4)
                        sprintf('sigr= %5.3f, cc= %5.3f\n',max(sigrpU),max(ccpU(:,end)))
                        %                     sigrp
                        %                     ccp(1:end-1,end)
                        %                     pause;
                        
                        k=k+1;
                        idp(k)=is;
                        ipol(k)=1;
                        scp(k)='U';
                        % tdp yields the full delay between the candidate P-wave and
                        % the remain S-waves (that have already been aligned and so
                        % should have almost the same delay if max(sigrp) is small).
                        tdp(k)=tdelpU(nx+1)-mean(tdelpU(1:nx));
                        flagup=1;
                        disp('Positive Polarity')
                        stalst(is,:)
                        seispU=shift(seispU,ndt,tdelpU);
                        
                        % Plot -------------
                        %                     h=figure(2);
                        %                     subplot(1,2,1)
                        %                     section(seispU,0,ndt,-1,[[stalst(istalstb(ix),:), repmat('.H',length(ix),1)];[stalst(is,:),'.V']]);
                        %                     title(['Time: ',num2str(tt(ie))]);
                        %                     subplot(1,2,2);
                        %                     section_psd(seispU,ndt,[[stalst(istalstb(ix),:), repmat('.H',length(ix),1)];[stalst(is,:),'.V']]);
                        %                     saveas(h,['PDetection_',date,'_t',num2str(fix(tt(ie))),'_sigr',num2str(fix(max(sigrpU)*1e3)),'ms_cc',num2str(fix(rmeanpU(end)*1e2)),'.jpg'])
                    end
                    
                    % Test negative polarities by simply negating same P channel.
                    %                 seispD=zscore([seis(ix,:);-zdcompp(is,:)]')';
                    seispD=[seis(ix,:);-zdcompp(is,:)];
                    [tdelpD,rmeanpD,sigrpD,ccpD,ccmaxpD]=mccc_coeff(seispD,ndt,cwinx,nlw,0);
                    
                    % Keep only if higher correlation than positive polarity
                    %                 if max(sigrpD) <= smax && rmeanpD(end)>ccmin && (rmeanpD(end)>rmeanpU(end)) % (mean(ccpD(1:end-1,end)>0.4))
                    if max(sigrpD) <= smax && max(ccpD(:,end))>ccmin && (rmeanpD(end)>rmeanpU(end)) % (mean(ccpD(1:end-1,end)>0.4))
                        sprintf('sigr= %5.3f, cc= %5.3f\n',max(sigrpD),max(ccpD(:,end)))
                        if flagup==0
                            k=k+1;
                        end
                        idp(k)=is;
                        ipol(k)=-1;
                        scp(k)='D';
                        tdp(k)=tdelpD(nx+1)-mean(tdelpD(1:nx));
                        disp('Negative Polarity')
                        stalst(is,:);
                        seispD=shift(seispD,ndt,tdelpD);
                        
                        % Plot -------------
                        %                     h=figure(2);
                        %                     subplot(1,2,1);
                        %                     section(seispD,0,ndt,-1,[[stalst(istalstb(ix),:), repmat('.H',length(ix),1)];[stalst(is,:),'.V']]);
                        %                     title(['Time: ',num2str(tt(ie))]);
                        %                     subplot(1,2,2);
                        %                     section_psd(seispD,ndt,[[stalst(istalstb(ix),:), repmat('.H',length(ix),1)];[stalst(is,:),'.V']]); title('Power spectral density')
                        %                     saveas(h,['PDetection_',date,'_t',num2str(fix(tt(ie))),'_sigr',num2str(fix(max(sigrpD)*1e3)),'ms_cc',num2str(fix(rmeanpD(end)*1e2)),'.jpg'])
                    end
                end
                
                % Log any detections. --------------------------------------
                if k > 0
                    
                    % Compute different temporal quantities for registering times in hyp format.
                    hr=fix(tt(ie)/3600);
                    mn=fix((tt(ie)-hr*3600)/60);
                    ssc=tt(ie)-hr*3600-mn*60;
                    dum=num2str(100+hr);
                    shr=dum(2:3);
                    dum=num2str(100+mn);
                    smn=dum(2:3);
                    time=['20',events(ie,1:6),shr,smn];
                    
                    % Print eventline data.
                    fprintf(fidl,'%s%5.2f                                                                                                 \n', time,ssc);
                    
                    % Determine absolute time from maximum of stack of S-envelopes (this should be consistent
                    % with Genevieve's earlier result).
                    twins2=[twins(1,:),twins(1,:)];
                    [dum,imx]=max(sum(abs(hilbert(seis(ix,:).').')));
                    
                    % S-times are sum of origin time, times of stack maximum, ccsearch alignment and onset of S-window.
                    ts=(imx-1)*ndt+twins2(ix)-tdel;
                    
                    for ik=1:k
                        
                        % P-times are sum of origin time, time of stack maximum and relative delays between P and S, and onset of P-window.
                        tp(ik)=(imx-1)*ndt+twinp(1,idp(ik))-tdp(ik);
                        fprintf(fidl,'%-4s PO ZBHZ IP%s%d%s%5.2f                                                                             \n',stalst(idp(ik),:),scp(ik),0,time,ssc+tp(ik));
                    end
                    
                    % Write S-wave times to file.
                    for jx=1:nx
                        scs='N';
                        if ix(jx) > ns
                            scs='E';
                        end
                        fprintf(fidl,'%-4s PO ZBHZ     %s            %5.2f%sSU%d                                                              \n',stalstb(ix(jx),:),time,ssc+ts(jx),scs,0);
                    end
                    % Print new line to separate detection from next detection.
                    fprintf(fidl,' \n');
                    
                    
                end
                % ---------------------------------------------------------
                
                
                % -------------------------------------------------------
                
                % Rotation for detection's hypocenter-----------------------------
                %         svcomps=zeros(size(nc));
                %         shcomps=zeros(size(nc));
                %         svcompp=zeros(size(nc));
                %         shcompp=zeros(size(nc));
                %         pcomps=zeros(size(nc));
                %         for isk=1:nsk
                %             [ svcomps(isk,:), shcomps(isk,:), pcomps(isk,:), rcomps(isk,:), tcomps(isk,:) ] = NEZtoSVSHP( phi(idkeep(isk)), pr(idkeep(isk)), nc(isk,:), ec(isk,:), zc(isk,:) );
                %         end
                %         secplot(nc,ec,svcomps,5,stalstk,0,'N','E','SV', 3);
                %         secplot(svcomps,shcomps,pcomps,6,stalstk,0,'SV','SH','P', 3);
                %
                %         seis= [nc;ec];
                %         [U,S,V]=svd(zscore(seis')');
                %         figure(10);plot(diag(S,0));title('Eigenvalues');
                %         figure(11);section(V(:,1:10)',0,ndt,-1,'k');title('10 first principal components')
                % ----------------------------------------------------------------
                
            end
            
            clear ordlst ncomp ecomp zcomp %ncall ecall zcall
            is=0;
            disp('---------------------------------------------------------')
            
        end
    end
    
    
end