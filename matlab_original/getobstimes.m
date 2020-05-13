clear all
format compact
addpath('/home/savardge/TOOLBOXES/grep/')
addpath('/home/savardge/TOOLBOXES/subtightplot/')
global ndt
ndt=0.025;
addpath('/home/savardge/hypfiles/VI/')
% cd('/home/savardge/hypfiles/NW/')
% cd('/home/savardge/hypfiles/old_Washington/NW_obs/')
events=dir('tt_04*.arc');
ne=size(events,1);
addpath('/home/savardge/')
masterstationlist='    ';
[b,a]=butter(4,[1.25 7]/(40/2),'bandpass');
stalst=['GOWB';'PFB ';'TWGB'; 'TWBB'; 'TSJB'; 'LZB '; 'TWKB'; 'MGCB'; 'KELB'; 'KLNB'; 'PGC ';'SHDB';'SSIB'; 'SILB'; 'SNB ';'KHVB';'SHVB';'VGZ ';'LOP3';'SOK4';'GLBC';'LCBC';'JRBC'];
% stalst=['GLBC';'JRBC';'KLNB';'LCBC';'LZB ';'MGCB';'PFB ';'PGC ';'SHDB';'SHVB';'SILB';'SNB ';'SSIB';'TSJB';'TWBB';'TWKB'; 'VGZ '];
% stalst=['OPC ';'SQM ';'BS11';'W020';'FACB';'DOSE';'W040';'PL11';'W060';'N050';'N060';'B001';'B007';'B013';'PA01';'DR01';'BH01';'CL02';'LC02';'TB01';'GC01';'OSD ';'STW ';'GNW '];
% ColorSet=varycolor(size(stalst,1));
% figure(4);set(gca,'ColorOrder',ColorSet);hold all

% for istat=4:size(stalst,1)
%     station=stalst(istat,:);
% ---------------------------------------
dir1='/home/savardge/hypfiles/VI/';
%     dir1='/home/savardge/hypfiles/NW/'
%         dir1='/home/savardge/hypfiles/old_Washington/NW_obs/'
dir2='/mnt/data4/data/bostock/CASC/Data/Stations/';
%     dir2='/mnt/data4/data/bostock/CAFE/Data/Stations/'
%         stalst=['GOWB';'PFB ';'TWGB'; 'TWBB'; 'TSJB'; 'LZB '; 'TWKB'; 'MGCB'; 'KELB'; 'PGC ';'SHDB';'SSIB'; 'SILB'; 'SNB ';'KHVB';'SHVB';'VGZ ';'LOP3';'SOK4';'GLBC';'LCBC';'JRBC'];
% stalst=['OPC ';'SQM ';'BS11';'W020';'FACB';'DOSE';'W040';'PL11';'B001';'B003';'B004';'B007';'B013';'PA01';'DR01';'BH01';'CL02';'LC02';'TB01';'GC01';'STW ';'GNW '];
ns=size(stalst,1);
[xs,ys]=read_crd(stalst);


for ie=1:ne
    
    
    event=events(ie).name(4:9)
    fname=['tt_',event,'.arc']
    
    
    [Aout,fb]=grep('-n ', 'GSC', fullfile(dir1,fname));
    
    freqStack  = zeros(257,1);
    % freqStack=struct('total',zeros(257,1),)
    countdetect=0;
    
    flag=0;
    is=0;
    fid=fopen(fullfile(dir1,fname),'r');
    for k=1:fb.mlines
        tline = fgetl(fid);
        
        if isempty(tline)
            flag=0;
            disp([num2str(is),' stations read'])
            
            masterstationlist = unique([masterstationlist; ordlst],'rows');
            % Plot waveforms
            %         figure(1);clf;
            %         section([ncomp;ecomp],0,ndt,-1,ordlst);
            if HorizErr<10
                countdetect=countdetect+1;
                
                %             % Exclude specific station
                %             stnbad='SSIB';
                %             ibad=find( strncmp(stnbad,cellstr(ordlst),4));
                %             ordlst(ibad,:)=[];
                %             ncomp(ibad,:)=[]; ecomp(ibad,:)=[]; zcomp(ibad,:)=[];
                %             is=is-1;
                
                % compute PSD and stack for N and E
                %                     for ia=1:is
                %
                %                         if strncmp(ordlst(ia,:),station,4)==1
                %                             %                 id=find( strncmp(ordlst(ia,:),cellstr(stalst),4));
                %
                %                             nfft=2^nextpow2(size(ncomp,2));
                %
                %                             % For N
                % %                             [p,f] = periodogram(ncomp(ia,:),hamming(size(ncomp,2)),nfft,1/ndt,'power');
                %                             [p,f] = periodogram(zscore(ncomp(ia,:)'),hamming(size(ncomp,2)),nfft,1/ndt,'power');
                %                             %                 figure(2); plot(f,p); xlim([1 8]); title('N')
                %                             if ~isreal(p)
                %                                 keyboard
                %                             end
                %
                %                             freqStack = freqStack + p;
                %
                %                             % For E
                %                             [p,f] = periodogram(zscore(ecomp(ia,:)'),hamming(size(ncomp,2)),nfft,1/ndt,'power');
                %                             if ~isreal(p)
                %                                 keyboard
                %                             end
                %                             freqStack = freqStack + p;
                %
                %                         end
                %                     end
                
                % Plot waveforms used in detection
%                 secplot(ncomp,ecomp,zcomp,1,ordlst,0,'N','E','Z', 3);
                
                %                 flagsave=input('Save ? 1 yes, 0 no ');
                %                 if flagsave==1
                %                     saveas(gcf,['/home/savardge/NW_detections/',event,'_',time,'__error',sprintf('%.0f',round(HorizErr)),'_filt.jpg'])
                %                 end
%                 secplot_psd(ncomp,ecomp,zcomp,10,ordlst,'N','E','Z', 3)
                
                
                ncomp_filt=zeros(size(ncomp)); ecomp_filt=zeros(size(ncomp)); zcomp_filt=zeros(size(ncomp));
                for is=1:size(ncomp,1)
                    ncomp_filt(is,:)=filtfilt(b,a,ncomp(is,:));
                    ecomp_filt(is,:)=filtfilt(b,a,ecomp(is,:));
                    zcomp_filt(is,:)=filtfilt(b,a,zcomp(is,:));
                end
                secplot(ncomp_filt,ecomp_filt,zcomp_filt,2,ordlst,0,'N','E','Z', 3);
                
                
                
                %                 [ PApprox, f, integ ] = regressSpectrum( [ncomp;ecomp] );
                %                 figure(20); imagesc(integ); axis xy ; colorbar
                
                secplot_st(ncomp_filt,ecomp_filt,zcomp_filt,20,ordlst);
                
                
            end
            
            
            %         % PLOT S-TRANSFORM---------------------------------------------
            %         figure(2);clf;
            %         nn=size(ordlst,1);
            %         dumlst=ordlst;
            %         k=1;
            %         for ii=1:nn
            %
            %             subtightplot(3,nn,k,0.03)
            %             [strans,t,f] = st(ncomp((ii),:),0,180,ndt,1);
            %             imagesc(t,f,abs(strans));axis xy;%colorbar;
            %             title(['Component N, Station: ',dumlst(ii,:)])
            %
            %             subtightplot(3,nn,k+nn,0.03)
            %             [strans,t,f] = st(ecomp((ii),:),0,180,ndt,1);
            %             imagesc(t,f,abs(strans));axis xy;%colorbar;
            %             title(['Component E, Station: ',dumlst(ii,:)])
            %
            %             subtightplot(3,nn,k+2*nn,0.03)
            %             [strans,t,f] = st(zcomp((ii),:),0,180,ndt,1);
            %             imagesc(t,f,abs(strans));axis xy;%colorbar;
            %             title(['Component Z, Station: ',dumlst(ii,:)])
            %
            %             k=k+1;
            %         end
            % ----------------------------------------------------
            
            % Plot for all stations according to predicted arrival-----------
            [xhyp,yhyp,~] = deg2utm(lat,-lon);
            xhyp=xhyp/1000;yhyp=yhyp/1000;
            zhyp=depth;
            xep=sqrt((xs-xhyp).^2+(ys-yhyp).^2); % distance station to template hypocenter horizontaly
            phi=atan2(ys-yhyp,xs-xhyp); % azimuth station to template
            [tpp,tss,~,~,~,~,pr]=ttimes_hess(zhyp,xep,1.73,phi);
            twinp=[tpp'-2.0;tpp'+3.0]+Otime;
            twins=[tss'-5.0;tss'+5.0]+Otime;
            iwinp=round(twinp/ndt);
            iwins=round(twins/ndt);
            nl=numel(iwins(1,1):iwins(end,1));
            ncall=zeros(size(stalst,1),nl);ecall=ncall;zcall=ncall;
            for ik=1:size(stalst,1)
                [ncall(ik,:), ecall(ik,:),zcall(ik,:)]=loc2win2(event,iwins(1,ik),iwins(end,ik),stalst(ik,:),dir2);
            end
            
            
            % Keep waveforms not empty
            idkeep=find((abs(sum(ncall,2))>0)+(abs(sum(ecall,2))>0));
            nc=ncall(idkeep,:);
            ec=ecall(idkeep,:);
            zc=zcall(idkeep,:);
            stalstk=stalst(idkeep,:);
            nsk=size(stalstk,1);
            nc_filt=zeros(size(nc)); ec_filt=zeros(size(nc)); zc_filt=zeros(size(nc));
            for is=1:nsk
                    nc_filt(is,:)=filtfilt(b,a,nc(is,:));
                    ec_filt(is,:)=filtfilt(b,a,ec(is,:));
                    zc_filt(is,:)=filtfilt(b,a,zc(is,:));
            end
            secplot(nc_filt,ec_filt,zc_filt,3,stalstk,0,'N','E','Z', 3);
            [stnproj,steproj,stzproj] = secplot_st(nc_filt,ec_filt,zc_filt,4,stalstk);
            
            % -------------------------------------------------------
            
            % Attempt alignment --------------------------------------
            uk=[nc;ec];
            [detect,structD]=secondpassdetect(uk,stalstk,4,Otime+tss');
            
            % -----------------------------------------------------------
            
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
            

            clear ordlst ncomp ecomp zcomp %ncall ecall zcall
            is=0;
            disp('---------------------------------------------------------')
            
        elseif ismember(k,fb.line) && strncmp(tline,'20',2) % Location header
            %                 disp('Header')
            flag=1;
            hr=str2double(tline(9:10));
            mn=str2double(tline(11:12));
            hrmin=hr*3600+mn*60;
            sec=str2double(tline(13:16))/100;
            Otime=hr*3600+mn*60+sec;
            sprintf('%s @ %d:%d:%.2f',tline(1:8),hr,mn,sec)
            time=sprintf('hr%d_min%d_sec%.2f',hr,mn,sec);
            lat=str2double(regexprep(tline(17:18), ' ', '0'))+(str2double(regexprep(tline(20:21), ' ', '0'))+str2double(regexprep(tline(22:23), ' ', '0'))/100)/60; % verify min/sec
            lon=str2double(regexprep(tline(24:26), ' ', '0'))++(str2double(regexprep(tline(28:29), ' ', '0'))+str2double(regexprep(tline(30:31), ' ', '0'))/100)/60; % verify min/sec
            depth=str2double(tline(32:36))/100;
            HorizErr=str2double(tline(86:89))/100;
            sprintf('lat: %.2f, lon: %.2f, depth: %.2f.  Horiz. error: %.2f',lat,lon,depth,HorizErr)
            
            
        elseif flag==1 && isletter(tline(1))  % Station line
            if (strcmp(tline(14:16),'IPU'))
                %                     disp('P time')
                Ptime=str2double(tline(30:34))/100;
            else
                is=is+1;
                %                     disp('S time')
                sprintf('%s',tline)
                ordlst(is,:)=tline(1:4);
                Stime(is)=str2double(tline(42:46))/100+hrmin;
                Epidistance(is)=str2double(tline(75:78));
                Azimuth(is)=str2double(tline(92:94));
                EmergenceAngle(is)=str2double(tline(79:81));
                sprintf('Epidistance: %d, Azimuth: %d, EmergenceAngle: %d',Epidistance(is),Azimuth(is),EmergenceAngle(is))
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
            
        end
    end
    
    % Save the PSD stack
    %         if exist('f')==1
    %             freqStack = freqStack./countdetect; % Normalize
    %             loglog(f,freqStack,'DisplayName',[event,'_',station]);xlim([1 8 ]); title('Frequency stack')
    %
    %             save(['PSDstack_',event,'_',station,'_norm.mat'],'f','freqStack','countdetect')
    %         end
    % -------------------------------------------
    
    %         clearvars -except ie ne events masterstationlist freqStack f station stalst ColorSet
    
    %     end
    
end