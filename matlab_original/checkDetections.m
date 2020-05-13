clear all
format compact
addpath(genpath('/home/savardge/TOOLBOXES/'))
addpath('/home/savardge/')
global ndt
ndt=0.025;
% SeaJade
dir1='/fiore2_shared/shared/savardge/hypfiles/SeaJade/';
dir2='/fiore2_shared/shared/Data5/SEAJADE/Data/Stations/';
% SVI
% dir1='/fiore2_shared/shared/savardge/hypfiles/SVI/Past2005/';
% dir2='/fiore2_shared/shared/Data5/CASC/Data/Stations/';

cd(dir1)
% events=dir('tt_*_svi.arc');
events=dir('tt_*_sj.arc');
ne=size(events,1);
masterstationlist='    ';

cd('/home/savardge/CSdetection')
% SeaJade
stalst=[ 'S01 ';'S02 ';'S03 ';'S04 ';'S05 ';'S06 ';'S07 ';'S08 ';'S09 ';'S10 ';'S11 ';'S12 ';'S13 ';'S14 ';'S15 ';'S16 ';'S17 ';'S18 ';'S19 ';'S20 ';'S21 ';'S22 ';'S23 ';'S24 ';'S25 ';'S26 ';'S27 ';'S28 ';'S29 ';'S30 ';'S31 ';'S32 ';'S33A';'S34A';'S35A']

% SVI
% stalst=['GOWB';'PFB ';'TWGB'; 'TWBB'; 'TSJB'; 'LZB '; 'TWKB'; 'MGCB'; 'KELB'; 'PGC ';'SHDB';'SSIB'; 'SILB'; 'SNB ';'KHVB';'VGZ ';'LOP3';'SOK4';'GLBC';'LCBC';'JRBC'];

ns=size(stalst,1);
ColorSet=varycolor(size(stalst,1));
figure(4);set(gca,'ColorOrder',ColorSet);hold all

for ie=1:ne
    
    event=events(ie).name(4:9)
%     fname=['tt_',event,'_svi.arc'];
    fname=['tt_',event,'_sj.arc'];
   
    
    ndt=0.025;
    [Aout,fb]=grep('-n ', 'GSC', fullfile(dir1,fname));
    flag=0;
    is=0;
    fid=fopen(fullfile(dir1,fname),'r');
    
    for k=1:fb.mlines % For each detection do
        tline = fgetl(fid);
        
        
        if ismember(k,fb.line) && strncmp(tline,'20',2)
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
                ordlst1(is,:)=tline(1:4);
                Stime(is)=str2double(tline(42:46))/100+hrmin;
%                 Epidistance(is)=str2double(tline(75:78));
%                 Azimuth(is)=str2double(tline(92:94));
%                 EmergenceAngle(is)=str2double(tline(79:81));
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
                [ncomp(is,:), ecomp(is,:),zcomp(is,:)]=loc2win2(event,istrt,iend,ordlst1(is,:),dir2);
            end
            
            
        elseif isempty(tline)
            disp('All info read in')
            %%%%%%%%%%%%%%%%%% PLOT WAVEFORMS %%%%%%%%%%%%%%%%%%
            flag=0;
            disp([num2str(is),' stations read'])
            
            masterstationlist = unique([masterstationlist; ordlst1],'rows');
            % Plot waveforms
            %         figure(1);clf;
            %         section([ncomp;ecomp],0,ndt,-1,ordlst1);
            if HorizErr<3
                
                secplot_psd(ncomp,ecomp,zcomp,10,ordlst1,'N','E','Z', 3)
                
                % Plot waveforms used in detection
                [b,a]=butter(4,[1.5 6.5]/(40/2),'bandpass');
            
                for is=1:size(ncomp,1)
                    ncomp(is,:)=filtfilt(b,a,ncomp(is,:));
                    ecomp(is,:)=filtfilt(b,a,ecomp(is,:));
                    zcomp(is,:)=filtfilt(b,a,zcomp(is,:));
                end 
                secplot(ncomp,ecomp,zcomp,1,ordlst1,0,'N filtered','E filtered','Z filtered', 3);
                
                
                
                % PLOT S-TRANSFORM---------------------------------------------
                disp('Plotting the S transfroms')
                figure(2);clf;
                nn=size(ordlst1,1);
                dumlst=ordlst1;
                k=1;
                for ii=1:nn
                    
                    subtightplot(3,nn,k,0.03)
                    [strans,t,f] = st(ncomp((ii),:),0,180,ndt,1);
                    imagesc(t,f,abs(strans));axis xy;%colorbar;
                    title(['Component N, Station: ',dumlst(ii,:)])
                    
                    subtightplot(3,nn,k+nn,0.03)
                    [strans,t,f] = st(ecomp((ii),:),0,180,ndt,1);
                    imagesc(t,f,abs(strans));axis xy;%colorbar;
                    title(['Component E, Station: ',dumlst(ii,:)])
                    
                    subtightplot(3,nn,k+2*nn,0.03)
                    [strans,t,f] = st(zcomp((ii),:),0,180,ndt,1);
                    imagesc(t,f,abs(strans));axis xy;%colorbar;
                    title(['Component Z, Station: ',dumlst(ii,:)])
                    
                    k=k+1;
                end
                % ----------------------------------------------------
                
                % Map ------------------------------------------------
                figure(20); clf;
                minlat=48.18;
                maxlat=49.5;
                minlon=-125.5;
                maxlon=-123;
                m_proj('albers equal-area', 'lat',[minlat,maxlat],'long',[minlon,maxlon],'rect','on');
                m_usercoast('CASCADIA.map','patch',[0.7 0.7 0.7],'edgecolor','black')
                m_grid
                hold on
                [x,y]=m_ll2xy(-lon,lat,'clip','point');
                plot(x,y,'ko','MarkerSize',10,'MarkerFaceColor','r')
%                 load /home/savardge/networkCC/templates_hypo.mat
%                 [xt,yt]=m_ll2xy(-hlon,hlat,'clip','point');
%                 plot(xt,yt,'kd','MarkerSize',10)
                %---------------------------------------------------
                % Plot for all stations according to predicted arrival-----------
%                   [xs,ys]=read_crd_SJ(stalst);
                [xs,ys]=read_crd(stalst);
%                 disp('PLotting all stations for aligned arrivals')
                [xhyp,yhyp,~] = deg2utm(lat,-lon);
                xhyp=xhyp/1000;yhyp=yhyp/1000;
                zhyp=depth;
                xep=sqrt((xs-xhyp).^2+(ys-yhyp).^2); % distance station to template hypocenter horizontaly
                phi=atan2(ys-yhyp,xs-xhyp); % azimuth station to template
                [tpp,tss,~,~,~,~,pr]=ttimes_hess(zhyp,xep,1.73,phi);
%                 twinp=[tpp'-2.0;tpp'+3.0]+Otime;
%                 twins=[tss'-5.0;tss'+5.0]+Otime;
%                 iwinp=round(twinp/ndt);
%                 iwins=round(twins/ndt);
%                 
%                 for ik=1:size(stalst,1)
%                     [ncall(ik,:), ecall(ik,:),zcall(ik,:)]=loc2win2(event,iwins(1,ik),iwins(end,ik),stalst(ik,:),dir2);
%                 end
%                 % Keep waveforms not empty
%                 idkeep=find((abs(sum(ncall,2))>0)+(abs(sum(ecall,2))>0));
%                 nc=ncall(idkeep,:);
%                 ec=ecall(idkeep,:);
%                 zc=zcall(idkeep,:);
%                 stalstk=stalst(idkeep,:);
%                 nsk=size(stalstk,1);
%                 secplot(nc,ec,zc,3,stalstk,0,'N','E','Z', 3);
                % -------------------------------------------------------
                
                
                max(Stime)-min(Stime)
                if max(Stime)-min(Stime)<60
                % Plot for all stations in same window, unaligned -----------
                disp('Plotting all stations')
                
                imid=fix(median(Stime))+1;
                
                iwins_min=fix((median(Stime)-30)/ndt)+1;
                iwins_max=fix((median(Stime)+30)/ndt)+1;
                tvec=(iwins_min-1:iwins_max-1).*ndt;
                
                for ik=1:size(stalst,1)
                    [ncall(ik,:), ecall(ik,:),zcall(ik,:)]=loc2win2(event,iwins_min,iwins_max,stalst(ik,:),dir2);
                end
                % Keep waveforms not empty
                idkeep=find((abs(sum(ncall,2))>0)+(abs(sum(ecall,2))>0));
                nc=ncall(idkeep,:);
                ec=ecall(idkeep,:);
                zc=zcall(idkeep,:);
                stalstk=stalst(idkeep,:);
                nsk=size(stalstk,1);
                [s1,s2,s3]=secplot(nc,ec,zc,3,stalstk,(iwins_min-1)*ndt,'N','E','Z', 3);
                t1a=min(tpp)+Otime-6;
                t2a=max(tss)+Otime+6;
                axes(s1); hold on; vline(t1a,'r--'); vline(t2a,'r--');
                axes(s2); hold on; vline(t1a,'r--'); vline(t2a,'r--');
                axes(s3); hold on; vline(t1a,'r--'); vline(t2a,'r--');
                t1=t1a-(iwins_min-1)*ndt;
                t2=t2a-(iwins_min-1)*ndt;
                % -------------------------------------------------------
                
                end
                
                
                % Attempt alignment --------------------------------------
                %         [tdel,sigr,cc]= perfalign( [nc;ec], 1 );
                %         secplot(shift(nc,ndt,tdel(1:nsk)),shift(ec,ndt,tdel(nsk+1:end)),zc,4,stalstk,0,'N','E','Z', 2);
                %         [ix,tdel,rmean,sigr,cc,cmax] = ccsearch([nc;ec],1,0.3);
                %         in=ix(find(ix<=nsk)); tdeln=tdel(ix<=nsk);
                %         ie=ix(find(ix>nsk))-nsk; tdele=tdel(ix>nsk);
                %         secplot(shift(nc(in,:),ndt,tdeln),shift(ec(ie,:),ndt,tdele),zc,4,stalstk,0,'N','E','Z', 2);
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
                
                %-------------- SELECT DETECTIONS FOR TEMPLATE --------------
%                 ikeep=input('Keep? (1 yes, 0 no) ');
%                 if ikeep==1
%                     ix=input('Enter indices for stations kept: ');
%                     ordlst=stalstk(ix,:);
%                     uk=[nc(ix,:);ec(ix,:);zc(ix,:)];
%                     Sdetect=struct('Latitude',lat,'Longitude',lon,'Depth',zhyp,'HorizontalError',HorizErr,'OrigStations',ordlst1,'origNcomp',ncomp,'origEcomp',ecomp,'origZcomp',zcomp);
%                     dirtemp=['/fiore2_shared/shared/Data5/CASC/savardge/NCC/Template/WF0/',event,'/'];
%                     if ~exist(dirtemp)
%                         mkdir(dirtemp)
%                     end
%                     sdfm=input('Enter template number: ');
%                     fnametemp=fullfile(dirtemp,['wf_A',sprintf('%03d',sdfm),'.mat'])
%                     save(fnametemp,'ordlst','uk','t1','t2','Sdetect')
%                     fidtemp=fopen('/fiore2_shared/shared/Data5/CASC/savardge/NCC/Template/WF0/Mytemplates.dat','a');
%                     fprintf(fidtemp,'%s %d %d %s %02d %d \n',sprintf('%03d',sdfm),lat, lon, event, fix(Otime/3600), Otime-3600*fix(Otime/3600))
%                     fclose(fidtemp)
%                     h=figure(20); print(h,fullfile(dirtemp,['map_A',sprintf('%03d',sdfm),'.fig']))
%                     h=figure(3); print(h,fullfile(dirtemp,['wf_A',sprintf('%03d',sdfm),'.fig']))
%                 end
                % ----------------------------------------------------------------
                
                
                
            end
            
            clear ordlst1 ordlst uk ncomp ecomp zcomp Stime ncall ecall zcall nc ec zc
            is=0;
            tline
            disp('---------------------------------------------------------')
            
        end
    end
    
    
end