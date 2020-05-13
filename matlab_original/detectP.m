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
stalst=['GOWB';'PFB '; 'TWBB'; 'TSJB'; 'LZB '; 'TWKB'; 'MGCB';  'KLNB'; 'PGC ';'SHDB';'SSIB'; 'SILB'; 'SNB ';'VGZ ';'GLBC';];
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
    dir2='/fiore2_shared/shared/Data5/CASC/Data/Stations/';
    % dir2='/mnt/data4/data/bostock/CAFE/Data/Stations/'
    
    [xs,ys]=read_crd(stalst);
    
    ndt=0.025;
    [Aout,fb]=grep('-n ', 'GSC', fullfile(dir1,fname));
    flag=0;
    is=0;
    fid=fopen(fullfile(dir1,fname),'r');
    
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
                secplot_psd(ncomp,ecomp,zcomp,10,ordlst,'N','E','Z', 3)
                
                
                %         % PLOT S-TRANSFORM---------------------------------------------
                disp('Plotting the S transfroms')
                figure(2);clf;
                nn=size(ordlst,1);
                dumlst=ordlst;
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
                
                % Plot for all stations according to predicted arrival-----------
                disp('PLotting all stations for aligned arrivals')
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
                secplot(nc,ec,zc,3,stalstk,0,'N','E','Z', 3);
                % -------------------------------------------------------
                
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
                keyboard
            end
            
            clear ordlst ncomp ecomp zcomp %ncall ecall zcall
            is=0;
            disp('---------------------------------------------------------')
            
        end
    end
    
    
end