function [detect,i1s,i2s,maxwin,structD,history]=lfedetect(uk,newlst,event,hr,tt,twin,history,dupflag)
% LFEDETECT: Search for LFE detection inside search window.
% INPUT:
%   - uk: traces organized by station (according to station list newlst ) and then by component: N first, than E: uk=[ncomp;ecomp]
%   - newlst: station list (char array)
%   - event: string YYMMDD
%   - hr: hour corresponding to start time of search window.
%   - tt: time in second for start of search window. ([0,3600])
%   - twin: matrix ns-by-ns with maximum lag to look for correlation (in
%   seconds)
%   - history: matrix recording result of previous windows (start and end
%   time, and correlation of detection if present). Used to avoid double
%   counts of events.
%   - flag: ns-by-1 array with value corresponding to associated cluster of
%   station. 0 if no station closer than 3 km


i1s=0;
i2s=0;
maxwin=0;
detect=0;
structD=struct('detect',[],'time',[],'stalst',[],'ts',[],'ssc',[],'sc',[],'wt',[],'dts',[]);
sigrthresh=0.3;

% global ndt ThreshCC weight
ndt=0.025;
ThreshCC=0.43;
weight=0;

% Detection Stages.
% 1. Envelope correlation 1 - checks cc coefficients to see what stations exhibit high similarity.
%    The key here is choosing the threshold for mean correlation coefficient. 0.35 seems to work OK.
% 2. Envelope correlation 2 - computes delays for only those stations with best cc coefficients,
% otherwise times are contaminated by those stations showing no correlation. Here we need to check
% that the twin is correctly defined and working correctly.
% 3. Waveform correlation  - uses the correlation/consistency analysis to compute delay modifcations
% using waveforms and allowing only consistent times.


% Determine if we use symmetric lag window, or we target and area of
% interest
if isstruct(twin)==1 %then assymmetric windows
    twinmin = twin.minim ;
    twinmax = twin.maxim ;
end

ns=size(newlst,1);
nc=uk(1:ns,:);
ec=uk(ns+1:2*ns,:);

% Envelopes.
ne=abs(hilbert(nc.').');
ee=abs(hilbert(ec.').');

%% Do first alignment with envelopes.
% [tdel,ac,ast,cc,tc]=mccc2([ne;ee],ndt,twin,0,weight);
if isstruct(twin)==1
    [~,~,~,cc,~]=mccc_target([ne;ee],twinmin,twinmax,weight);
else
    %
    %     [~,~,~,cc,~]=mccc_env([ne(:,1:end-1);ee(:,1:end-1)],4.025,-twin./2,twin/2,0);
    if size(twin,1) ~= size(ne,1)*2
        keyboard
    end
    %     [~,~,~,cc,~]=mccc22([ne;ee],twin);
    [~,~,~,cc,~]=mccc2([ne;ee],twin,weight);
end

% Search for stations that exhibit high correlation coefficients.
[ir1,ic1]=find(cc>ThreshCC);
ix2=unique([ir1;ic1])';
ian=ix2((ix2<ns+1));
iae=ix2((ix2>ns))-ns;

if (isempty([ian,iae]) || length(unique([ian,iae]))<4 )
    % disp('[ian,iae] empty or not enough stations.')
    return
    
else
%     %
%     % Look for stations closer than 3 km in selection.
%     % Reorganize cc matrix in rows
%     ccrows=triu(cc)+triu(cc,1)';
%     ccrows(sub2ind(size(cc),1:ns,(1:ns)+ns))=zeros(1,ns);
%     ccrows(sub2ind(size(cc),(1:ns)+ns,1:ns))=zeros(1,ns);
%     
%     % Get station index from same cluster
%     ia=[ian,iae];
%     ix2_corr=ix2; % initialize new index vector
%     ia_corr=ia; % initialize new index vector
%     flagia=dupflag(ia);
%     u=unique(flagia);
%     n=histc(flagia,u);
%     ndup=u(n>1);
%     ndup=ndup(ndup>0);
%     needupdate=0;
%     if ~isempty(ndup) % there are duplicates
%         for ik=1:length(ndup)
%             idup=ia(flagia==ndup(ik));
%             idup2=ix2(flagia==ndup(ik));
%             if length(unique(idup))>1
%                 needupdate=1; % will need to update ian and iae later
%                 %                 disp('Selected stations close to each other: ')
%                 %                 newlst(unique(idup),:)
%                 cctmp=ccrows(idup,:);
%                 cctmp(:,idup2)=0;
%                 %                 figure;imagesc(cctmp);axis xy
%                 maxcctmp=max(cctmp,[],2);
%                 [~,ikeep] = max(maxcctmp);
%                 %                 disp('Station kept: ')
%                 %                 newlst(idup(ikeep),:)
%                 
%                 iremove = idup2(idup~=idup(ikeep));
%                 ia_corr(ismember(ix2_corr,iremove))=[];
%                 ix2_corr(ismember(ix2_corr,iremove))=[];
%                 %                 disp('New station list: ')
%                 %                 newlst(ia_corr,:)
%                 
%             end
%         end
%         if needupdate==1
%             ian=ix2_corr((ix2_corr<ns+1));
%             iae=ix2_corr((ix2_corr>ns))-ns;
%             
%             %             newlst([ian,iae],:)
%         end
%     end
    
    % Reapply mccc to best traces to avoid delay time contamination from bad ones.
    if isstruct(twin)==1
        twin.minim=twin.minim([ian,iae+ns],[ian,iae+ns]);
        twin.maxim=twin.maxim([ian,iae+ns],[ian,iae+ns]);
        [tdel0,~,~,~]=mccc_target([ne(ian,:);ee(iae,:)],seis,twin.minim,twin.maxim,0);
        %         tdel0= perfalign( [ne(ian,:);ee(iae,:)],twin );
    else
        
        [tdel0,~,~,~]=mccc2([ne(ian,:);ee(iae,:)],twin([ian,iae+ns],[ian,iae+ns]),0);
        %         [tdel0,~,~,~]=mccc22([ne(ian,:);ee(iae,:)],twin([ian,iae+ns],[ian,iae+ns]));
        %         [tdel0,~,~,~]=mccc_env([ne(ian,1:end-1);ee(iae,1:end-1)],4.025,-twin([ian,iae+ns],[ian,iae+ns])./2,twin([ian,iae+ns],[ian,iae+ns])/2,0);
        %         tdel0= perfalign( [ne(ian,:);ee(iae,:)],twin([ian,iae+ns],[ian,iae+ns]) );
        %         tdel0= perfalign( [nc(ian,:);ec(iae,:)],twin([ian,iae+ns],[ian,iae+ns]) );
    end
    if ~isempty(ian)
        nc(ian,:)=shift(nc(ian,:),ndt,tdel0(1:length(ian)));
        ne(ian,:)=shift(ne(ian,:),ndt,tdel0(1:length(ian)));
    end
    if ~isempty(iae)
        ec(iae,:)=shift(ec(iae,:),ndt,tdel0(length(ian)+1:end));
        ee(iae,:)=shift(ee(iae,:),ndt,tdel0(length(ian)+1:end));
    end
end



%% Compute mask from stack of envelopes so that
% final delays are computed only on highest SNR signal.
mask=zeros(size([ne(ian,:);ee(iae,:)]));
esum=mean([ne(ian,:);ee(iae,:)]);
[maxwin,ie]=max(esum);
tabs=(ie-1)*ndt;
i1=max(1,ie-fix(4/ndt));
i2=min(size(mask,2),ie+fix(4/ndt));
mask(:,i1:i2)=1;


%% Avoid double count
% Check for overlap with previous window
i1s=tt+(i1*ndt); %i1s=round(i1s.*1e3)./(1e3);
i2s=tt+(i2*ndt); %i2s=round(i2s.*1e3)./(1e3);

hind=find(history(:,1)==1);
exit=0;
for hi=1:length(hind)
    ind=hind(hi);
    hi1s=history(ind,2);
    hi2s=history(ind,3);
    hismax=history(ind,4);
    %     disp('Previous values')
    %     disp([hi1s,hi2s,hismax])
    %     disp('Current values: ')
    %     disp([i1s,i2s,maxwin])
    v1=i1s:ndt:i2s;
    v2=hi1s:ndt:hi2s;
    
    if (hi1s>i1s)
        overlap=length(find(v1>=hi1s))*ndt;
    elseif (hi1s<i1s)
        overlap=length(find(v2>=i1s))*ndt;
    elseif hi1s==i1s
        overlap=length(v1)*ndt;
    else
        continue
    end
    
    if overlap>4 % Double count detected (overlap of 4 seconds)
        %         disp('Double count...')
        %         history
        %         keyboard;
        
        if hismax>maxwin % Previous detection better
            detect=0;
            exit=1;
            break
            %             disp('Previous detection better')
            % return
        elseif maxwin>=hismax %% THIS detection better than previous
            history(ind,1)=0;
            
            %             disp('Removing detection flag of: ')
            %             disp(history(hi,:))
            %             pause;
        end
        
    end
    
end
if exit==1
    %     disp('Next window.')
    %     pause;
    return
end


%% Fine adjustment with waveforms.
% Allow for a +/-1.0 s
% adjustment (twin=2.0 s) and select that portion of
% data that are characterized by no more than 0.3 second rms residual.
% [iy,tdel,rmean,sigr,cc,cmax]=ccsearch(mask.*[nc(ian,:);ec(iae,:)],ndt,2,0.1);
% [iy,tdel,~,sigr,~,~]=ccsearch(mask.*[nc(ian,:);ec(iae,:)],ndt,2,0.1);
[iy,tdel,~,sigr,~,~]=ccsearch([nc(ian,i1:i2);ec(iae,i1:i2)],2,sigrthresh);

if isempty(iy)
    %     disp('No station passed the ccsearch criteria')
    %     pause(2)
    return
else
    iw=[ian';iae'+ns];
    ix=iw(iy);
    
    
    % Determine which measurements come from North/East components
    % and shift if plotting desired.
    [ii]=find(ix < ns+1);
    ixn=ix(ii);
    tdeln=tdel(ii);
    [ii]=find(ix > ns);
    ixe=ix(ii)-ns;
    tdele=tdel(ii);
    
    if ~isempty(ixn)
        %         ne(ixn,:)=shift(ne(ixn,:),ndt,tdeln);
        nc(ixn,:)=shift(nc(ixn,:),ndt,tdeln);
    end
    if ~isempty(ixe)
        %         ee(ixe,:)=shift(ee(ixe,:),ndt,tdele);
        ec(ixe,:)=shift(ec(ixe,:),ndt,tdele);
    end
    
    % Determine number of stations corresponding to picks.
    ixz=unique([ixn;ixe]);
    nu=size(ixz,1);
    
    % Data for location. Accept detection only if we have traveltimes from 4
    % or more stations.
    if nu >= 4
        detect=1;
        %         disp('Detection')
        % Detection --------------------------------------------------------
        sc=zeros(1,nu);
        tdels=zeros(1,nu);
        ts=zeros(1,nu);
        wt=zeros(1,nu);
        for it=1:nu
            
            iq=find([ixn;ixe] == ixz(it));
            if size(iq,1) == 2
                sc(it)='B';
            elseif iq <= size(ixn,1)
                sc(it)='N';
            else
                sc(it)='E';
            end
            
            % Waveform delay (average for repeated stations).
            tdels(it)=mean(tdel(iq));
            
            % Compute total time from envelope and waveform delays.
            ir=([ian,iae] == ixz(it));
            ts(it)=tabs-tdels(it)-mean(tdel0(ir));
            
            % Base pick quality on residuals for hyp2000 weighting. Set
            % weight 0 = sigr < 0.1 s
            % weight 1 = 0.1 s < sigr < 0.2 s
            % weight 2 = 0.2 s < sigr < 0.3 s
            % weight 3 - reserved for bogus P time, set to 0 in hyp 2000.
            wt(it)=fix(mean(sigr(iq))/0.1);
            
        end % end loop over nu
        
        
        %% Save detection param. in struct for writing
        dum=num2str(100+hr);
        shr=dum(2:3);
        dum=num2str(100+fix(tt/60));
        smn=dum(2:3);
        ssc=tt-60*fix(tt/60);
        time=['20',event(1:6),shr,smn];
        structD=struct('detect',detect,'time',time,'stalst',newlst(ixz,:),'ts',ts,'ssc',ssc,'sc',sc,'wt',wt,'dts',tdels);
        %         disp('Stations included in this detection: ')
        %         newlst(ixz,:)
     
        % Plot --------------------------------------------------------------
        %         if nu>5
        %             dirfig='/home/savardge/Detections_fig/';
        % %         [[newlst(ixn,:),repmat('.N',length(ixn),1)];[newlst(ixe,:),repmat('.E',length(ixe),1)]]
%         figure(1);
%         clf;
%         section([nc(ixn,:);ec(ixe,:)],hr*3600+tt,ndt,-1,[[newlst(ixn,:),repmat('.N',length(ixn),1)];[newlst(ixe,:),repmat('.E',length(ixe),1)]]);
        %         cc(sub2ind(size(cc),ir1,ic1))
        %         ccrows([ixn; ixe+ns],[ixn; ixe+ns])
        %         keyboard
        %         saveas(gcf,[event,'_',num2str(fix(hr*3600+tt)),'_Finalwf.jpg'])
        %         %             %         secplot([ne(ixn,:);ee(ixe,:)],[nc(ixn,:);ec(ixe,:)],0,1,[newlst(ixn,:);newlst(ixe,:)],0,'Envelope','waveform','',2);
        
        %         section_psd([nc(ixn,:);ec(ixe,:)],0.025,[[newlst(ixn,:),repmat('.N',length(ixn),1)];[newlst(ixe,:),repmat('.E',length(ixe),1)]]);xlim([1 10])
        
        %             saveas(gcf,fullfile(dirfig,[event,'_',num2str(fix(hr*3600+tt)),'_Finalpsd.fig']))
%         figure(2);
%         clf;
%         section([ne(ian,:);ee(iae,:)],hr*3600+tt,ndt,1,[[newlst(ian,:),repmat('.N',length(ian),1)];[newlst(iae,:),repmat('.E',length(iae),1)]])
%         vline(tabs+(hr*3600+tt),'r-.')
%         vline(i1*ndt+(hr*3600+tt),'g-.')
%         vline(i2*ndt+(hr*3600+tt),'g-.')
%         keyboard;
        %                     saveas(gcf,fullfile(dirfig,[event,'_',num2str(fix(hr*3600+tt)),'_Envelopes.fig']))
        %             clf;
        %             secplot(nc,ec,hr*3600+tt,11,newlst,0,'','','', 2)
        %             saveas(gcf, fullfile(dirfig,[event,'_',num2str(fix(hr*3600+tt)),'_Initialset.fig']))
        %         end
        %
        %                 seis=[nc(ixn,:);ec(ixe,:)];
        %
        %                 for is=1:size(seis,1)
        %                     [p,f] = periodogram(seis(is,:),hamming(size(seis,2)),256,1/ndt,'power');
        %                     Cxy=zeros(size(seis,1)-1,length(f));
        %                     k=1;
        %                     figure(2);clf;
        %                     for js=1:size(seis,1)
        %                         if js~=is
        %                              [Cxy(k,:),F] = mscohere(seis(is,:),seis(js,:),[],[],256,1/ndt);
        %                              Cxy(k,:)=Cxy(k,:)'.*p;
        %                              figure(2);hold on;plot(F,Cxy); xlim([1 10])
        %                              k=k+1;
        %                         end
        %                     end
        %                     stack(is,:)=sum(Cxy,1);
        %                     figure(10);plot(F,stack(is,:));xlim([1 10])
        %                     pause;
        %                 end
        %         figure(11);plot(F,sum(stack,1))
        %         %
        %                 for is=1:size(seis,1)
        %                     struct_PSD.Frequency=F;
        %                     struct_PSD.Power= mag2db(stack(is,:));
        %                     [Signal] = fEqualizer(seis(is,:), 1/ndt, struct_PSD);
        %                     seis2(is,:)=zscore(Signal);
        %                 end
        %
        %                 figure(5);section(seis2,0,ndt,-1,[newlst(ixn,:);newlst(ixe,:)])
        %         figure(6);section_psd(seis2,0.025,[newlst(ixn,:);newlst(ixe,:)]);
        
        
        %title('Detection-Waveform shift')
        %         secplot(nc(ixz,:),ec(ixz,:),zc(ixz,:),4,newlst(ixz,:),0,'N','E','Z',3);
        %         disp('Displaying detection')
        % PLOT S-TRANSFORM
        %                 figure(12);clf;
        %                 nn=ns %length(ixz);
        %                 dumlst=newlst(ixz,:);
        % %                 dumc=[nc(ixn,:);ec(ixe,:)];
        %                 k=1;
        %                 for ii=1:nn
        %         %             subplot(nn,3,k)
        %                     subtightplot(nn,2,k)
        % %                     [strans,t,f] = st(nc(ixz(ii),:),0,180,ndt,1);
        %                     [strans,t,f] = st(nc(ii,:),0,180,ndt,1);
        %                     imagesc(t,f,abs(strans));axis xy;colorbar;
        %                     title(['Component N, Station: ',dumlst(ii,:)])
        %
        %         %             subplot(nn,3,k+1)
        %                     subtightplot(nn,2,k+1)
        % %                     [strans,t,f] = st(ec(ixz(ii),:),0,180,ndt,1);
        %                     [strans,t,f] = st(nc(ii,:),0,180,ndt,1);
        %                     imagesc(t,f,abs(strans));axis xy;colorbar;
        %                     title(['Component E, Station: ',dumlst(ii,:)])
        %
        %                     k=k+2;
        %                 end
        
        %         pause;
        
        
        %         subplot(1,3,1)
        %         vline(i1*ndt,'r');
        %         vline(i2*ndt,'r');
        %         subplot(1,3,2)
        %         title(['Date: ',event, ' Hr: ',num2str(hr),' Se: ',num2str(tt)]);
        %         vline(i1*ndt,'r');
        %         vline(i2*ndt,'r');
        %         subplot(1,3,3)
        %         vline(i1*ndt,'r');
        %         vline(i2*ndt,'r');
        
        %         flag=plotspec([nc(ixn,:);ec(ixe,:)],1);
        %         if flag==1
        %             detect=0;
        % %             disp('Low frequency noise: next!')
        %         end
        
    else
        %         disp('Not enough stations...')
        %         pause;
    end % Detection + location
    
end % Detection.

end
