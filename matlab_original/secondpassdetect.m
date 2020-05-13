function [detect,structD]=secondpassdetect(uk,newlst,twin,tabs)

addpath('/home/savardge/TOOLBOXES/graph_connected_components/')
detect=0;
structD=struct('stalst',[],'ts',[],'sc',[],'wt',[],'dts',[]);


% global ndt ThreshCC weight
ndt=0.025;
ThreshCC=0.50;
weight=0;
ns=size(newlst,1);
nc=uk(1:ns,:);
ec=uk(ns+1:2*ns,:);

% Envelopes.
ne=abs(hilbert(nc.').');
ee=abs(hilbert(ec.').');
[stnproj,steproj,~] = secplot_st(nc,ec,ec,[],[]);
ne=stnproj;
ee=steproj;

%% Do first alignment with envelopes.

[~,~,~,cc,~]=mccc2([ne;ee],twin,weight);

% Search for stations that exhibit high correlation coefficients.
[ir1,ic1]=find(cc>ThreshCC);
ix2=unique([ir1;ic1])';
ian=ix2((ix2<ns+1));
iae=ix2((ix2>ns))-ns;

if (isempty([ian,iae]) || length(unique([ian,iae]))<4 )
    disp('[ian,iae] empty or not enough stations.')
    return
    
else
    % Check consistency
    nval = length(ir1);
    nnode = max([ir1; ic1]);
    adj = sparse(ir1, ic1, ones(nval,1), nnode, nnode);
    lbl = graph_connected_components(adj);
    grp = accumarray(lbl', (1:nnode)', [max(lbl) 1], @(x) {x});
    [maxsize, maxidx] = max(cellfun('size', grp, 1));
    
    ix2=grp{maxidx}';
    ian=ix2((ix2<ns+1));
    iae=ix2((ix2>ns))-ns;
    
    if maxsize < 4
        disp('Not enough stations connected.')
        return
    end
    [tdel0,~,~,~]=mccc2([ne(ian,:);ee(iae,:)],twin,0);
    
    if ~isempty(ian)
        nc(ian,:)=shift(nc(ian,:),ndt,tdel0(1:length(ian)));
        ne(ian,:)=shift(ne(ian,:),ndt,tdel0(1:length(ian)));
    end
    if ~isempty(iae)
        ec(iae,:)=shift(ec(iae,:),ndt,tdel0(length(ian)+1:end));
        ee(iae,:)=shift(ee(iae,:),ndt,tdel0(length(ian)+1:end));
    end
end


%% Fine adjustment with waveforms.
secplot(nc(ian,:),ec(iae,:),0,21,'',0,'N','E','Z', 2);
[iy,tdel,~,sigr,~,~]=ccsearch([nc(ian,:);ec(iae,:)],4,2);

if isempty(iy)
        disp('No station passed the ccsearch criteria')
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
        disp('Detection')
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
            ts(it)=tabs(ixz(it))-tdels(it)-mean(tdel0(ir));
            
            % Base pick quality on residuals for hyp2000 weighting. Set
            % weight 0 = sigr < 0.1 s
            % weight 1 = 0.1 s < sigr < 0.2 s
            % weight 2 = 0.2 s < sigr < 0.3 s
            % weight 3 - reserved for bogus P time, set to 0 in hyp 2000.
            wt(it)=fix(mean(sigr(iq))/0.1);
            
        end % end loop over nu
        
        
        %% Save detection param. in struct for writing
        
        structD=struct('stalst',newlst(ixz,:),'ts',ts,'sc',sc,'wt',wt,'dts',tdels);
        %         disp('Stations included in this detection: ')
        %         newlst(ixz,:)
        
        % Plot --------------------------------------------------------------
        secplot(nc(ixz,:),ec(ixz,:),0,10,newlst(ixz,:),0,'N','E','Z', 2);
        
        %         if nu>5
        %             dirfig='/home/savardge/Detections_fig/';
        % %         [[newlst(ixn,:),repmat('.N',length(ixn),1)];[newlst(ixe,:),repmat('.E',length(ixe),1)]]
%         figure(1);
%         clf;
%         section([nc(ixn,:);ec(ixe,:)],0,ndt,-1,[[newlst(ixn,:),repmat('.N',length(ixn),1)];[newlst(ixe,:),repmat('.E',length(ixe),1)]]);
%         pause;
        %         saveas(gcf,[event,'_',num2str(fix(hr*3600+tt)),'_Finalwf.jpg'])
        %         %             %         secplot([ne(ixn,:);ee(ixe,:)],[nc(ixn,:);ec(ixe,:)],0,1,[newlst(ixn,:);newlst(ixe,:)],0,'Envelope','waveform','',2);
        %         figure(2);
        %         clf;
        %         section_psd([nc(ixn,:);ec(ixe,:)],0.025,[[newlst(ixn,:),repmat('.N',length(ixn),1)];[newlst(ixe,:),repmat('.E',length(ixe),1)]]);xlim([1 10])
        
        
        
    else
        disp('Not enough stations...')
        sigr
    end % Detection + location
    
end % Detection.
pause;
end
