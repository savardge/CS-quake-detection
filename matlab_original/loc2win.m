function [nc, ec,zc]=loc2win(event,tstrt,tdur,ordlst,dir2,flagshow)

% Define time window
iti=(round((tstrt)/ndt)+1:round((tstrt+tdur)/ndt));
offset=(iti(1)-1)*4;
nump=length(iti);

event0=[event,'0000'];


for is=1:ns
    stat=deblank(ordlst(is,:));
    ni=zeros(1,nl);
    ei=zeros(1,nl);
    zi=zeros(1,nl);
    
    
    % Open files.
    fid=fopen([dir2,stat,'/',event0,'/ncomp.bin'],'r','b');
%         fid=fopen([dir2,event0,'/',stat,'/ncomp.bin'],'r','b');
    if fid ~= -1
        fseek(fid,offset,'bof');
        ni(1,1:nump)=fread(fid,[1,nump],'float32');
        fclose(fid);
    end
    fid=fopen([dir2,stat,'/',event0,'/ecomp.bin'],'r','b');
%         fid=fopen([dir2,event0,'/',stat,'/ecomp.bin'],'r','b');
    if fid ~= -1
        fseek(fid,offset,'bof');
        ei(1,1:nump)=fread(fid,[1,nump],'float32');
        fclose(fid);
    end
    fid=fopen([dir2,stat,'/',event0,'/zcomp.bin'],'r','b');
%         fid=fopen([dir2,event0,'/',stat,'/zcomp.bin'],'r','b');
    if fid ~= -1
        fseek(fid,offset,'bof');
        zi(1,1:nump)=fread(fid,[1,nump],'float32');
        fclose(fid);
    end
    umax=max(abs([ni(1,:),ei(1,:),zi(1,:),eps]));
    ni(1,:)=ni(1,:)/umax;
    ei(1,:)=ei(1,:)/umax;
    zi(1,:)=zi(1,:)/umax;
    
    nc(is,:)=ni(1,:);
    ec(is,:)=ei(1,:);
    zc(is,:)=zi(1,:);
    
end % station loop.

% Normalize.
for is=1:ns
    amax=abs(max([nc(is,:),ec(is,:),zc(is,:)]));
    nc(is,:)=nc(is,:)/amax;
    ec(is,:)=ec(is,:)/amax;
    zc(is,:)=zc(is,:)/amax;
end


%% DIsplay

if flagshow==1
    % dum=num2str(tt+10000);
    % stt=dum(2:5);
    figure(2);
    subplot(1,3,1);
    section(nc,0,ndt);
    ylabel(' ');
    ylim([0,ns+1])
    %         xlim([ta,tb]);
    ax1=gca;
    set(ax1,'YTickLabel',ordlst,'YTick',(1:ns),'FontSize',14)
    subplot(1,3,2);
    section(ec,0,ndt);
    ylabel(' ');
    ylim([0,ns+1])
    %         xlim([ta,tb]);
    ax1=gca;
    set(ax1,'YTickLabel',ordlst,'YTick',(1:ns),'FontSize',14)
    % title([event0(1:6),'   ',shr,'   ',stt])
    subplot(1,3,3);
    section(zc,0,ndt);
    ylabel(' ');
    ylim([0,ns+1])
    %         xlim([ta,tb]);
    ax1=gca;
    set(ax1,'YTickLabel',ordlst,'YTick',(1:ns),'FontSize',14)
    %         pause
end


end