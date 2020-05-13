function [ ncomp,ecomp,zcomp , ordlst, HorizErr, Otime, StimeN,StimeE] = extractDetectionSP( event, yr,month,day, hr, mn, sec, shw, dir1 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
format compact
addpath('/home/savardge/')
addpath('/home/savardge/TOOLBOXES/grep/')
addpath('/home/savardge/TOOLBOXES/subtightplot/')


ndt=0.025;
% dir1='/home/savardge/hypfiles/VI/interETS/';
cd(dir1)
% cd('/home/savardge/hypfiles/NW/')
% cd('/home/savardge/hypfiles/old_Washington/NW_obs/')
dir2='/mnt/data4/data/bostock/CASC/Data/Stations/';
% dir2='/mnt/data4/data/bostock/CAFE/Data/Stations/'
% event=[sprintf('%02d',yr),sprintf('%02d',month),sprintf('%02d',day)]; event=event(3:end);
fname=['tt_',event,'_withP_relax.arc'];
masterstationlist='    ';
expr=['20',event,sprintf('%02d',hr),sprintf('%02d',mn),sprintf('%02.f',floor(sec))]
[~,fb]=grep('-n ', expr, fullfile(dir1,fname));
changed=0;
if isempty(fb.line)
    m=1;
    while isempty(fb.line)
        m
        if m>20 && changed==0;
            if day==30
                event2=[num2str(yr),sprintf('%02d',month+1),sprintf('%02d',1)]; event2=event2(3:end);
            else
                event2=[num2str(yr),sprintf('%02d',month),sprintf('%02d',day+1)]; event2=event2(3:end);
            end
            fname=['tt_',event2,'_svi.arc']
            m=0; changed=1;
        end
        if m>20 && changed==1;
            event2=[num2str(yr),sprintf('%02d',month),sprintf('%02d',day-1)]; event2=event2(3:end);
            fname=['tt_',event2,'_svi.arc']
            m=0; changed=2;
        end
        tt=(mn*60+sec)+m;
        if floor(tt/3600)>0 && hr==23
            expr=['20',event,sprintf('%02d',0),sprintf('%02d',floor((tt-3600*floor(tt/3600))/60)),sprintf('%02.f',floor(tt-60*floor((tt-3600*floor(tt/3600))/60))-3600*floor(tt/3600))]
        end
        expr=['20',event,sprintf('%02d',hr+floor(tt/3600)),sprintf('%02d',floor((tt-3600*floor(tt/3600))/60)),sprintf('%02.f',floor(tt-60*floor((tt-3600*floor(tt/3600))/60))-3600*floor(tt/3600))]
        [~,fb]=grep('-n ', expr, fullfile(dir1,fname));
        fb.line;
        
        if isempty(fb.line)
            tt=(mn*60+sec)-m;
            expr=['20',event,sprintf('%02d',hr+floor(tt/3600)),sprintf('%02d',floor((tt-3600*floor(tt/3600))/60)),sprintf('%02.f',floor(tt-60*floor((tt-3600*floor(tt/3600))/60))-3600*floor(tt/3600))]
            [~,fb]=grep('-n ', expr, fullfile(dir1,fname));
            fb.line;
        end
        m=m+1;
        if m>30
            keyboard
        end
    end
end
fb.line
fid=fopen(fullfile(dir1,fname),'r');
flag=0;
issn=0;isse=0;isp=0;
Ptime=[];
StimeN=[];
StimeE=[];

% keyboard

for k=1:fb.line
    tline = fgetl(fid);
end
for k=1:50
    
    if isempty(tline) && flag==1
        flag=0;
        
        
        if exist('ordlstSN','var') && exist('ordlstSE','var')
            ordlst=unique([ordlstSN;ordlstSE;ordlstP],'rows');
        elseif exist('ordlstSN','var')
            ordlst=unique([ordlstSN;ordlstP],'rows');
        elseif exist('ordlstSE','var')
            ordlst=unique([ordlstSE;ordlstP],'rows');
        end
        
        % istrt=fix((min(Ptime)-2)/ndt)+1;
        % iend=fix((max([StimeN(:);StimeE(:)])+2)/ndt)+1;
        % Compute times and windows for stations. Take
        % z-coordinate downward.
        [xs,ys]=read_crd(ordlst);
        [xhyp,yhyp]=deg2utm(lat,-lon); xhyp=xhyp.*1e-3; yhyp=yhyp.*1e-3;
        xep=sqrt((xs-xhyp).^2+(ys-yhyp).^2);
        phi=atan2(ys-yhyp,xs-xhyp);
        % Use R=1.73 to be consistent with hyp2000 location
        [tpp,tss,~,~,~,~]=ttimes_hess(depth,xep,1.73,phi);
        twins=[tss'-max(tss-tpp)-2;tss'+2];
        iwins=round(twins/0.025);
        npad=max([iwins(2,:)-iwins(1,:)])+fix(1/ndt);
        
%         lowert=min(tpp)-1; uppert=max(tss)+1; halflen=round(roundn((uppert-lowert)/2,-1)/0.025);
%         iwins=[ (round(tss'./0.025) -halflen) ;(round(tss'./0.025) + halflen) ];
%         npad=max([iwins(2,:)-iwins(1,:)])+fix(1/ndt);
        [ncomp,ecomp,zcomp]=read_data(ordlst,event,Otime,iwins,npad);
        if sum(ncomp(1,:))==0
            keyboard
        end
        % Extraction and Filtering
%         [b,a]=butter(4,[1.5 7]/(40/2),'bandpass');
        t=zeros(size(ncomp));
        for k=1:size(ordlst,1)
%             [ncomp(k,:), ecomp(k,:),zcomp(k,:)]=loc2win2(event,iwins(1,k),iwins(2,k),ordlst(k,:),dir2);
%             ncomp(k,:)=filtfilt(b,a,ncomp(k,:));
%             ecomp(k,:)=filtfilt(b,a,ecomp(k,:));
%             zcomp(k,:)=filtfilt(b,a,zcomp(k,:));
            t(k,:)=(twins(1,k):ndt:((npad-1)*ndt+twins(1,k)))+Otime;
        end

        % Plot waveforms used in detection
        if shw==1
           
            
            tsN=zeros(size(ordlst,1),1); tp=tsN; tsE=tsN;
            if exist('ordlstSN','var')
                [flagSN,inds]=ismember(cellstr(ordlst),cellstr(ordlstSN)); %tsN(inds)=StimeN;
            end  
            if exist('ordlstSE','var')
                [flagSE,inds]=ismember(cellstr(ordlst),cellstr(ordlstSE)); %tsE(inds)=StimeE;
            end
            if exist('flagSN','var') && exist('flagSE','var')
                flagS=(flagSN+flagSE)>0;
            elseif exist('flagSN','var')
                flagS=flagSN;
            elseif exist('flagSE','var')
                flagS=flagSE;
            end
            [flagP,indp]=ismember(cellstr(ordlst),cellstr(ordlstP)); %tp(indp)=Ptime;
            
            secplot_pretty(t,ncomp,ecomp,zcomp,1,ordlst, tss+Otime,tss+Otime, tpp+Otime,flagS,flagP);
%             secplot_pretty(t,ncomp,ecomp,zcomp,1,ordlst, tsN,tsE, tp);
            
        end
        
        disp('---------------------------------------------------------')
        break
    elseif strncmp(tline,'20',2) % Location header
        %         disp('Header')
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
        VertErr=str2double(tline(90:93))/100;
        sprintf('lat: %.2f, lon: %.2f, depth: %.2f.  Horiz. error: %.2f  Vert. error: %.2f',lat,lon,depth,HorizErr,VertErr)
        
        
    elseif flag==1 && isletter(tline(1))  % Station line
        if (strcmp(tline(14:15),'IP'))
            %             disp('P time')
            %             sprintf('%s',tline)
            isp=isp+1;
            ordlstP(isp,:)=tline(1:4);
            Ptime(isp)=str2double(tline(30:34))/100+hrmin;
        else
            %             iss=iss+1;
            %             disp('S time')
            %             sprintf('%s',tline)
            
            component=tline(47);
            if component=='N'
                issn=issn+1;
                ordlstSN(issn,:)=tline(1:4);
                StimeN(issn)=str2double(tline(42:46))/100+hrmin;
            elseif component=='E'
                isse=isse+1;
                ordlstSE(isse,:)=tline(1:4);
                StimeE(isse)=str2double(tline(42:46))/100+hrmin;
            elseif component=='B'
                isse=isse+1;
                issn=issn+1;
                ordlstSN(issn,:)=tline(1:4);
                ordlstSE(isse,:)=tline(1:4);
                StimeE(isse)=str2double(tline(42:46))/100+hrmin;
                StimeN(issn)=str2double(tline(42:46))/100+hrmin;
            end
            
        end
        
    end
    tline = fgetl(fid);
end

cd('/home/savardge/Locations');

fclose(fid);

end

