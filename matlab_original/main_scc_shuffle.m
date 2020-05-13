clear all
close all
% addpath('/home/savardge') for fiore

global ndt ThreshCC weight

% Read in times identified for each cluster.
tlen=18
tinc=0.5
ThreshCC=0.4
Tlow=60
weight=0 % 0=unweighted mccc, 1=weighted mccc

% Define stations which you want to include (organized more or less from N to S).
%SVI'GOWB'
ordlst=['GOWB';'PFB ';'TWGB'; 'TWBB'; 'TSJB'; 'LZB '; 'TWKB'; 'MGCB'; 'KELB'; 'PGC ';'SHDB';'SSIB'; 'SILB'; 'SNB ';'KHVB';'SHVB';'VGZ ';'LOP3';'SOK4';'GLBC';'LCBC';'JRBC'];
% ordlst=['YOUB';'OZB ';'SOKB';'PFB ';'VGZ '; 'SNB ';'LZB ';'SOKB';'B011';'B027';'B028'];

% ordlst=['PFB ';'SOKB';'VGZ ';'LZB ';'SNB ';'PGC ';'YOUB';'BPCB';'A04A';'BMSB';'LRIV';'OPC ';'SQM ']; % 2013

% NVI
% load('/home/savardge/ordlstNVI.mat')

% load('ordlstNVI.mat')


% NW
% ordlst=['B04A';'LRIV';'OPC ';'SQM ';'BS11';'W020';'FACB';'W030';'DOSE';'C04A';'GNW ';'W040';'PL11';'W060';'N050';'N060'];

ns=size(ordlst,1);

%% Causality constraint: 
% get windows on correlation functions for mccc using Genevieve's routine.
[twin, newordlst, pairs, lia] = gettwin( ordlst, 0.3,'double' );
if newordlst ~= ordlst
    error(1,'Problem in GETTWIN')
    return
end

% Parameter initialization.
ndt=0.025;
nl=round((3600+tlen)/ndt);

% dir2='/media/genevieve/Seagate_Expansion_Drive/datafiore/'
dir2='/mnt/data4/data/bostock/CASC/Data/Stations/'
% dir2='/mnt/data4/data/bostock/NVI/Data/Stations/';
% dir2='/mnt/data4/data/bostock/CAFE/Data/Stations/'; %NW

% Display parameters.
ta=0;
tb=24;
pflag=0;

% For each date
% 2004
% events=['040711';'040710';'040709';'040708';'040707';'040723';'040724';'040725'];
%events=['040704';'040705';'040706';'040707';'040708';'040709';'040710';'040711';'040712';'040713';'040714';'040715';'040716';'040717';'040718';'040719';'040720';'040721';'040722';'040723';'040724';'040725';'040726';'040727'];
% events=['040715';'040717';'040718';'040719';'040720';'040721';'040722';'040723';'040724';'040725';'040726';'040727'];
% events=['040724';'040725';'040726';'040727'];
% events=['040712'; '040715'; '040719']
% 2005
% events=[ '050910'; '050911'; '050912'; '050913'; '050914'; '050915'; '050916'; '050917'; '050918'; '050919'; '050920'; '050921'; '050922'; '050923'; '050924'; '050925'];
% events=[ '050913'; '050914'; '050915'; '050916'];
% 2003
events=['030301'; '030302'; '030303'; '030304'; '030305'; '030306'; '030307'; '030308'; '030309'; '030310'; '030311'; '030312'];
%'030223'; '030224'; '030225'; '030226'; '030227'; '030228';'030227'; '030228';
% 2012
% events=['120830' ;' 120831' ;' 120901' ]; %' 120902' ;' 120903' ;' 120904' ;' 120905' ;' 120906' ;' 120907' ;' 120908' ;' 120909' ;' 120910' ;'1209011000 120912' ;' 120913' ;' 120914' ;' 120915' ;' 120916' ;' 120917' ;' 120918' ;' 120919' ;' 120920' ;'120921' ;' 120922' ;' 120923' ;' 120924' ;' 120925'];
% 2013
% events=['130913';'130914';'130915';'130916';'130917';'130918';'130919';'130920';'130921';'130922';'130923';'130924';'130925';'130926';'130927';'130928';'130929';'130930'; 
%     events=['131001';'131002';'131003';'131004';'131005';'131006';'131007';'131008';'131009';'131010';'131011';'131012';'131013'];

% NVI
% events=['060905';'060906';'060907';'060908';'060909';'060910';'060911'];
% events=['070301';'070302';'070303'];
% load('/home/savardge/eventsNVI.mat')

%NW
% events=['100807'; '100808'; '100809'; '100810'; '100811'; '100812'; '100813'; '100814'; '100815'; '100816'; '100817'; '100818'; '100819'; '100820'; '100821'; '100822'; '100823'; '100824'];
%events=['110804'; '110805'; '110806'; '110807'; '110808'; '110809'; '110810'; '110811'; '110812'; '110813'; '110814'; '110815'; '110816'; '110817'; '110818'; '110819'; '110820'; '110821'; '110822'; '110823'];
% events=['110808'; '110809'; '110810'; '110811'; '110812'; '110813'; '110814'; '110815'; '110816'; '110817'; '110818'; '110819'; '110820'; '110821'; '110822'; '110823'];
% events='110810'

ne=size(events,1);

% Pick random shifts for each station to decorrelate
% stashiftN=tlen.*(randi(ns,ns,1)-round(ns/2))+100;
% stashiftE=tlen.*(randi(ns,ns,1)-round(ns/2))+100;
stashiftN=120.*(-round(ns/2):round(ns/2)); stashiftN=stashiftN(randperm(ns));
stashiftE=fliplr(stashiftN);
% stashiftE=120.*fliplr([-round(ns/2):round(ns/2)]); %stashiftE=stashiftE(randperm(ns));


for ie=1:ne
    numdetect=0;
    tot=tic;
    event=events(ie,:)
    
    % Open hyp2000 traveltime file.
    fidl=fopen(['tt_',event,'_shuffle_w',num2str(weight),'.obs'],'w');
%     fidl=fopen('030303_shuffle_w0','w')
    
    numMem=ceil(tlen/tinc);
    history=zeros(numMem,4);
    S=struct('detect',[],'time',[],'stalst',[],'ts',[],'ssc',[],'sc',[],'wt',[]);
    S=repmat(S,numMem,1);
    
 
    % For each hour, read in hour's worth of data
    % to save on i/o.
    for hr=0:23
        tic;
        dum=num2str(hr+100);
        shr=dum(2:3)
        % for each station in list do
        nhour=zeros(ns,nl);
        ehour=zeros(ns,nl);
       
        for is=1:ns
            stat=deblank(ordlst(is,:));
            ni=zeros(1,nl);
            ei=zeros(1,nl);
            
            % For each detection in cluster do.
            event0=[event,'0000'];
            
            tt1=hr*3600;
            iti=(round((tt1-0)/ndt)+1:round((tt1+3600)/ndt));
            % Check for last hour and modify nump
            if hr == 23
                offset=(iti(1)-1)*4;
                nump=3600/ndt;
            else
                offset=(iti(1)-1)*4;
                nump=nl;
            end;
            
            % Open files.
            fid=fopen([dir2,stat,'/',event0,'/ncomp.bin'],'r','b');
%                 fid=fopen([dir2,event0,'/',stat,'/ncomp.bin'],'r','b');
            if fid ~= -1
                fseek(fid,offset,'bof');
                ni(1,1:nump)=fread(fid,[1,nump],'float32');
                fclose(fid);
            end
            fid=fopen([dir2,stat,'/',event0,'/ecomp.bin'],'r','b');
%                 fid=fopen([dir2,event0,'/',stat,'/ecomp.bin'],'r','b');
            if fid ~= -1
                fseek(fid,offset,'bof');
                ei(1,1:nump)=fread(fid,[1,nump],'float32');
                fclose(fid);
            end

            umax=max(abs([ni(1,:),ei(1,:),eps]));
            ni(1,:)=ni(1,:)/umax;
            ei(1,:)=ei(1,:)/umax;
            
            nhour(is,:)=ni(1,:);
            ehour(is,:)=ei(1,:);

            
        end % station loop.
        
        % SHIFT
        nhour=shift(nhour,ndt,stashiftN);
        ehour=shift(ehour,ndt,stashiftE);
        
        % For each time window.
        tt=0; %170 for hr 22 
        for kk=1:fix(3600/tinc)
            dum=num2str(1000+kk);
           
            t=3600*hr+tt:ndt:3600*hr+tt+tlen-ndt;
            
            % Extract window around desired time.
            iti=(round((tt-0)/ndt)+1:round((tt+ tlen)/ndt));
            
            ncomp=nhour(:,iti);
            ecomp=ehour(:,iti);
           
            
            uk=[ncomp;ecomp];
            
            % Plot according to several options.
         
            if pflag == 1                
                secplot(ncomp,ecomp,[],1,ordlst,0,'N','E','Z', 2)
            end
            
            
            % Clean uk from bad stations
            [ newuk,newlst,newtwin ] = cleanwf( ordlst, uk, Tlow );
            if isempty(newuk) || size(newuk,1)<8
                disp('empty or less than 4 stations ok')
               
                tt=tt+tinc;
                continue
            end
            
            % Test for detection/location.
            [detect,i1s,i2s,maxwin,structD,history]=lfedetect(newuk,newlst,event,hr,tt,newtwin,history);
            
            % Print to file "exiting" window if detect=1
            if history(1,1)==1
%                 disp('writing a detection')
                numdetect=numdetect+1;
                writeHYPO( fidl,S(1).time,S(1).ssc, S(1).ts, S(1).stalst, S(1).sc, S(1).wt );
            end
            
            % update memory matrix and struct
            history=[ history(2:end,:) ; [detect,i1s,i2s,maxwin]];
            S(1:end-1) = S(2:end);
            S(end)=structD;
            
            % Update to next window.
            tt=tt+tinc;
            
        end % End window loop
        toc
        numdetect
    end % hour loop.
    fclose(fidl);
    numdetect
    toc(tot) % 
end % event loop.
