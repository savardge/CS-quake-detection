clear all

% To run in parallel
numworkers=12;
% filecell={'main_scc_parr.m','ccsearch.m','cleanwf.m','combine.m','deg2utm.m','gettwin.m','lfedetect.m','mccc2.m','mccc22.m','perfalign.m','shift.m','graph_connected_components.m','writeHYPO.m','Dsta_NVISVINW_AofA.mat'}
filecell={'main_scc_parr.m','ccsearch.m','cleanwf.m','combine.m','deg2utm.m','gettwin.m','lfedetect.m','lfedetect2.m','mccc2.m','mccc22.m','perfalign.m','shift.m','graph_connected_components.m','writeHYPO.m','Dsta_NVISVINW_AofA.mat','Dstations_NCal.mat','Dstations_SWS.mat'}

matlabpool('open',numworkers,'AttachedFiles',filecell)
% parpool('local',numworkers,'AttachedFiles',filecell)

% global ThreshCC weight ndt

% Choose Region: SVI, NVI, NW, SW, OR, NCal
% -----------------------------------------------------------------------
region='SVI'
disp('Don''t forget to check Dstations in gettwin.m and parameters in lfedetect.m!')

% PARAMETERS
% -----------------------------------------------------------------------
tlen=18; %18
tinc=0.5;
% ThreshCC=0.4 %0.4
% weight=0; % 0=unweighted mccc, 1=weighted mccc
ndt=0.025;
nl=round((3600+tlen)/ndt);

autostationlist = 0; % 1 if station selection chosen automatically from availability, 0 else

% To plot raw waveforms
pflag=0;

% If targetting region
flagtarget = 0;
% lat=dms2deg(48,26,0);
% lon=dms2deg(236,40,0);
% [xx,yy,~] = deg2utm(lat,lon);
% epi(1)=xx.*1e-3; epi(2)=yy.*1e-3;
% epi(3)=32; %34
% Lx=10; Ly=10; Lz=4;

% Directory of data on fiore
if strcmp(region,'SVI') % Southern Vancouver Island
    dir2='/fiore2_shared/shared/Data5/CASC/Data/Stations/';
elseif strcmp(region,'NVI') % Northern Vancouver Island
    dir2='/mnt/data4/data/bostock/NVI/Data/Stations/';
elseif strcmp(region,'NW') % North Washington
    dir2='/mnt/data4/data/bostock/CAFE/Data/Stations/'; %NW
elseif strcmp(region,'OR') % Oregon
    dir2='/fiore2_shared/shared/Data5/Oregon/Data/Stations/';
elseif strcmp(region,'NCal') % Northern California
    dir2='/fiore2_shared/shared/Data5/NORC/Data/Stations/';
elseif strcmp(region,'SWS') % SW Shikoku
    dir2='/fiore2_shared/shared/Data5/SWShikoku/Data/Stations/';
elseif strcmp(region, 'KII')
    dir2='/fiore2_shared/shared/Data5/KiiChannel/Data/Stations/';
elseif strcmp(region, 'CVI')
    dir2='/nas_data/Data5/CVI/Data/Stations_time_correct/';
end


% Display parameters.

% DEFINE DATES
% -----------------------------------------------------------------------
% For each date
% 2004
% events=['040711';'040710';'040709';'040708';'040707';'040723';'040724';'040725'];
% events=['040704';'040705';'040706';'040707';'040708';'040709';'040710';'040711';'040712';'040713';'040714';'040715';'040716';'040717';'040718';'040719';'040720';'040721';'040722';'040723';'040724';'040725';'040726';'040727'];
% events=['040715';'040717';'040718';'040719';'040720';'040721';'040722';'040723';'040724';'040725';'040726';'040727'];
% events=['040724';'040725';'040726';'040727'];
% events=['040710';'040723';'040709';'040724';'040708';'040725';'040707';'040726';'040704';'040705';'040706';'040727'];
% events=['040715';'040717';'040718';'040719';'040720'] %'040711' %'040722' %'040721' %;;;]
events=['040720';'040711']

% 2005
% events=[ '050910'; '050911'; '050912'; '050913'; '050914'; '050915'; '050916'; '050917'; '050918'; '050919'; '050920'; '050921'; '050922'; '050923'; '050924'; '050925'];
% events=[ '050913'; '050914'; '050915'; '050916'];
% 2003
% events=[ '030223'; '030224';'030225'; '030226'; '030227'; '030228';'030301'; '030302'; '030303'; '030304'; '030305'; '030306'; '030307'; '030308'; '030309'; '030310'; '030311'; '030312';'030313';'030314';'030315';'030316'];
%;
% 2012
% events=['120830' ;' 120831' ;' 120901' ]; %' 120902' ;' 120903' ;' 120904' ;' 120905' ;' 120906' ;' 120907' ;' 120908' ;' 120909' ;' 120910' ;'1209011000 120912' ;' 120913' ;' 120914' ;' 120915' ;' 120916' ;' 120917' ;' 120918' ;' 120919' ;' 120920' ;'120921' ;' 120922' ;' 120923' ;' 120924' ;' 120925'];
% 2013'130914';'130915';'130916';'130917';
% events=['130918';'130919';'130920';'130921';'130922';'130923';'130924';'130925';'130926';'130927';'130928';'130929';'130930']; 
%     events=['131001';'131002';'131003';'131004';'131005';'131006';'131007';'131008';'131009';'131010';'131011';'131012';'131013'];
% events='040712'

% load availdayspast2006.mat

% NVI
% events=['060905';'060906';'060907';'060908';'060909';'060910';'060911'];
% events=['060905';'060906';'060907';'060908';'060909';'060910';'060911';'070228';'070301';'070302';'070303';'070613';'070614';'070615';'070616';'070617';'070618'];
% load('/home/savardge/eventsNVI.mat')
% events='080917'

%NW
% events=['100807'; '100808'; '100809'; '100810'; '100811'; '100812'; '100813'; '100814'; '100815';'100816'; '100817'; '100818' ;'100819'; '100820'; '100821'; '100822'; '100823'; '100824'];
% events=['110804'; '110805'; '110806'; '110807'; '110808'; '110809'; '110810'; '110811'; '110812'; '110813'; '110814'; '110815'; '110816'; '110817'; '110818'; '110819'; '110820'; '110821'; '110822'; '110823'];
% events=['110808'; '110809'; '110810'; '110811'; '110812'; '110813'; '110814'; '110815'; '110816'; '110817'; '110818'; '110819'; '110820'; '110821'; '110822'; '110823'];
% events='100820'

% Strait Juan de Fuca
% events=['080510';'080511';'080512';'080513';'080514';'080515';'080516';'080517';'080518';'080519';'080520';'080521';'080522';'080523'];
% events=['080512';'080513';'080514';'080515';'080516';'080517';'080518';'080519';'080520';'080521';'080522';'080523'];
% events=[];
% 
% load('phi_list.mat');
% load('dt_list.mat');

% SW
% events=['070112';'070113';'070114';'070115';'070116';'070117';'070118';'070119';'070120';'070121';'070122';'070123';'070124';'070125';'070126';'070127';'070128';'070129';'070130';'070131'];
% 080501 to 080527

% 
% k=1;
% for k1=1:12
%     for k2=1:30
%     	events(k,:)=['03',sprintf('%02d',k1),sprintf('%02d',k2)];
%     	k=k+1;
%     end
% end
% events=events(18:end,:);
% events=['050924';'050925';'050926';'050926';'050927';'050928';'050929';'050930';'050931';events]
% events=['060117';'060118';'060119';'060120';'060121';'060122';'060123';'060124';'060125';'060126';'060127';'060128';'060129';'060130';'060131';events];

% OR
% events=['090909']
% events=events(6:end,:);

% Northern California
% 2008 episode: 04/05 to 04/25
% clear events
%  k=1;
% for k1=4
%     for k2=17:25
%     	events(k,:)=['08',sprintf('%02d',k1),sprintf('%02d',k2)];
%     	k=k+1;
%     end
% end
% events='080410'

% SW Shikoku
% k=1;
% for k1=10
%     for k2=1:31
%     	events(k,:)=['14',sprintf('%02d',k1),sprintf('%02d',k2)];
%     	k=k+1;
%     end
% end

% Kii Channel
% A=dir('/fiore2_shared/shared/Data5/KiiChannel/Data/Events'); 
% events=char({A.name}'); events=events(3:end,1:6);

% CVI (SeaJade)
% A=dir([dir2,'MWAB']); 
% events=char({A.name}'); events=events(3:end,1:6);
events=repmat('      ',15,1);
for k=17:31
%     sprintf('%02d',k)
	events(k-16,:)=['1408',sprintf('%02d',k)]
end


nev=size(events,1);


% DEFINE STATIONS
% -----------------------------------------------------------------------
if autostationlist==0
    if strcmp(region,'SVI')
%         ordlst=['GOWB';'PFB ';'TWGB'; 'TWBB'; 'TSJB'; 'LZB '; 'TWKB'; 'MGCB'; 'KELB'; 'KLNB';'PGC ';'SHDB';'SSIB'; 'SILB'; 'SNB ';'KHVB';'SHVB';'VGZ ';'LOP3';'SOK4';'GLBC';'LCBC';'JRBC';'YOUB'];
        % ordlst=['YOUB';'OZB ';'SOKB';'PFB ';'VGZ '; 'SNB ';'LZB ';'SOKB';'B011';'B027';'B028'];
        
%         ordlst=['KELB';'LZB ';'MGCB';'PFB ';'PGC ';'SILB';'SNB ';'SSIB';'TSJB';'TWBB';'TWGB';'TWKB';'VGZ ']
%         ordlst = ['GOWB';'KELB';'KHVB';'MGCB';'PFB ';'PGC';'SHDB';'SHVB';'SILB';'SNB ';'SSIB';'TWBB';'TWGB';'TWKB';'VGZ ']%2003 used list
%         ordlst = ['GOWB';'KELB';'KLNB';'KHVB';'LCBC';'LZB ';'MGCB';'PFB ';'PGC ';'SHDB';'SHVB';'SILB';'SNB ';'SSIB';'TSJB';'TWBB';'TWGB';'TWKB';'VGZ ']; % 2004 used list
        % ordlst=['PFB ';'SOKB';'VGZ ';'LZB ';'SNB ';'PGC ';'YOUB';'BPCB';'A04A';'BMSB';'LRIV';'OPC ';'SQM ']; % 2013
        % ordlst=['SSIB';'SNB ';'GOWB';'SILB';'PGC ';'KELB';'MGCB';'TWKB'; 'LZB ';'TSJB';'TWBB';'TWGB';'PFB ';'LCBC';'JRBC';'SOKB';'GLBC';'SHVB';'KHVB';'VGZ '];
    % 2003
%     ordlst=['PFB ';'TWGB';'TWBB';'TWKB';'MGCB';'KELB';'SHDB';'SSIB';'GOWB';'SNB ';'VGZ ';'SHVB';'KHVB';'PGC ';'SILB'];
    
    % 2004
    ordlst=['PFB ';'TWGB';'TWBB';'TSJB';'LCBC';'LZB ';'TWKB';'MGCB';'KLNB';'SSIB';'SNB ';'GOWB';'SHDB';'PGC ';'KHVB';'SILB';'SHVB';'VGZ '];
    
    % 2005
%     ordlst=['PFB ';'TWBB';'TSJB';'LZB ';'TWKB';'MGCB';'KLNB';'SHDB';'SSIB';'SNB ';'SILB';'PGC ';'SHVB';'VGZ ';'SOK4';'JRBC';'LCBC'];
    
    % inter-ETS
    
    % 2013
%     ordlst=['PFB ';'SOKB';'VGZ ';'LZB ';'SNB ';'PGC ';'YOUB';'BPCB';'A04A';'BMSB';'LRIV';'OPC ';'SQM ']
    

    elseif strcmp(region,'NVI')
        % NVI
        load('/home/savardge/ordlstNVI.mat')

    elseif strcmp(region,'NW')
        % NW
        % ordlst=['B04A';'LRIV';'OPC ';'SQM ';'BS11';'W020';'FACB';'W030';'DOSE';'C04A';'GNW ';'W040';'PL11';'W060';'N050';'N060';'B001';'B007';'B013';'B943'];
        % ordlst=['OPC ';'SQM ';'BS11';'W020';'FACB';'DOSE';'GNW ';'W040';'PL11';'W060';'N050';'N060';'B001';'B007';'B943'];
        % ordlst=['OPC ';'SQM ';'BS11';'W020';'FACB';'DOSE';'W040';'PL11';'W060';'N050';'N060';'B001';'B007';'B013';'PA01';'DR01';'BH01';'CL02';'LC02';'TB01';'GC01';'OSD ';'STW ';'GNW '];
%         ordlst=['OPC ';'SQM ';'BS11';'W020';'FACB';'DOSE';'GNW ';'W040';'PL11';'W060';'N050';'N060';'B001';'B007';'B013';'PA01';'PA02';'PA03';'PA04';'PA06';'PA07';'PA09';'PA11';'PA12';'PA13';'DR01';'DR02';'DR03';'DR04';'DR05';'DR06';'DR07';'DR08';'DR10';'DR12';'BH01';'BH03';'BH04';'BH05';'BH06';'BH07';'BH08';'BH09';'BH10';'BH11';'CL02';'CL03';'CL04';'CL05';'CL06';'CL07';'CL08';'CL09';'CL10';'CL11';'CL12';'CL13';'CL14';'CL16';'CL18';'CL19';'LC02';'LC03';'LC04';'LC05';'LC06';'LC07';'LC08';'LC09';'LC10';'LC11';'LC12';'LC13';'LC14';'TB01';'TB04';'TB05';'TB06';'TB07';'TB08';'TB09';'TB10';'TB11';'TB12';'TB13';'TB14';'GC01';'GC02';'GC04';'GC05';'GC06';'GC07';'GC08';'GC09';'GC10';'GC11';'GC12';'GC13';'GC14';'OSD ';'STW '];
%         ordlst=['OPC ';'BS11';'W020';'FACB';'DOSE';'GNW ';'W040';'PL11';'W060';'N050';'N060';'B001';'B007';'B013';'PA07';'DR08';'BH08';'CL12';'LC13';'TB07';'GC14';'OSD ';'STW '];
%         ordlst=['OPC ';'BS11';'W020';'FACB';'DOSE';'W040';'PL11';'B001';'B003';'B004';'B007';'B013';'PA01';'DR01';'BH01';'CL02';'LC02';'TB01';'GC01';'STW ';'GNW '];
        % Additional stations AofA:  PA??,  DR??, LC??, BH??, CL??, TB??, SL??
        
        % south only
%         ordlst=['OSD ';'B013';'DOSE';'HDW ';'W040';'GNW ';'PL11';'W060';'OHC ';'B014';'NLWA';'OSR ';'WISH';'B017';'B018']
        
        
        % Strait Juan de Fuca
        % ordlst=['SOKB'; 'VGZ ';'B007';'B001';'OPC ';'SQM ';'BS11';'W020'];
    elseif strcmp(region,'SW')
        % SW
        % ordlst=['B001';'B007';'B013';'B943';'B04A';'BLN ';'BS01';'BS02';'BS03';'BS11';'C04A';'FACA';'FACB';'FACC';'FACD';'GNW ';'HDW ';'N010';'N015';'N020';'N025';'N030';'N040';'N050';'N060';'N070';'N080';'N085';'N090';'OPC ';'PL01';'PL02';'PL03';'PL04';'PL05';'PL11';'S010';'S015';'S020';'S030';'S040';'S050';'S055';'S060';'S065';'S070';'SMW ';'SQM ';'STW ';'W020';'W040';'W060';'W070';'WISH'];
    elseif strcmp(region,'OR')
        % OR
        % load Dsta_OR.mat
        % ordlst=char(Dstn)
    elseif strcmp(region,'NCal')
%         ordlst=['KSXB';'ME01';'YBH ';'ME03';'ME02';'KRMB';'BO39';'ME04';'ME29';'ME09';'ME10';'ME30';'ME13';'KHBB';'ME27';'WDC ';'ME16';'ME24';'ME28';'ME08';'ME26';'ME12']
     elseif strcmp(region,'SWS')
        ordlst=['GHKH';'HIYH';'IKKH';'IKTH';'INOH';'KWBH';'MISH';'NAKH';'OOTH';'OOZH';'SJOH';'SSKH';'TBEH';'TBRH';'TSMH';'TSSH';'TSYH';'UWAH';'YNDH']
     elseif strcmp(region,'KII')
          ordlst=[ 'OYT ';'TOS1';'TOS2';'MJY ';'NSM ';'OZK ';'ANK1';'BND ';'YST ';'ASAI';'SUMI';'HAGH';'IZSH';'KAMH';'KKGH';'KSMH';'MHRH';'MROH';'SGUH';'SNTH';'TOKH';'YMSH'];
    elseif strcmp(region,'CVI')
        ordlst=['FLYN';'HSNT';'ZBLS';'MLSP'; 'MTCH';'NMSH';'PLMP';'TAHS';'FHRB';'EDB ';'NTKA';'ETB ';'GDR ';'WOSB'];
    end
end

% Causality constraint: 
% -----------------------------------------------------------------------
% get windows on correlation functions for mccc using Genevieve's routine.
% if flagtarget==1
%     [twinmin,twinmax]=targettwin(ordlst,epi, Lx, Ly, Lz);
%     
%     twin=struct('minim',twinmin,'maxim',twinmax);
% else
%     [twin, newordlst, pairs, lia] = gettwin( ordlst, 0.3,'double' );
%     if size(newordlst,1) ~= size(ordlst,1)
%         error(1,'Problem in GETTWIN: Missing station information.')
%         return
% %     end
% end


parfor ie=1:nev
% for ie=1:nev
    
    tot=tic;
    event=events(ie,:)
    
%     if autostationlist==1
%         
%         % Define boundaries
%         tollon=[-124.5,-123];
%         tollat=[48.18,49.15];
%         
%         % Extract stations according to availability
%         ordlst = checkAvail(event,tollat,tollon, 0);
%         
%         if isempty(ordlst)
%             continue
%         end
%     end
    
    ns=size(ordlst,1);
    numdetect=0;

    % Open hyp2000 traveltime file.
%     if exist(['tt_',event,'_',lower(region),'.obs'],'file')~=0 % event already scanned
%         disp(['Event ',event,' already scanned. Pass.'])
%         continue
%     else
        fidl=fopen(['tt_',event,'_',lower(region),'.obs'],'w'); % a for append, w to overwrite
%     end
%     fidl=fopen('dumdum4.obs','a');
    
    
    numMem=ceil(tlen/tinc);
    history=zeros(numMem,4);
    S=struct('detect',[],'time',[],'stalst',[],'ts',[],'ssc',[],'sc',[],'wt',[],'dts',[]);
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
       
        % Extract waveforms
        % -----------------------------------------------------------------
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
            end
            
            % Open files.
            fid=fopen([dir2,stat,'/',event0,'/ncomp.bin'],'r','b');
%             fid=fopen([dir2,stat,'/',event0,'/ncomp_g.bin'],'r','b');
            if fid ~= -1
                fseek(fid,offset,'bof');
                ni(1,1:nump)=fread(fid,[1,nump],'float32');
                fclose(fid);
            end
            
            fid=fopen([dir2,stat,'/',event0,'/ecomp.bin'],'r','b');
%             fid=fopen([dir2,stat,'/',event0,'/ecomp_g.bin'],'r','b');
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
        
        % For each time window.
        % -----------------------------------------------------------------
        tt=0;% 170 for hr 22 
        for kk=1:fix(3600/tinc)
            dum=num2str(1000+kk);
           
            t=3600*hr+tt:ndt:3600*hr+tt+tlen-ndt;
%             disp(['Start time of window: ',num2str(t(1))])
            % Extract window around desired time.
            iti=(round((tt-0)/ndt)+1:round((tt+ tlen)/ndt));
            
            ncomp=nhour(:,iti);
            ecomp=ehour(:,iti);
            
            uk=[ncomp;ecomp];
            
            

            % Clean uk from bad stations
            % --------------------------------------------------------------
            [ newuk,newlst] = cleanwf( ordlst, uk); %, Tlow );
            if isempty(newuk) || size(newuk,1)<8
%                 disp('empty or less than 4 stations ok')
                tt=tt+tinc;
                continue
            end
            new_ns=size(newlst,1);  
            
            % MANUAL FILTERING
            % ------------------------------------------------------------
%             [b,a]=butter(4,[1.5 6.5]/(40/2),'bandpass');
%             
%             for is=1:size(newuk,1)
%                 newuk(is,:)=filtfilt(b,a,newuk(is,:));
%             end  
            
            % Plot if desired
            % -------------------------------------------------------------
            if pflag == 1                
                secplot(newuk(1:new_ns,:),newuk(new_ns+1:2*new_ns,:),0,11,ordlst,0,'N raw','E raw','Linearized', 2);
            end
            
            % Get lag windows
            % -------------------------------------------------------------
%             [lia,locb]=ismember(cellstr(newlst),cellstr(ordlst));
%             locb=locb(lia);
%             if flagtarget==1
%                 newtwin=twin;
%                 newtwin.minim=twin.minim([locb;(locb+ns)],[locb;(locb+ns)]);
%                 newtwin.maxim=twin.maxim([locb;(locb+ns)],[locb;(locb+ns)]);
%             else
            [newtwin, newordlst, pairs, lia, dupflag] = gettwin( newlst, 0.3,'double',region );
%                 newtwin=twin([locb;(locb+ns)],[locb;(locb+ns)]);
%             end
            if size(newtwin,1)/2 ~= size(newordlst,1)
                keyboard
            end
            
           
            
            % Test for detection/location.
           % --------------------------------------------------------------           
            [detect,i1s,i2s,maxwin,structD,history]=lfedetect(newuk,newlst,event,hr,tt,newtwin,history,dupflag);
%             [detect,i1s,i2s,maxwin,structD,history]=lfedetect(newuk,newlst,event,hr,tt,newtwin,history,dupflag);
            
            
            % Print to file "exiting" window if detect=1
            % -------------------------------------------------------------
            if history(1,1)==1
                
                numdetect=numdetect+1;
                writeHYPO( fidl,S(1).time,S(1).ssc, S(1).ts, S(1).stalst, S(1).sc, S(1).wt );
%                 Allstruct(numdetect)=S(1);
            end
            
            % update memory matrix and struct
            % -------------------------------------------------------------
            history=[ history(2:end,:) ; [detect,i1s,i2s,maxwin]];
            S(1:end-1) = S(2:end);
            S(end)=structD; % entering info from last lfedetect call
            
            % Update to next window.
            tt=tt+tinc;
            
        end % End window loop
        toc
        sprintf('event %s, hour %d, #detections: %d',event,hr,numdetect)
    end % hour loop.
    fclose(fidl);
    sprintf('event %s, #detections: %d',event,numdetect)


    sprintf('Time required to scan day %s : ',event)
    toc(tot)
end % event loop.

matlabpool('close');
% delete(gcp('nocreate'))