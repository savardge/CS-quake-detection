clear all
close all

addpath('/home/savardge/m_map')

region='SVI'

global blat blon bsta ndt

ndt=0.025;
% Read in regional stations.
mk_casc_struct;
blat=[casc.lat];
blon=[casc.lon];
for id=1:length(blon)
    bsta(id,1:4)=[casc(id).name,blanks(4-length(casc(id).name))];
end

% Map parameters
fig1=figure(10); clf; hold on
set(fig1,'DefaultaxesFontSize',14);
set(fig1,'DefaultlineMarkerSize',10);
if strcmp(region,'NVI')
    m_proj('albers equal-area','lat',[49.5,51],'long',[-128,-125.5],'rect','on'); % NVI
elseif strcmp(region,'SVI')
    m_proj('albers equal-area', 'lat',[48.1,49.2],'long',[-124.5,-123.0],'rect','on');
elseif strcmp(region,'NW')
    m_proj('albers equal-area', 'lat',[46.0,49.0],'long',[-126.5,-121.0],'rect','on');
end
m_usercoast('CASCADIA.map','patch',[0.7 0.7 0.7],'edgecolor','black')
m_grid

% Get archive file
event='040718'
fname=['tt_',event,'.arc'];
dir='/home/savardge/hypfiles/'
fid=fopen(fullfile(dir,fname),'r');

countl=0;
countd=0;
TAB=[];

tline=fgetl(fid); % gte first line -> will be header
countl=countl+1;

while ~feof(fid)
    
    %% Header
    countd=countd+1 % detection count
    tline;
    
    
    sdate=tline(3:8);
    hr=str2double(tline(9:10));
    mn=str2double(tline(11:12));
    hrmin=hr*3600+mn*60;
    sec=str2double(tline(13:16))/100;
    Otime=hr*3600+mn*60+sec; %Origin time
    lat=str2double(tline(17:18))+(str2double(tline(20:21))+str2double(tline(22:23))/100)/60; % verify min/sec
    lon=str2double(tline(24:26))++(str2double(tline(28:29))+str2double(tline(30:31))/100)/60; % verify min/sec
    depth=str2double(tline(32:36))/100;
    % numPS=str2double(tline(40:42)) % # of P S times with final weight over 0.1
    % maazimgap=str2double(tline(43:45))
    % DistnearestSta=str2double(tline(46:48))
    % RMSres=str2double(tline(49:52)) % RMS travel time residual
    % AzimMaxErr=str2double(tline(53:55)) % azimuth of largest principal error
    % (deg E of N)
    % DipMaxErr=str2double(tline(56:57)) % Dip of largest principal error
    % SizeMaxErr=str2double(tline(58:61)) % size of largest principal error (km)
    % SzMinErr=str2double(tline(67:70))
    % NumS=str2double(tline(83:85)) %# of S times weight over 0.1
    HorizErr=str2double(tline(86:89))/100;  % Horizontal error (km)
    TAB=[TAB HorizErr];
    VertErr=str2double(tline(90:93))/100; % Vertical error (km)
    % NumPfirstmot=str2double(tline(94:96)) % # of P first motions
    NumReadings=str2double(tline(119:121)); % # of P ans S readings
    
    disp('Header read')
    
    tline=fgetl(fid); % first station line
    countl=countl+1;
    
    %% Station lines
    
%     Ssec=zeros(NumReadings-1,1);
%     SemergAngle=zeros(NumReadings-1,1);
%     Scomp=zeros(NumReadings-1,1);
    Sstalst=repmat('    ',NumReadings-1,1);
%     Sloc=struct('Psec',[],'Psta',[],'Ssec',Stimes,'Scomp',Scomp,'Sstalist',Sstalst,'lat',lat,'lon',lon,'depth',depth,'time',time);
	Sloc=struct('Psec',[],'Psta',[],'Ssec',[],'Scomp',[],'Sstalist',Sstalst,'lat',lat,'lon',lon,'depth',depth,'OrigTime',Otime);
    is=1;
    while ~isempty(tline)
    %for k=1:NumReadings
        sta=tline(1:4);
        if (strcmp(tline(14:16),'IPU'))

            Sloc.Psta=sta;
            Sloc.Ptime=str2double(tline(30:34))/100+hrmin;
            % Ptimeres=str2double(tline(35:38));
            Pepidist=str2double(tline(75:78));
            Sloc.PemergAngle=str2double(tline(79:81));
            Pazim=str2double(tline(92:94));
            
        else

            Sloc.Sstalst(is,:)=sta;
            Sloc.Stime(is)=str2double(tline(42:46))/100+hrmin;
            Sloc.Scomp(is)=tline(47); %E, N or B
            Sloc.Sepidist(is)=str2double(tline(75:78));
            Sloc.SemergAngle(is)=str2double(tline(79:81));
            Sloc.Sazim(is)=str2double(tline(92:94));
            Sloc.Simportance(is)=str2double(tline(105:108));
            
            is=is+1;
        end
        
        tline=fgetl(fid); % get next station line
        countl=countl+1;
    end
    
    %% DISPLAY location
    HorizErr
   if HorizErr < 5
        figure(fig1); clf; hold on
        if strcmp(region,'NVI')
             m_proj('albers equal-area','lat',[49.5,51],'long',[-128,-125.5],'rect','on'); % NVI
         elseif strcmp(region,'SVI')
             m_proj('albers equal-area', 'lat',[48.1,49.2],'long',[-124.5,-123.0],'rect','on');
         elseif strcmp(region,'NW')
             m_proj('albers equal-area', 'lat',[46.0,49.0],'long',[-126.5,-121.0],'rect','on');
         end

        m_usercoast('CASCADIA.map','patch',[0.7 0.7 0.7],'edgecolor','black')
        m_grid

        % add LFE location
        [-Sloc.lon, Sloc.lat]
        [x,y]=m_ll2xy(-Sloc.lon,Sloc.lat,'clip','point');
        sk=plot(x,y,'r*'); hold on
        set(sk,'MarkerSize',10)
    %     text(x+0.0002,y,['( ' num2str(-Sloc.lon) ' , ' num2str( Sloc.lat) ' )']);

        % add stations
        
        for k=1:size(Sloc.Sstalst,1)
            idx=strmatch(Sloc.Sstalst(k,:),bsta);
            [x,y]=m_ll2xy(blon(idx),blat(idx),'clip','point');
            bh=plot(x,y,'kv');
            set(bh,'MarkerSize',7.0,'MarkerFaceColor','r');
            text(x+0.0002,y,Sloc.Sstalst(k,:));
%             text(x+0.0002,y-0.0003,num2str(Sloc.Sepidist(k)));
%             text(x+0.0002,y+0.0004,num2str(Sloc.Sazim(k)));
%             text(x-0.001,y,num2str(Sloc.Simportance(k)));
            
        end

        %% Display waveforms
        med=(max(Sloc.Stime)-Sloc.Ptime)/2+Sloc.Ptime;
        istrt=fix((med-10)/ndt);
        iend=fix((med+10)/ndt);
        
%         [istrt, iend];
        ns=size(Sloc.Sstalst,1);
        [ncomp, ecomp,zcomp]=loc2win(event,hr,istrt,iend,ns,Sloc.Sstalst);
        
       plotspec([ncomp;ecomp],1)
%         for is=1:ns
%             if ~all(isnan(ncomp(3,:)),2)
%             [t_trans,phi_trans,err_t_trans,err_phi_trans,dflag]=split(ncomp(is,:),ecomp(is,:),0.025,Sloc.Sazim(is),1,Sloc.Sstalst(is,:))
%             pause;
%             end
%         end
        
%         nscomp=zeros(size(ncomp));
%         escomp=zeros(size(ncomp));
%         [val, id, I2I1, I3I1, theta, psi]  = tremortest( ncomp, ecomp, zcomp , Sloc.Sstalst )
%         for is=1:ns
%             [ phi, dt, usn, use ] = splitparam( ncomp(is,:), ecomp(is,:) );
%             nscomp(is,:)=usn;
%             escomp(is,:)=use;
%             pause(3);
%         end
%         secplot(nscomp,escomp,0,5,[Sloc.Sstalst;Sloc.Sstalst],0,'Corrected N wf','Corrected E wf','',2);
%         pause;
   pause;
   end
  
   
   
   tline=fgetl(fid); % Skip one (empty) line to next header


end

disp(['Number of lines: ' num2str(countl)])
disp(['Number of locations: ' num2str(countd)])

figure; hist(TAB,100);title('Distribution of horizontal error')
