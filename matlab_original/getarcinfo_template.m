function S = getarcinfo_template(template)
% GETARCINFO_TEMPLATE Get info on template from .arc HYPO2000 file
% Returns a struct with location, station parameters and P and S arrivals
% for a given template given in input (e.g. '002')


% Get archive file
% template='002';
fname=['ttimessvi_',template,'.arc'];
dir='/home/savardge/refine/templates/';

% Get # of stations
fid=fopen(fullfile(dir,fname),'r');
tline=fgetl(fid);
countl=0;
while ~feof(fid)
    tline=fgetl(fid);
    if isempty(tline)
        break
    end
    countl=countl+1;
end 
fclose(fid);
ns = countl;
%%%%%%%%%%%%%%%%

fid=fopen(fullfile(dir,fname),'r');

tline=fgetl(fid); % gte first line -> will be header

while ~feof(fid)
    %% Header
%     tline;
   
    
%     sdate=tline(3:8);
    hr=str2double(tline(9:10));
    mn=str2double(tline(11:12));
    hrmin=hr*3600+mn*60;
    sec=str2double(tline(13:16))/100;
    Otime=hr*3600+mn*60+sec; %Origin time
    lat=dms2deg(str2double(tline(17:18)),(str2double(tline(20:21))+str2double(tline(22:23))/100),0); % verify min/sec
    lon=dms2deg(str2double(tline(24:26)),(str2double(tline(28:29))+str2double(tline(30:31))/100),0); % verify min/sec
    depth=str2double(tline(32:36))/100;
%     numPS=str2double(tline(40:42)); % # of P S times with final weight over 0.1
%     maazimgap=str2double(tline(43:45));
%     DistnearestSta=str2double(tline(46:48));
%     RMSres=str2double(tline(49:52)); % RMS travel time residual
%     AzimMaxErr=str2double(tline(53:55)) ;% azimuth of largest principal error(deg E of N)
%     DipMaxErr=str2double(tline(56:57)); % Dip of largest principal error
%     SizeMaxErr=str2double(tline(58:61)); % size of largest principal error (km)
%     SzMinErr=str2double(tline(67:70));
%     NumS=str2double(tline(83:85)); %# of S times weight over 0.1
    HorizErr=str2double(tline(86:89))/100;  % Horizontal error (km)
    VertErr=str2double(tline(90:93))/100; % Vertical error (km)
%     NumPfirstmot=str2double(tline(94:96)); % # of P first motions
%     NumReadings=str2double(tline(119:121)); % # of P ans S readings
    
%     disp('Header read')
    
%     disp('1st station line')
    tline=fgetl(fid) ;% first station line
    
    
    %% Station lines
    % Initialization
    Stime=zeros(ns,1); % S arrival time in sec
    Ptime=zeros(ns,1); % P arrival time in sec
    epidist=zeros(ns,1); % Epidistance in km
    emergAngle=zeros(ns-1,1); % Emergence angle, from nadir
    azim=zeros(ns,1); % Azimuth angle, degrees East of North
    Scomp=repmat(' ',ns,1); % North, East or Both
    stalst=repmat('    ',ns,1); % station list
    Pflag=zeros(ns,1); % 1 if P reading available, 0 if not
    Sflag=zeros(ns,1); % 1 if S reading available, 0 if not
    
    S=struct('Ptime',Ptime,'Stime',Stime,'Pflag',Pflag,'Sflag',Sflag,'Epidistance',epidist,'EmergenceAngle',emergAngle,'Azimuth',azim,'Scomp',Scomp,'stalst',stalst,'lat',lat,'lon',lon,'depth',depth,'OriginTime',Otime,'HorizontalError',HorizErr,'VerticalError',VertErr);
    
    for is=1:ns
    
        %for k=1:NumReadings
        sta=tline(1:4);
        S.stalst(is,:)=sta;
        S.Epidistance(is)=str2double(tline(75:78));
        S.Azimuth(is)=str2double(tline(92:94));
        S.EmergenceAngle(is)=str2double(tline(79:81));
        
        % P arrival
        S.Ptime(is)=str2double(tline(30:34))/100+hrmin;
        % Ptimeres=str2double(tline(35:38));
        S.Pflag(is)=str2double(tline(17));
        
        % S arrival
        S.Stime(is)=str2double(tline(42:46))/100+hrmin;
        S.Scomp(is)=tline(47); %E, N or B
        S.Sflag(is)=str2double(tline(50));
        
        
%         disp('Next station line');
        tline=fgetl(fid); % get next station line
        
    end
    
    
    tline=fgetl(fid);
end