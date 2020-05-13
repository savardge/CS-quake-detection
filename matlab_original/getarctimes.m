clear all; clf;
format compact
addpath('/home/savardge/TOOLBOXES/grep/')
addpath('/home/savardge/TOOLBOXES/subtightplot/')
global ndt
ndt=0.025;
addpath('/home/savardge/CSdetection/')
% cd('/home/savardge/hypfiles/NW/')
% cd('/home/savardge/hypfiles/old_Washington/NW_obs/')
addpath('/home/savardge/')
stalst=['GOWB';'PFB ';'TWGB'; 'TWBB'; 'TSJB'; 'LZB '; 'TWKB'; 'MGCB'; 'KELB'; 'KLNB'; 'PGC ';'SHDB';'SSIB'; 'SILB'; 'SNB ';'KHVB';'SHVB';'VGZ ';'LOP3';'SOK4';'GLBC';'LCBC';'JRBC'];
% stalst=['GLBC';'JRBC';'KLNB';'LCBC';'LZB ';'MGCB';'PFB ';'PGC ';'SHDB';'SHVB';'SILB';'SNB ';'SSIB';'TSJB';'TWBB';'TWKB'; 'VGZ '];
% stalst=['OPC ';'SQM ';'BS11';'W020';'FACB';'DOSE';'W040';'PL11';'W060';'N050';'N060';'B001';'B007';'B013';'PA01';'DR01';'BH01';'CL02';'LC02';'TB01';'GC01';'OSD ';'STW ';'GNW '];
ns=size(stalst,1);
[xs,ys]=read_crd(stalst);

% dir1='/home/savardge/hypfiles/VI/';
%     dir1='/home/savardge/hypfiles/NW/'
dir1='/home/savardge/DetectP/';
dir2='/mnt/data4/data/bostock/CASC/Data/Stations/';
%     dir2='/mnt/data4/data/bostock/CAFE/Data/Stations/'
A=dir(fullfile(dir1,'tt_040718_withP_relax.arc'));
events={A.name}'; events=char(events); events=events(:,4:9);
ne=size(events,1);
veclat=[];veclon=[];vecdepth=[]; vecerrh=[]; vecerrv=[];

fout=fopen('/home/savardge/DetectP/KeptDetections.txt','w');

figure(2); clf; hold on;
load /home/savardge/AlexRtemplates/MGBstraightline.mat
load /home/savardge/Locations/templatehypos.mat
[ s, n, l ] = xy2sn( taby(:), tabx(:), cxinterpstr, cyinterpstr);
distalongdipTemp=-n;
distalongstrikeTemp=(1-s).*l;
hold on; scatter(distalongdipTemp,-tabdepth,40,'k','fill','d')
plot(rvals,czinterpdip,'r-',rvals,czinterpdipMcCrory,'g-');
axis([-40 40 -60 0])
ylabel('Depth [km]')
xlabel('Distance along dip [km]')

for ie=1:ne
    
    event=events(ie,:)
    fname=['tt_',event,'_withP_relax.arc']
    
    
    [Aout,fb]=grep('-n ', 'GSC', fullfile(dir1,fname));
    
    fid=fopen(fullfile(dir1,fname),'r');
    for k=1:fb.mlines
        tline = fgetl(fid);
        if ismember(k,fb.line) && strncmp(tline,'20',2)
            
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
            
            if HorizErr <3 && VertErr < 3
                [ ncomp,ecomp,zcomp , ordlst, HorizErr, Otime, Stime] = extractDetectionSP( event, num2str(['20',event(1:2)]),num2str(event(3:4)),num2str(event(5:6)), hr, mn, sec, 1, dir1 );
                
                keep=input('Good(1), ambiguous(2) or bad? ');
                if keep==1
                    verdict=1;
                    veclat=[veclat, lat];
                    veclon=[veclat, lon];
                    vecdepth=[veclat, depth]; 
                    vecerrh=[vecerrh,HorizErr];
                    vecerrv=[vecerrv,VertErr];
                    fprintf(fout,'%d \t %d \t %5.2f \t %5.2f \t %5.2f \t %5.2f \t %5.2f \n',ie,verdict,lat,lon,depth,HorizErr,VertErr);
                    [ dstr, ddip ] = ll2StrikeDip( lat, -lon );
                    figure(2); hold on; plot(ddip,-depth,'ro')
                    
                elseif keep==2
                    verdict=2;
                    veclat=[veclat, lat];
                    veclon=[veclat, lon];
                    vecdepth=[veclat, depth]; 
                    vecerrh=[vecerrh,HorizErr];
                    vecerrv=[vecerrv,VertErr];
                    fprintf(fout,'%d \t %d \t %5.2f \t %5.2f \t %5.2f \t %5.2f \t %5.2f \n',ie,verdict,lat,lon,depth,HorizErr,VertErr);
                    [ dstr, ddip ] = ll2StrikeDip( lat, -lon );
                    figure(2); hold on; plot(ddip,-depth,'bo')
                    
                else
                    verdict=0;
                    veclat=[veclat, lat];
                    veclon=[veclat, lon];
                    vecdepth=[veclat, depth]; 
                    vecerrh=[vecerrh,HorizErr];
                    vecerrv=[vecerrv,VertErr];
                    fprintf(fout,'%d \t %d \t %5.2f \t %5.2f \t %5.2f \t %5.2f \t %5.2f \n',ie,verdict,lat,lon,depth,HorizErr,VertErr);
                    [ dstr, ddip ] = ll2StrikeDip( lat, -lon );
                    figure(2); hold on; plot(ddip,-depth,'ko')
                end
%                 keyboard
                clear ncomp ecomp zcomp ordlst HorizErr Otime Stime
            end
        end
    end
end
fclose(fid);
fclose(fout);