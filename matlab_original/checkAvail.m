function availStations = checkAvail(event,tollat,tollon, flagmap)
% Get available stations for specified date (event) within the boundaries
% latitude in range [tollat(1),tollat(2)] and longitude in range
% [tollon(1),tollon(2)].
%  Maps the stations if flagmap=1

% Get list of available stations
dirEvent=['/mnt/data4/data/bostock/CASC/Data/Events/',event,'0000'];

d = dir(dirEvent);
isub = [d(:).isdir]; %# returns logical vector
availStations = {d(isub).name}';
% Remove . and ..
availStations(ismember(availStations,{'.','..','.AppleDouble'})) = [];

if ~isempty(availStations)
    
    % Extract station info
    addpath('/home/savardge/Locations')
    mk_casc_struct;
    blat=[casc.lat];
    blon=[casc.lon];
    for id=1:length(blon)
        bsta(id,1:4)=[casc(id).name,blanks(4-length(casc(id).name))];
    end
    [lia,locb]=ismember(availStations,bsta);

    disp('Stations not in mk_casc_struct: ')
    disp(availStations(lia==0))
    availStations=availStations(lia==1);
    blat=blat(locb(lia)); blon=blon(locb(lia)); bsta=bsta(locb(lia),:);
    
    % Select stations within boundaries
    IN = find(inpolygon(blon,blat,tollon,tollat));
    availStations=availStations(IN,:);
    blat=blat(IN); blon=blon(IN);
    
    % Convert cell array to char array
    availStations=char(availStations);
    
    % Map available stations
    if flagmap==1
        addpath('/home/savardge/m_map')
        minlat=min(blat);
        maxlat=max(blat);
        minlon=min(blon);
        maxlon=max(blon);
        fig1=figure(100);
        clf
        set(fig1,'DefaultaxesFontSize',[14]);
        set(fig1,'DefaultlineMarkerSize',[10]);
        m_proj('albers equal-area', 'lat',[minlat,maxlat],'long',[minlon,maxlon],'rect','on');    
        m_usercoast('CASCADIA.map','patch',[0.7 0.7 0.7],'edgecolor','black')
        m_grid
        m_ruler([.1 .5],.1)
        slab_contours;
        m_line(A40(:,1)-360,A40(:,2),'Color','k');
        m_line(A30(:,1)-360,A30(:,2),'Color','k');
        m_line(A20(:,1)-360,A20(:,2),'Color','k');
        hold on
        [xs,ys]=m_ll2xy(blon,blat,'clip','point');
        bh=plot(xs,ys,'kv','MarkerFaceColor','r');
        title('Available network')
%         pause;
%         close(gcf)
    end
else
    disp(['No station available for this date: ', event])
    availStations=[];
end


end

