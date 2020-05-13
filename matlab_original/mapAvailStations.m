addpath('/home/savardge/TOOLBOXES/grep/')
addpath('/home/savardge/TOOLBOXES/m_map/')

% Date list
% listing=dir('/mnt/data4/data/bostock/CASC/Data/Events/09*');
listing=dir('/mnt/data4/data/bostock/CAFE/Data/Events/1108*');
events=char({listing.name}'); events=events(:,1:6);
ne=size(events,1);

% For each date, extract and map the available stations
% MAP LIMITS
region='NW '
if strcmp(region,'NVI')
    minlat=49.5;
    maxlat=50.5;
    minlon=-127.5;
    maxlon=-125.5;
elseif strcmp(region,'SVI')
    minlat=47.8;
    maxlat=49.5;
    minlon=-125;
    maxlon=-122.5;
    
elseif strcmp(region,'NW ')
    minlat=46.8;
    maxlat=48.5;
    minlon=-124.0;
    maxlon=-122.1;
elseif strcmp(region,'SW ')
    minlat=46.3;
    maxlat=48.5;
    minlon=-124.0;
    maxlon=-122.1;
end

for ie=1:ne
    event=events(ie,:)
    [ stalst, lat, lon ] = getStationsFromDate( event, '/mnt/data4/data/bostock/CAFE/Data/Events/' );
    stalst
    fig1=figure(1);
    clf
    set(fig1,'DefaultaxesFontSize',[12]);
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
    [xs,ys]=m_ll2xy(-lon,lat,'clip','point');
    bh=plot(xs,ys,'kv'); set(bh,'MarkerSize',6.0,'MarkerFaceColor','r');
    text(xs+0.0001,ys,stalst)
    title(events(ie,:))
    pause;
    
end