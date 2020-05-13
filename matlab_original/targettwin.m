function [mindts,maxdts] = targettwin(stalst, epi, Lx, Ly, Lz)
% TARGETTWIN Compute lag window for a parallelepiped target volume centered on
% epi (x, y, z) and with half-width Lx, Ly and Lz in the x, y and z directions.
% Z coordinated downward

load Dsta_NVISVINW.mat
   
ns=size(stalst,1);
[lia,locb]=ismember(cellstr(stalst),Dstn);
stalst=stalst(lia,:); locb=locb(lia);

% Grid
nx=10; ny=10;
[X,Y,Z]=meshgrid(linspace(epi(1)-Lx,epi(1)+Lx,nx),linspace(epi(2)-Ly,epi(2)+Ly,ny),[epi(3)-Lz,epi(3)+Lz]);

% Compute predicted travel time for an hypocenter at each of the 8 corners
% of the target volume using GSC 1D model
travel=zeros(ns,2,nx*ny*2);
% for each station do
for is=1:ns 
     % for each node of the grid do
     ic=1;
     for ix=1:nx
         for iy=1:ny
             for iz=1:2
                 [ix,iy,iz];
                 xhyp=X(ix,iy,iz);
                 yhyp=Y(ix,iy,iz);
                 zhyp=Z(ix,iy,iz);
                 xs=Dx(locb(is));
                 ys=Dy(locb(is));
                 xep=sqrt((xs-xhyp).^2+(ys-yhyp).^2); % distance station to hypocenter horizontaly
                 phi=atan2(ys-yhyp,xs-xhyp); % azimuth station to hypocenter
                 
                 % Use R=1.73 to be consistent with hyp2000 location
                 [tpp,tss,~,~,~,~]=ttimes_hess(zhyp,xep,1.73,phi);
                 travel(is,1,ic)=tpp;
                 travel(is,2,ic)=tss;
                 
                 ic=ic+1;
             end
         end
    end
    
end

% for is=1:ns
%     tmp=reshape(travel(is,2,:),2,ny,nx);
%     figure(4);imagesc(squeeze(tmp(1,:,:)))
%     pause;
% end

% Look at each pair of station and compute maximum time delay between them
pairs=combine(ns,2);
np = size(pairs,1); % # of pairs
mindts = zeros(ns,ns); % lower bound for S
maxdts = zeros(ns,ns); % upper bound for S
mindtp = zeros(ns,ns); % lower bound for P
maxdtp = zeros(ns,ns); % upper bound for P

for ip=1:np
    i1=pairs(ip,1); i2=pairs(ip,2);
    dts=travel(i1,2,:)-travel(i2,2,:);
    mindts(i1,i2)=min((dts));
    maxdts(i1,i2)=max((dts));
    
    dtp=travel(i1,1,:)-travel(i2,1,:);
    mindtp(i1,i2)=min((dtp));
    maxdtp(i1,i2)=max((dtp));
    
end
% figure;imagesc(mindts);colorbar

% value for if search on same station (+/- 0.3 s)
mindts(logical(eye(ns)))=0.025;
maxdts(logical(eye(ns)))=0.3;

mindts=mindts+triu(mindts,-1)'; % make twin symmetric
maxdts=maxdts+triu(maxdts,-1)'; % make twin symmetric

% if strcmp(mode,'double') % If ordlst is [ordlst ; ordlst]
    mindts=repmat(mindts,2,2);
    maxdts=repmat(maxdts,2,2);
% end

end