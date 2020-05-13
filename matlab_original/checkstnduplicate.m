load Dsta_OR.mat

ns=size(Dstn,1);

thresh=3;
flag=zeros(ns,1);
count=1;
for k=1:ns
%     disp(Dstn(k))
    if flag(k)==0
        dist = sqrt(  (Dx - Dx(k)).^2 + (Dy - Dy(k)).^2 ); % in km
        idum = find( dist < thresh & dist > 0);
        if ~isempty(idum)
%             disp('Close stations: ')
            Dstn(idum)
            flag([k;idum])=count;
            count=count+1;
        end
       
    end
end
mymap = varycolor(count);
figure;scatter(Dx,Dy,30,mymap(flag+1,:));hold on;text(Dx,Dy,Dstn)




% D=pdist([Dx,Dy]);
% [ix,iy]=find( D<3 & D>0 );
% dupsta1 = Dstn(ix,:);
% dupsta2 = Dstn(iy,:);
% dupsta=unique([dupsta1;dupsta2]);
% [~,inds]=ismember(dupsta,Dstn);
