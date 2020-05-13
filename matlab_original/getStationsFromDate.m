function [ stalst, lat, lon ] = getStationsFromDate( event , dir1)
%GETSTATIONSFROMDATE Get all stations available at specified date

% dir1='/mnt/data4/data/bostock/CASC/Data/Events/';

listing = dir([dir1,event,'0000']);
n1=size(listing,1);
stalst=[]; lat=[]; lon=[];
staexcl=[];
for i1=1:n1
   if ~ismember(listing(i1).name,{'.','..','.AppleDouble'} )
       
%        [~,P] = grep('-l','-Q',sprintf('%-4s',listing(i1).name),'/home/savardge/hypfiles/VI/svi.sta');
       [~,P] = grep('-l','-Q',sprintf('%-4s',listing(i1).name),'/home/savardge/hypfiles/NW/wash_sort.sta');
       if P.lcount ~= 0
           stalst=[stalst ; sprintf('%-4s',listing(i1).name)];
           lat = [lat ; str2double(P.result(2,16:17)) +  (str2double(P.result(2,18:25))/60)];
           lon = [lon ; str2double(P.result(2,26:29)) +  (str2double(P.result(2,31:37))/60)];
       else
            staexcl=[staexcl ; sprintf('%-4s',listing(i1).name)];
       end
   end
end

disp('Station info not found for: ')
disp(staexcl)

end

