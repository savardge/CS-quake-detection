function [xs,ys]=read_crd(stalst)

% Function to determine station coordinates (UTM km)

% load Dsta_NVISVINW_AofA
load Seajadestn
ns=size(stalst,1);
for is=1:ns
   ik=find(strcmp(deblank(stalst(is,:)),Dstn));
   if isempty(ik) || size(ik,1) ~= 1
      stalst(is)
      disp([stalst(is,:),' STATION COORDINATES NOT FOUND'])
   end
   xs(is)=Dx(ik);
   ys(is)=Dy(ik);
end
xs=xs';
ys=ys';
return
