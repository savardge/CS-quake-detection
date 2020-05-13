function [tpp,tss]=ttimes_gsc(hdepth,xep,r)

% Function TTABLE GSC computes traveltimes TPP and TSS for LFE at 
% depth HDEPTH using GSC 1-D location model for stations at 
% distances XDEP. R is Vp/Vs ratio

% GSC MODEL 
depth=[0.00 1.00 6.00 30.00 45.00];
vp=[5.00  6.00  6.70  7.10 7.75 ];
vs=vp/r;
up=1./vp;
us=1./vs;

% Determine layer in which event resides.
id=find(depth>hdepth);
np=100;
% Define vertical layer transit thicknesses.
dd=[depth(1:id(1)-1),hdepth];
dz=diff(dd);
nz=length(dz);

% Set ray parameter range from pure vertical to 
% XPER percent above horizontal ray.
xper=0.001;

% P quantities.
ppmin=0;
ppmax=up(id(1)-1)-xper*(up(id(1)-1)-ppmin);
pp=ppmin+[0:np-1]*(ppmax-ppmin)/(np-1);

% P loop. 
for ip=1:np
  xp(ip)=0;
  tp(ip)=0;
  for iz=1:nz
    xp(ip)=xp(ip)+pp(ip)*dz(iz)/sqrt(up(iz)^2-pp(ip)^2);
    tp(ip)=tp(ip)+up(iz)^2*dz(iz)/sqrt(up(iz)^2-pp(ip)^2);
  end
end

% S quantities.       
psmin=0;
psmax=us(id(1)-1)-xper*(us(id(1)-1)-psmin);
ps=ppmin+[0:np-1]*(psmax-psmin)/(np-1);

% S loop. 
for ip=1:np
  xs(ip)=0;
  ts(ip)=0;
  for iz=1:nz
    xs(ip)=xs(ip)+ps(ip)*dz(iz)/sqrt(us(iz)^2-ps(ip)^2);
    ts(ip)=ts(ip)+us(iz)^2*dz(iz)/sqrt(us(iz)^2-ps(ip)^2);
  end
end

% Interpolate times.
tpp=interp1(xp,tp,xep,'spline');
tss=interp1(xs,ts,xep,'spline');
 
return
