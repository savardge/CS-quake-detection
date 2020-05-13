function [tpp,tss,Jp,Js,JJp,JJs,pr,theta_p,theta_s]=ttimes_hess(hdepth,xep,r,phi)
%function [tpp,tss,Jp,Js,JJp,JJs,Hp,Hs]=ttimes_hess(hdepth,xep,r,phi)

% FUNCTION [TPP,TSS]=TTIMES_HESS(HDEPTH,XEP,R)
% Function TTIMES_HESS computes traveltimes TPP and TSS for LFE at 
% depth HDEPTH using GSC 1-D velocity model for stations at 
% distances XEP. R is Vp/Vs ratio. In addition it computes the
% terms required for the computation of the Jacobian and Hessian 
% for a non-linear earthquake location.

% GSC MODEL 
depth=[0.00 1.00 6.00 30.00 45.00];
vp=[5.00  6.00  6.70  7.10 7.75 ];
%vp=[6.00  6.00  6.00  6.00 6.00 ];
vs=vp/r;
up=1./vp;
us=1./vs;

% Determine layer in which event resides.
if hdepth>= depth(end)
    id=length(depth);
else
    id=find(depth>hdepth);
end
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
  dxdpp(ip)=0;
  for iz=1:nz
    xp(ip)=xp(ip)+pp(ip)*dz(iz)/sqrt(up(iz)^2-pp(ip)^2);
    tp(ip)=tp(ip)+up(iz)^2*dz(iz)/sqrt(up(iz)^2-pp(ip)^2);
    dxdpp(ip)=dxdpp(ip)+dz(iz)*(1+pp(ip)^2/(up(iz)^2-pp(ip)^2))/sqrt(up(iz)^2-pp(ip)^2);
  end
end

% S quantities.       
psmin=0;
psmax=us(id(1)-1)-xper*(us(id(1)-1)-psmin);
ps=psmin+[0:np-1]*(psmax-psmin)/(np-1);

% S loop. 
for ip=1:np
  xs(ip)=0;
  ts(ip)=0;
  dxdps(ip)=0;
  for iz=1:nz
    xs(ip)=xs(ip)+ps(ip)*dz(iz)/sqrt(us(iz)^2-ps(ip)^2);
    ts(ip)=ts(ip)+us(iz)^2*dz(iz)/sqrt(us(iz)^2-ps(ip)^2);
    dxdps(ip)=dxdps(ip)+dz(iz)*(1+ps(ip)^2/(us(iz)^2-ps(ip)^2))/sqrt(us(iz)^2-ps(ip)^2);
  end
end

% Interpolate times at xep.
tpp=interp1(xp,tp,xep,'spline');
tss=interp1(xs,ts,xep,'spline');

% Interpolate d(xep)/dp derivatives (or should 
% we interpolate 1/(d(xep)/dp).
dpdxp=interp1(xp,1./dxdpp,xep,'spline');
dpdxs=interp1(xs,1./dxdps,xep,'spline');

% Produce Jacobian/Hessian quantities.
pr=interp1(xp,pp,xep,'spline');
Jpx=-pr.*cos(phi);
Jpy=-pr.*sin(phi);
Jpz=sqrt(up(nz)^2-pr.^2);
Jp=[Jpx';Jpy';Jpz'];
JJp=[Jp*Jp'];

%%Hessian.
%Hpxx=dpdxp.*cos(phi).^2;
%Hpyy=dpdxp.*sin(phi).^2;
%Hpxy=-dpdxp.*cos(phi).*sin(phi);
%Hpzz=dpdxp.*pr.^2./(up(nz)^2-pr.^2);
%Hpxz=dpdxp.*cos(phi).*pr./sqrt(up(nz)^2-pr.^2);
%Hpyz=dpdxp.*sin(phi).*pr./sqrt(up(nz)^2-pr.^2);
%%Hp=[Hpxx,Hpxy,Hpxz;Hpxy,Hpyy,Hpyz;Hpxz,Hpyz,Hpzz]';

%rr=sqrt(xep.^2+hdepth^2);
%Jpx2=-xep.*cos(phi)./(6*rr);
%Jpy2=-xep.*sin(phi)./(6*rr);
%Jpz2=hdepth./(6*rr);

%Hpxx2=(1-(xep.*cos(phi)).^2./rr.^2)./(6*rr);
%Hpyy2=(1-(xep.*sin(phi)).^2./rr.^2)./(6*rr);
%%Hpzz2=(1-hdepth^2./rr.^2)./(6*rr);
%Hpzz2=(1-hdepth^2./rr.^2)./(6*rr);
%Hpxz2=xep.*cos(phi)*hdepth./(6*rr.^3);
%Hpyz2=xep.*sin(phi)*hdepth./(6*rr.^3);
%Hpxy2=-xep.^2.*sin(phi).*cos(phi)./(6*rr.^3);
%% Toy example properties
%Jp=[Jpx2';Jpy2';Jpz2'];
%JJp=[Jp*Jp'];
%Hp=[Hpxx2,Hpxy2,Hpxz2;Hpxy2,Hpyy2,Hpyz2;Hpxz2,Hpyz2,Hpzz2]';
%keyboard

%nJpx=xp./sqrt(xp.^2+hdepth^2)/6;
%dpdxx2=(1./r-xp.^2./r.^3)/6;
%dpdxx3=(sqrt(up(nz)^2-pp.^2)*6).^3/(6*hdepth);
%r=sqrt(xp.^2+hdepth^2);
%dpdxx=(1./sqrt(xp.^2+hdepth^2)-(xp.^2)./(xp.^2+hdepth^2).^1.5)/6;
%hpxz0=dpdxx.*pp./sqrt(up(nz)^2-pp.^2);
%hpxz1=(hdepth.*xp./r.^3)/6;
%%plot(1./dxdpp)
%plot(pp,hpxz0)
%hold on
%plot(pp,hpxz1,'r--')
%plot(pr,dpdxp.*pr./sqrt(up(nz)^2-pr.^2),'ko')
%plot(dpdxx,'r-')
%plot(dpdxx2,'k--')
%plot(dpdxx3,'g--')
%hold off

%keyboard


% S-waves.
pr=interp1(xs,ps,xep,'spline');
Jsx=-pr.*cos(phi);
Jsy=-pr.*sin(phi);
Jsz=sqrt(us(nz)^2-pr.^2);
Js=[Jsx';Jsy';Jsz'];
JJs=[Js*Js'];


Hsxx=-dpdxs.*cos(phi).^2;
Hsyy=-dpdxs.*sin(phi).^2;
Hsxy=dpdxs.*cos(phi).*sin(phi);
Hszz=dpdxs.*pr.^2./(us(nz)^2-pr.^2);
Hsxz=-dpdxs.*cos(phi).*pr./sqrt(us(nz)^2-pr.^2);
Hsyz=-dpdxs.*sin(phi).*pr./sqrt(us(nz)^2-pr.^2);
Hs=[Hsxx,Hsxy,Hsxz;Hsxy,Hsyy,Hsyz;Hsxz,Hsyz,Hszz]';

% Incidence angle:
% GSC Velocity profile: top depth | P velocity  | S velocity
Vprof = [  0  5    2.89;
           1  6    3.46;
           6  6.7  3.87;
          30  7.1  4.10;
          45  7.75 4.47;
          65  8.1  4.67];
      
alpha_surf = Vprof(1,2); % P velocity at surface
beta_surf = Vprof(1,3); % S velocity at surface
theta_p = asind( pr .* alpha_surf ) ;
theta_s = asind( pr .* beta_surf ) ;
 
return
