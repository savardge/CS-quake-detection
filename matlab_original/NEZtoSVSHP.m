function [ svcomp, shcomp, pcomp, rcomp, tcomp ] = NEZtoSVSHP( azim, p, ncomp, ecomp, zcomp )
%NEZtoSVSHP Rotate the horizontal components from N/E/Z to SV/SH/P
%   INPUT:
%      azim: azimuth, measured clockwise from North, in deg
%      theta: incidence angle, measured from the vertical (nadir), and from the
%      instrument to the source, in deg
%       p : ray parameter
%       ncomp, ecomp, zcomp : traces for N, E and Z component
%   OUTPUT:
%       pcomp: component of the P wave
%       svcomp: component of the SV wave
%       shcomp: component of the SH wave


% Rotate N/E to R/T ---------------------------------------------
% Rotate in fast / slow reference frame
% A = sparse(2*ns,2*ns);
% 
% diagIdx11 = sub2ind(size(A), 1:2:2*ns, 1:2:2*ns);
% A(diagIdx11) = cosd(azim);
% 
% diagIdx21 = sub2ind(size(A), 2:2:2*ns, 1:2:2*ns);
% A(diagIdx21) = -sind(azim);
% 
% diagIdx12 = sub2ind(size(A), 1:2:2*ns, 2:2:2*ns);
% A(diagIdx12) = sind(azim);
% 
% diagIdx22 = sub2ind(size(A), 2:2:2*ns, 2:2:2*ns);
% A(diagIdx22) = cosd(azim);
% 
% d = zeros(2*ns,size(ncomp,2));
% d(1:2:2*ns,:) = ncomp;
% d(2:2:2*ns,:)= ecomp;
% 
% dum = A*d;
% rcomp = dum(1:2:2*ns,:);
% tcomp = dum(2:2:2*ns,:);


rot_mat = [cos(azim),sin(azim);-sin(azim),cos(azim)];
newcomp = rot_mat * [ncomp; ecomp];
rcomp=newcomp(1,:);
tcomp=newcomp(2,:);


% Rotate R/Z to SV/P ---------------------------------------------
% See Bostock 1998 for transfer matrix M definition
% [P, SV, SH] = M * [R, T, Z]
% P and S velocities at the surface (km/s) (from GSC 1D model)
alpha = 5; 
beta =  2.89;
eta_alpha = sqrt( (alpha).^(-2) - p.^2 );
eta_beta = sqrt( (beta).^(-2) - p.^2 );

M = [ p*beta*beta/alpha,  0,  ((beta*p)^2 - 0.5 )/ (alpha*eta_alpha) ; 
        (0.5 - (beta*p)^2 )/(beta*eta_beta),  0,  - p*beta
        0, 0.5, 0];


% Apply the transfer matrix in the (frequency,slowness) domain
nt=length(ncomp)*2;
dum = M*[fft(rcomp,nt) ; fft(tcomp,nt) ; fft(zcomp,nt)];
pcomp = ifft(dum(1,:),nt); pcomp=pcomp(1:nt/2);
svcomp = ifft(dum(2,:),nt); svcomp=svcomp(1:nt/2);
shcomp = ifft(dum(3,:),nt); shcomp=shcomp(1:nt/2);

% % Accprding to Shearer p.200 (matrix applied to R and Z only)
% M = [ p*beta*beta/alpha, ((beta*p)^2 - 0.5 )/ (alpha*eta_alpha) ; 
%         (0.5 - (beta*p)^2 )/(beta*eta_beta),  - p*beta];
% 
% 
% % Apply the transfer matrix in the (frequency,slowness) domain
% dum = M*[rcomp; zcomp];
% pcomp = dum(1,:); 
% svcomp = dum(2,:); 
% 
% shcomp=tcomp;

end




% *** Old code ***

% % Back azimuths. (in deg)
% % The rule to remember to determine a back azimuth is as follows: Less than 180 degrees, add 180 degrees. More than 180 degrees, subtract 180 degrees.
% bazim=azim; % Back-azimuth
% bazim(azim < 180) = bazim(azim < 180) + 180;
% bazim(azim > 180) = bazim(azim > 180) - 180;

% %% Rotate traces from Z/N/E to P/SV/SH
% 
% % Roation matrix
% rot_mat = [cos(theta), -sin(theta)*sin(bazim), -sin(bazim)*cos(bazim) ; sin(theta), cos(theta)*sin(bazim), cos(theta)*cos(bazim) ; 0, -cos(bazim), sin(bazim)];
% 
% newcomp = rot_mat * [zcomp; ecomp; ncomp];
% pcomp=newcomp(1,:);
% svcomp=newcomp(2,:);
% shcomp=newcomp(3,:);
