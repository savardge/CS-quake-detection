function  [val, id, I2I1, I3I1, theta, psi]  = tremortest( un, ue, uz , ordlst )
%TREMORTEST Test for geometric attributes inc ase of anisotropy
%
%   INPUT:
%   - un, ue, uz : ns-by-nt trace matrix for N, E and Z components
%   respectively.
%
%   CALCULATED ATTRIBUTES:
%   - Ray azimuth psi
%   - Ray angle theta
%   - Eigenvalues of covariance matrix: I1, I2, I3
%   - Eigenvectors of " " : v1,v2,v3
%
%   OUTPUT: 
%   - val: if true= tremor/LFE passed the test, 
%           otherwise, skip to next window
%   - id : index of stations that passed the test
%   - psi: ray azimuth for each station that passed
%
% Ref: Bostock & Christensen 2012

% un=un(:,800:1600)
% ue=ue(:,800:1600)
% uz=uz(:,800:1600)

ns=size(un,1);
N=size(un,2);
I2I1=zeros(ns,1);
I3I1=zeros(ns,1);
psi=zeros(ns,1);
theta=zeros(ns,1);

% % remove the mean
% un=un-repmat(mean(un,2),1,N);
% ue=ue-repmat(mean(ue,2),1,N);


%% Calculate attributes
% For each station do:
for is=1:ns
    
    % Compute 3-C covariance matrix C0 (0-lag correlation) 
    C0=cov([un(is,:)',ue(is,:)',uz(is,:)']);

    ordlst(is,:);
    
    % Diagonalize C0 into eigenvectors & eigenvalues
    [V, I]=eigs(C0);
    v1=V(:,1); v2=V(:,2); v3=V(:,3);
    I1=I(1,1); I2=I(2,2) ; I3=I(3,3);
    I2I1(is)=I2/I1;
    I3I1(is)=I3/I1;
    
%     disp(['I2/I1 = ' num2str(I2I1(is))])
%     disp(['I3/I1 = ' num2str(I3I1(is))])
    
%     figure;
%     quiver3([0,0,0],[0,0,0],[0,0,0],[v1(1)*I1,v2(1)*I2,v3(1)*I3],[v1(2)*I1,v2(2)*I2,v3(2)*I3],[v1(3)*I1,v2(3)*I2,v3(3)*I3])
%     quiver3([0,0,0],[0,0,0],[0,0,0],[v1(1),v2(1),v3(1)],[v1(2),v2(2),v3(2)],[v1(3),v2(3),v3(3)])
%     xlim([-1 1])
%      ylim([-1 1])
%      zlim([-1 1])
%      pause;

    % Estimate ray azimuth psi from v3
    psi(is)=270-atand(abs(v3(2)/v3(1)));
%     disp(['Back azimuth = ' num2str(psi(is))])
    
    % Estimate ray angle theta from v3
    theta(is)=atand(abs(sqrt((v3(1)^2)+(v3(2)^2))/v3(3)));
%     disp(['Ray angle = ' num2str(theta(is))])
    
end

%% Tremor test:  
% ( I2/I1 ~ 1 AND I3/I1 ~ 0 )   AND  ( theta < ~30 ) 
% = wavefield dominated by near vertical S-wave propagation

T1=0.5; % lower bound for I2/I1
T2=0.2; % upper bound for I3/I1
Tang=30; % upper bound for ray angle theta

test=(I2I1>0.5)&(I3I1<0.2)&(theta<20);
id=find(test);
numsta=sum(test);

disp(['# of stations that passed tremor test: ' num2str(numsta)])

figure(10); hold on
for is=1:ns
%     if numsta>3
%         plot(is,I2I1(is),'c*',is,I3I1(is),'m*')
%     else
%         plot(is,I2I1(is),'b*',is,I3I1(is),'r*')
%     end
    
    text(is,1,ordlst(is,:),'Rotation',50)
end
% set(gca,'XTick',1:ns,'XTickLabel',ordlst,'Rotation',45)
hline(T1,'b')
hline(T2,'r')


if numsta>=1 % More than 4 stations pass the test
    val=1;
else
    val=0;
end


end

