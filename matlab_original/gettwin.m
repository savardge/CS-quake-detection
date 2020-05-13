function [ twin, newordlst, pairs, lia , dupflag] = gettwin( ordlst , default,mode, region)
% function [ twin, newordlst, pairs, lia ] = gettwin( ordlst )
%GETTWIN get width of time window for mccc for each station
%   
if strcmp(region,'SVI') || strcmp(region,'NVI') || strcmp(region,'NW') % Vancouver Island and North Washington
    load Dsta_NVISVINW_AofA.mat
elseif strcmp(region,'OR') % Oregon
    load Dsta_OR.mat
elseif strcmp(region,'NCal') % North California
    load Dstations_NCal.mat
elseif strcmp(region,'SWS') % SW Shikoku
    load Dstations_SWS.mat
elseif strcmp(region,'KII') % SW Shikoku
    load Dstations_KII.mat
elseif strcmp(region,'CVI') % Central Vancouver Island
    load Seajadestn.mat
end

ns=size(ordlst,1);
[lia,locb]=ismember(cellstr(ordlst),cellstr(char(Dstn)));
if ~isempty(find(lia==0,1))
    error(['Station info unavailable for: ',ordlst(lia==0,:)])
end
newordlst=ordlst(lia,:); locb=locb(lia);
ns=size(newordlst,1);
% "Close station" flag
dupflag = clusterflag(locb);

% Compute length of time window
pairs=combine(size(newordlst,1),2);
nnz=sqrt((Dx(locb(pairs(:,1)))-Dx(locb(pairs(:,2)))).^2+(Dy(locb(pairs(:,1)))-(Dy(locb(pairs(:,2))))).^2);
D=sparse(pairs(:,1),pairs(:,2),nnz); % Matrix of distance between stations
D(max(max(pairs)),max(max(pairs)))=0; % pivot for square sparse matrix
v=3.6; % km/s

twin=1.5.*(D./v);

maxangle=45*pi/180;
minz=25;
ir=find(D<minz*tan(maxangle));
twin(ir)=((minz/cos(maxangle))-sqrt((minz*tan(maxangle)-D(ir)).^2+(minz).^2))./v;

% MULTIPLY TWIN BY 2 TO BE CONSISTENT WITH DEFINITION OF LAG USED IN MCCC.
twin=twin*2;

% Make twin symmetric
twin=twin+triu(twin,-1)';

% Put time on diagonal for same station search
twin=twin+diag(ones(1,ns).*default,0); % 0.3 seconds search if same station

if strcmp(mode,'double') % If ordlst is [ordlst ; ordlst]
    twin=repmat(twin,2,2);
end

end

