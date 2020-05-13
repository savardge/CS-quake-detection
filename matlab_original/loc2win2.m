function [ncomp, ecomp,zcomp]=loc2win2(event,istrt,iend,ordlst,dir2)
%  Extracts N, E and Z waveforms for station ordlst and start and end
%  time samples specified in directory dir2 for day "event". Normalize.

global ndt

nl=length(istrt:iend); % maximum number of samples to read
tmax = (24*3600)./ndt ; % index of last sample of the day

stat=deblank(ordlst);
ni=zeros(1,nl);
ei=zeros(1,nl);
zi=zeros(1,nl);
event0=[event,'0000'];

if istrt > tmax
    istrt=istrt-tmax; iend=iend-tmax;
    iti=istrt:iend;
elseif iend > tmax % window reach beyond end-of-day
    iti=istrt:tmax;
else
    iti=istrt:iend;
end

offset=(iti(1)-1)*4; % 0 if starting at sample #1
nump=length(iti); % number of samples to actually read in


% Open files.
fid=fopen([dir2,stat,'/',event0,'/ncomp.bin'],'r','b');
% fid=fopen([dir2,stat,'/',event0,'/ncomp_g.bin'],'r','b');
if fid ~= -1
    fseek(fid,offset,'bof');
    tmp=fread(fid,[1,nump],'float32');
    ni(1,1:length(tmp))=tmp;
    fclose(fid);
else
%     disp(['ncomp.bin for station ',stat,' can''t be opened'])
%     [dir2,stat,'/',event0,'/ncomp.bin']
end
fid=fopen([dir2,stat,'/',event0,'/ecomp.bin'],'r','b');
% fid=fopen([dir2,stat,'/',event0,'/ecomp_g.bin'],'r','b');
if fid ~= -1
    fseek(fid,offset,'bof');
    tmp=fread(fid,[1,nump],'float32');
    ei(1,1:length(tmp))=tmp;
    fclose(fid);
else
%     disp(['ecomp.bin for station ',stat,' can''t be opened'])
end
fid=fopen([dir2,stat,'/',event0,'/zcomp.bin'],'r','b');
% fid=fopen([dir2,stat,'/',event0,'/zcomp_g.bin'],'r','b');
if fid ~= -1
    fseek(fid,offset,'bof');
    tmp=fread(fid,[1,nump],'float32');
    zi(1,1:length(tmp))=tmp;
    fclose(fid);
else
%     disp(['zcomp.bin for station ',stat,' can''t be opened'])
end

umax=max(abs([ni(1,:),ei(1,:),zi(1,:)]));
% umax=1;
ncomp=ni(1,:)/umax;
ecomp=ei(1,:)/umax;
zcomp=zi(1,:)/umax;


end