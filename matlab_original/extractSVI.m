fl=getAllFiles('/mnt/data4/data/bostock/CASC/Data/Stations');
h=char(fl);
ordlst='    ';
events='      ';
for i=1:size(fl,1)
    f=h(i,:);
    i;
    if strcmp(f(44),'.')==1
        continue
    end
    if strcmp(f(47),' ')==1
        ordlst=[ordlst ; [f(44:46),' ']];
    elseif strcmp(f(48),' ')==1
        ordlst=[ordlst ; f(44:47)];
    elseif strcmp(f(47),'/')==1
        events=[events; f(48:53)];
    elseif strcmp(f(48),'/')==1
        events=[events; f(49:54)];
    end
end
events=unique(cellstr(events)); ordlst=unique(cellstr(ordlst));
events=char(events); ordlst=char(ordlst);