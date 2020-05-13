fl=getAllFiles('/mnt/data4/data/bostock/CASC/Data/Stations');
events=fl;
ordlst=fl;
for i=1:size(fl,1)
    f=fl{i};
    if strcmp(f(44),'.')==1 
        continue
    end
    if strcmp(f(46),'/')==1
        events{i}=f(47:52);
        ordlst{i}=f(43:45);
    else
        events{i}=f(48:53);
        ordlst{i}=f(43:46);
    end
end
events=unique(events); ordlst=unique(ordlst);
events=char(events); ordlst=char(ordlst);