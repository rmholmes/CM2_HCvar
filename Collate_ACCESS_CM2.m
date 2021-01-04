% Collate .mat files created by Process_ACCESS_CM2.m

base = '/scratch/e14/rmh561/access-cm2/HCvar/';
name = 'PIcontrol';

files = dir([base 'CM2_' name '_*.mat']);

% Fix order:
nums = [];
for fi=1:length(files)
    fname = files(fi).name;
    if (length(fname)==26)
        num = str2num(fname(15:22));
        nums = cat(1,nums,num);
    end
end

nums = sort(nums);

load(sprintf([base 'CM2_' name '_%08d.mat'],nums(1)));
CINfields = fieldnames(CIN);
Tvfields = fieldnames(Tv);
Zvfields = fieldnames(Zv);
Yvfields = fieldnames(Yv);

for fi = 2:length(nums)
    fname = sprintf([base 'CM2_' name '_%08d.mat'],nums(fi))
    next = load(fname);

    time = cat(1,time,next.time);
    DT_A = cat(1,DT_A,next.DT_A);

    for vi = 1:length(CINfields)
        eval(['CIN.' CINfields{vi} ' = cat(1,CIN.' CINfields{vi} ...
              ',next.CIN.' CINfields{vi} ');']);
    end

    for vi = 1:length(Tvfields)
        eval(['Tv.' Tvfields{vi} ' = cat(2,Tv.' Tvfields{vi} ...
              ',next.Tv.' Tvfields{vi} ');']);
    end

    for vi = 1:length(Zvfields)
        eval(['Zv.' Zvfields{vi} ' = cat(2,Zv.' Zvfields{vi} ...
              ',next.Zv.' Zvfields{vi} ');']);
    end

    for vi = 1:length(Yvfields)
        eval(['Yv.' Yvfields{vi} ' = cat(2,Yv.' Yvfields{vi} ...
              ',next.Yv.' Yvfields{vi} ');']);
    end
end

clear next files CINfields Tvfields Yvfields Zvfields fname ...
    fi vi;

save([base 'CM2_' name '_ALL.mat']);


    



