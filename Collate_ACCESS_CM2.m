% Collate .mat files created by Process_ACCESS_CM2.m

base = '/scratch/e14/rmh561/access-cm2/HCvar/';
name = 'PIcontrol_';
Tyz = 1;

if (~Tyz)

files = dir([base 'CM2_' name '_*.mat']);

% Fix order:
nums = [];
for fi=1:length(files)
    fname = files(fi).name;
    num = str2num(fname(end-11:end-4));
    nums = cat(1,nums,num);
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

else

files = dir([base 'CM2_' name '_*_Tyz.mat']);

% Fix order:
nums = [];
for fi=1:length(files)
    fname = files(fi).name;
    num = str2num(fname(end-15:end-8));
    nums = cat(1,nums,num);
end

nums = sort(nums);
 
% Tyz collate:
load(sprintf([base 'CM2_' name '_%08d_Tyz.mat'],nums(1)));
Tyzfields = fieldnames(TyzS);

for fi = 2:length(nums)
    fname = sprintf([base 'CM2_' name '_%08d_Tyz.mat'],nums(fi))
    next = load(fname);

    for vi = 1:length(Tyzfields)
        eval(['TyzS.' Tyzfields{vi} ' = cat(3,TyzS.' Tyzfields{vi} ...
              ',next.TyzS.' Tyzfields{vi} ');']);
    end
end

clear next files Tyzfields fname fi vi;

save([base 'CM2_' name '_Tyz_ALL.mat']);

end
