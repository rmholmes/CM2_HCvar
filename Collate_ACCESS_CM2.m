% Collate .mat files created by Process_ACCESS_CM2.m into a single
% .mat file. 

base = '/scratch/e14/rmh561/access-cm2/HCvar/';
name = 'PIcontrolTb05_';
Tyz = 0;
MOCyz = 0;
TAU = 1;

if (~Tyz & ~MOCyz & ~TAU)

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

elseif (Tyz)

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

elseif (MOCyz)

files = dir([base 'CM2_' name '_*_MOCyz.mat']);

% Fix order:
nums = [];
for fi=1:length(files)
    fname = files(fi).name;
    num = str2num(fname(end-17:end-10));
    nums = cat(1,nums,num);
end

nums = sort(nums);
 
% MOCyz collate:
load(sprintf([base 'CM2_' name '_%08d_MOCyz.mat'],nums(1)));
MOCyzfields = fieldnames(MOCyzS);

for fi = 2:length(nums)
    fname = sprintf([base 'CM2_' name '_%08d_MOCyz.mat'],nums(fi))
    next = load(fname);

    for vi = 1:length(MOCyzfields)
        eval(['MOCyzS.' MOCyzfields{vi} ' = cat(3,MOCyzS.' MOCyzfields{vi} ...
              ',next.MOCyzS.' MOCyzfields{vi} ');']);
    end
end

clear next files MOCyzfields fname fi vi;

save([base 'CM2_' name '_MOCyz_ALL.mat']);

elseif (TAU)

files = dir([base 'CM2_' name '_*_taux.mat']);

% Fix order:
nums = [];
for fi=1:length(files)
    fname = files(fi).name;
    num = str2num(fname(end-16:end-9));
    nums = cat(1,nums,num);
end

nums = sort(nums);
 
% taux collate:
load(sprintf([base 'CM2_' name '_%08d_taux.mat'],nums(1)));
tauxfields = fieldnames(tauxS);

for fi = 2:length(nums)
    fname = sprintf([base 'CM2_' name '_%08d_taux.mat'],nums(fi))
    next = load(fname);

    for vi = 1:length(tauxfields)
        eval(['tauxS.' tauxfields{vi} ' = cat(2,tauxS.' tauxfields{vi} ...
              ',next.tauxS.' tauxfields{vi} ');']);
    end
end

clear next files tauxfields fname fi vi;

save([base 'CM2_' name '_taux_ALL.mat']);

end
