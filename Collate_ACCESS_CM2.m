% Collate .mat files created by Process_ACCESS_CM2.m

base = '/scratch/e14/rmh561/access-cm2/HCvar/';
name = 'PIcontrol';

files = dir([base 'CM2_' name '_*.mat']);

load([base files(1).name]);
CINfields = fieldnames(CIN);
Tvfields = fieldnames(Tv);
Zvfields = fieldnames(Zv);
Yvfields = fieldnames(Yv);

for fi = 2:length(files)
    fname = [base files(fi).name]
    next = load([base files(fi).name]);

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


    



