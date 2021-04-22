% Post-processing for HC variability project

%%%%%% OPTIONS %%%
clear all;
baseMAT = '/srv/ccrc/data03/z3500785/access-cm2/';

PI_or_his = 1;
% 1 = PI-control, 0 = historical simualtion,
% 2 = hisNAT e1, 3 = hisNAT e2, 4 = hisNAT e3
% 5 = hisAER e1, 6 = hisGHG e1

% Streamline post-processing:
doBUDGET = 1;
doTyz = 0;

if (PI_or_his == 1)
% $$$     load([baseMAT 'CM2_PIcontrol_ALL.mat']);
% $$$     load('CM2_PIcontrol_SH_ALL.mat');
    load([baseMAT 'CM2_PIcontrolTb05__ALL.mat']);
    saveNAME = 'PIcontrolPP_Tb05.mat';
elseif (PI_or_his == 0)
    load([baseMAT 'CM2_hisTb05__ALL.mat']);
    saveNAME = 'hisPP_Tb05.mat';
elseif (PI_or_his == 2)
    load([baseMAT 'CM2_hisNATe1Tb05__ALL.mat']);
    saveNAME = 'hisNATe1PP_Tb05.mat';
% $$$ elseif (PI_or_his == 3)
% $$$     load([baseMAT 'CM2_hisNATe2__ALL.mat']);
% $$$     saveNAME = 'hisNATe2PP.mat';
% $$$ elseif (PI_or_his == 4)
% $$$     load([baseMAT 'CM2_hisNATe3__ALL.mat']);
% $$$     saveNAME = 'hisNATe3PP.mat';
% $$$ elseif (PI_or_his == 5)
% $$$     load([baseMAT 'CM2_hisAERe1__ALL.mat']);
% $$$     saveNAME = 'hisAERe1PP.mat';
% $$$ elseif (PI_or_his == 6)
% $$$     load([baseMAT 'CM2_hisGHGe1__ALL.mat']);
% $$$     saveNAME = 'hisGHGe1PP.mat';
end

saveMAT = 1;

%%%%%%

% Define a new percentile grid:
dP = 0.25;
Pe = 0:dP:100;
P = (Pe(2:end)+Pe(1:end-1))/2;
PL = length(P);

% Time vector:
tL = length(time);
time = time/365.25; % time in years

% $$$ % Coarse grain temperature bins (assumes linear):
% $$$ n = 1; % Coarse grain factor
% $$$ TLn = floor(TL/n);
% $$$ TLr = mod(TL,n);
% $$$ convm = zeros(TL,TLn);
% $$$ cnt = 1;
% $$$ for i=1:TLn
% $$$     convm(cnt:(cnt+n-1),i) = 1;
% $$$     cnt = cnt+n;
% $$$ end
% $$$ if (TLr>0)
% $$$     convm((end-TLr+1):end,end) = 1;
% $$$ end
% $$$ dT = n*dT;
% $$$ Te = -3:dT:(-3+TLn*dT);
% $$$ T = (Te(2:end)+Te(1:end-1))/2;
% $$$ TL = length(T);
% $$$ fns = fieldnames(Tv);
% $$$ for vi=1:length(fns)
% $$$     eval(['Tv.' fns{vi} ' = (Tv.' fns{vi} '''*convm)'';']);
% $$$ end

% Fix advection budget term by residual:
typs = {'T','Z','Y'};
vars = {'temp_submeso', 'temp_vdiffuse_diff_cbt', 'temp_nonlocal_KPP', ...
        'temp_vdiffuse_sbc','frazil_3d','sw_heat','temp_rivermix', ...
        'neutral_diffusion_temp','neutral_gm_temp', 'temp_vdiffuse_k33', 'mixdownslope_temp', ...
        'temp_sigma_diff','sfc_hflux_pme','temp_eta_smooth'};    
if (doBUDGET)
for ty = 1:length(typs)
    eval([typs{ty} 'v.temp_advection = ' typs{ty} 'v.temp_tendency;']);
    for vi=1:length(vars)
        eval([typs{ty} 'v.temp_advection = ' typs{ty} 'v.temp_advection ' ...
                       ' - ' typs{ty} 'v.' vars{vi} ';']);
    end
end

% group budget terms:
for vi = 1:length(typs)
    eval([typs{vi} 'v.TEN = ' typs{vi} 'v.temp_tendency;']);
    eval([typs{vi} 'v.ADV = ' typs{vi} 'v.temp_advection+' typs{vi} 'v.temp_submeso+' typs{vi} 'v.neutral_gm_temp;']);
    eval([typs{vi} 'v.ADVGM = ' typs{vi} 'v.neutral_gm_temp;']);
    eval([typs{vi} 'v.FOR = ' typs{vi} 'v.temp_vdiffuse_sbc+' typs{vi} 'v.frazil_3d+' typs{vi} 'v.sw_heat' ...
          '+' typs{vi} 'v.temp_rivermix+' typs{vi} 'v.sfc_hflux_pme' ...
          '+' typs{vi} 'v.temp_eta_smooth;']);
    eval([typs{vi} 'v.RMIX = ' typs{vi} 'v.neutral_diffusion_temp+' typs{vi} 'v.temp_vdiffuse_k33+' typs{vi} 'v.mixdownslope_temp' ...
          '+' typs{vi} 'v.temp_sigma_diff;']);
    eval([typs{vi} 'v.VMIX = ' typs{vi} 'v.temp_vdiffuse_diff_cbt+' typs{vi} 'v.temp_nonlocal_KPP;']);
end
vars = {vars{:},'temp_tendency','temp_advection'};
else
    vars = {vars{:},'temp_tendency','pme_river'};
end
for vi =1:length(vars)
    eval(['Tv = rmfield(Tv,''' vars{vi} ''');']);
    eval(['Yv = rmfield(Yv,''' vars{vi} ''');']);
    eval(['Zv = rmfield(Zv,''' vars{vi} ''');']);
end

% Do cumulative sums:
fnames = fieldnames(Zv);
for vi =1:length(fnames)
    eval(['Zv.' fnames{vi} '_c = cat(1,zeros(1,tL),cumsum(Zv.' ...
          fnames{vi} ',1));']);
end
fnames = fieldnames(Yv);
for vi =1:length(fnames)
    eval(['Yv.' fnames{vi} '_c = cat(1,zeros(1,tL),cumsum(Yv.' ...
          fnames{vi} ',1));']);
end
fnames = fieldnames(Tv);
for vi =1:length(fnames)
    eval(['Tv.' fnames{vi} '_c = cat(1,cumsum(Tv.' ...
          fnames{vi} ',1,''reverse''),zeros(1,tL));']);
end

OHC = squeeze(Zv.H_c(end,:))';

% $$$ % Time-integrate budget terms (prior to remapping to P):
% $$$ bvars = {'TEN_c','ADV_c','FOR_c','RMIX_c','VMIX_c'};
% $$$ for ti =1:length(typs)
% $$$     for vi=1:length(bvars)
% $$$         eval([typs{ti} 'v.' bvars{vi} ' = cumsum(' typs{ti} ...
% $$$               'v.' bvars{vi} '.*repmat(DT_A'',[length(' typs{ti} ...
% $$$               'v.' bvars{vi} '(:,1)) 1]),2);']);
% $$$     end
% $$$ end    


%%% Remap to percentile

% Define percentiles using volume:
Vtot = Zv.V_c(end,:)';
Atot = Tv.A_c(1,:)';
Tv.P   = 100*Tv.V_c./repmat(Vtot',[TL+1 1]);
Zv.P   = 100*Zv.V_c./repmat(Vtot',[zL+1 1]);
Yv.P   = 100*Yv.V_c./repmat(Vtot',[yL+1 1]);
Tv.Pa  = 100*Tv.A_c./repmat(Atot',[TL+1 1]);

NEWold = 0;
bvars = {'TEN_c','ADV_c','ADVGM_c','FOR_c','RMIX_c','VMIX_c'};

if (NEWold)

    %%% New approach: Interpolate HC
    % Interpolate onto high-resolution grids for smoothness:
    zfac = 1;
% $$$     Zv.H_c = interp1(Zv.H_c,1:1/zfac:(zL+1));Zv.V_c = interp1(Zv.V_c,1:1/zfac:(zL+1));
    Tfac = 1;
% $$$     Tv.H_c = interp1(Tv.H_c,1:1/Tfac:(TL+1));Tv.V_c = interp1(Tv.V_c,1:1/Tfac:(TL+1));
    yfac = 1;
% $$$     Yv.H_c = interp1(Yv.H_c,1:1/yfac:(yL+1));Yv.V_c = interp1(Yv.V_c,1:1/yfac:(yL+1));

% $$$     ZvP.Hp = zeros(PL+1,tL);
% $$$     YvP.Hp = zeros(PL+1,tL);
% $$$     TvP.Hp = zeros(PL+1,tL);
% $$$     TvP.Hap = zeros(PL+1,tL);
    for ti = 1:tL
        ZvP.Hp(:,ti) = interp1((Zv.P(:,ti)+(1:zL*zfac+1)'/1e10),Zv.H_c(:,ti),Pe,'linear');
        YvP.Hp(:,ti) = interp1((Yv.P(:,ti)+(1:yL*yfac+1)'/1e10),Yv.H_c(:,ti),Pe,'linear');
        TvP.Hp(:,ti) = interp1((Tv.P(:,ti)+(1:TL*Tfac+1)'/1e10),Tv.H_c(:,ti),Pe,'linear');
        TvP.Hap(:,ti) = interp1((Tv.Pa(:,ti)+(1:TL*Tfac+1)'/1e10),Tv.A_c(:,ti),Pe,'linear');
    end
    ZvP.Hp(1,:) = 0;
    YvP.Hp(1,:) = 0;
    TvP.Hp(1,:) = 0;
    TvP.Hap(1,:) = 0;
    ZvP.Hp(end,:) = Zv.H_c(end,:);
    YvP.Hp(end,:) = Yv.H_c(end,:);
    TvP.Hp(end,:) = Yv.H_c(1,:);
    TvP.Hap(end,:) = Tv.A_c(1,:);

    % Calculate temperatures from HC:
    ZvP.Tp = (ZvP.Hp(2:end,:)-ZvP.Hp(1:(end-1),:))/rho0/Cp*100/dP./repmat(Vtot',[PL 1]);
    YvP.Tp = (YvP.Hp(2:end,:)-YvP.Hp(1:(end-1),:))/rho0/Cp*100/dP./repmat(Vtot',[PL 1]);
    TvP.Tp = (TvP.Hp(2:end,:)-TvP.Hp(1:(end-1),:))/rho0/Cp*100/dP./repmat(Vtot',[PL 1]);
    TvP.Tap = (TvP.Hap(2:end,:)-TvP.Hap(1:(end-1),:))/rho0/Cp*100/dP./repmat(Atot',[PL 1]);

else
    %%% Old approach: Interpolate temperatures
    % Calculate z-temperature:
    Zv.T = Zv.H/rho0/Cp./Zv.V;
    Zv.Te = cat(1,Zv.T(1,:),(Zv.T(2:end,:)+Zv.T(1:end-1,:))/2,Zv.T(end,:));

    % Calculate y-temperature:
    Yv.T = Yv.H/rho0/Cp./Yv.V;
    Yv.Te = cat(1,Yv.T(1,:),(Yv.T(2:end,:)+Yv.T(1:end-1,:))/2,Yv.T(end,:));

    % Interpolate variables to P levels:
% $$$     TvP.Tp = zeros(PL,tL);
% $$$     TvP.Tap = zeros(PL,tL);
% $$$     ZvP.Tp = zeros(PL,tL);
% $$$     YvP.Tp = zeros(PL,tL);
    for ti = 1:tL
        TvP.Tp(:,ti) = interp1(Tv.P(:,ti)+(1:TL+1)'/1e10,Te,P,'linear');
        TvP.Tp(1,ti) = Te(end);
        TvP.Tap(:,ti) = interp1(Tv.Pa(:,ti)+(1:TL+1)'/1e10,Te,P,'linear');
        TvP.Tap(1,ti) = Te(end);
        ZvP.Tp(:,ti) = interp1(Zv.P(:,ti)+(1:zL+1)'/1e10,Zv.Te(:,ti),P,'linear');
        YvP.Tp(:,ti) = interp1(Yv.P(:,ti)+(1:yL+1)'/1e10,Yv.Te(:,ti),P,'linear');
    end
end

if (doBUDGET)
    for ti = 1:tL
        for vi=1:length(bvars)
            eval(['ZvP.' bvars{vi} '(:,ti) = interp1(Zv.P(:,ti)+(1:zL+1)''/1e10,Zv.' bvars{vi} '(:,ti),Pe,''linear'');']);
            eval(['YvP.' bvars{vi} '(:,ti) = interp1(Yv.P(:,ti)+(1:yL+1)''/1e10,Yv.' bvars{vi} '(:,ti),Pe,''linear'');']);
            eval(['TvP.' bvars{vi} '(:,ti) = interp1(Tv.P(:,ti)+(1:TL+1)''/1e10,Tv.' bvars{vi} '(:,ti),Pe,''linear'');']);
            eval(['ZvP.' bvars{vi} '(1,ti) = 0;']);eval(['ZvP.' bvars{vi} '(end,ti) = Zv.' bvars{vi} '(end,ti);']);
            eval(['TvP.' bvars{vi} '(1,ti) = 0;']);eval(['TvP.' bvars{vi} '(end,ti) = Tv.' bvars{vi} '(1,ti);']);
            eval(['YvP.' bvars{vi} '(1,ti) = 0;']);eval(['YvP.' bvars{vi} '(end,ti) = Yv.' bvars{vi} '(end,ti);']);
        end
    end

% $$$     %%%% Check/plot budgets:
% $$$     colors = {'m','b','k','r',[0 0.5 0]};     
% $$$     figure;
% $$$     set(gcf,'Position',[1921           1        1920        1005]);
% $$$     subplot(1,3,1);
% $$$     bvars = {'TEN_c','ADV_c','FOR_c','RMIX_c','VMIX_c'};
% $$$     for gi=1:length(bvars)
% $$$         eval(['var = mean(ZvP.' bvars{gi} ',2);']);
% $$$         plot(var/1e15,Pe,'-','color',colors{gi},'linewidth',2);
% $$$         hold on;
% $$$     end
% $$$     plot([0 0],[0 100],'--k');
% $$$     ylabel('Percentile');
% $$$     set(gca,'ydir','reverse');
% $$$     ylim([0 100]);
% $$$     xlim([-2 2]);
% $$$     xlabel('Upward Vertical heat transport (PW)');
% $$$     legend('$\partial\mathcal{H}_z/\partial t$','$\mathcal{A}_z$','$\mathcal{F}_z$',...
% $$$            '$\mathcal{M}_z^{neutral}$',['$\' ...
% $$$                         'mathcal{M}_z^{vertical}$']);
% $$$     set(gca,'Position',[0.06   0.1400    0.2580    0.8150]);
% $$$     
% $$$     subplot(1,3,2);
% $$$     bvars = {'TEN_c','ADV_c','FOR_c','RMIX_c','VMIX_c'};
% $$$     for gi=1:length(bvars)
% $$$         eval(['var = mean(TvP.' bvars{gi} ',2);']);
% $$$         plot(var/1e15,Pe,'-','color',colors{gi},'linewidth',2);
% $$$         hold on;
% $$$     end
% $$$     plot([0 0],[0 100],'--k');
% $$$     set(gca,'ydir','reverse');
% $$$     ylim([0 100]);
% $$$     set(gca,'yticklabel',[]);
% $$$     xlabel('Cold-to-warm Diathermal heat transport (PW)');
% $$$     legend('$\partial\mathcal{H}_\Theta/\partial t$','$\mathcal{M}_\Theta^{numerical}$','$\mathcal{F}_\Theta$',...
% $$$            '$\mathcal{M}_\Theta^{neutral}$','$\mathcal{M}_\Theta^{vertical}$');
% $$$     set(gca,'Position',[0.3661    0.1400    0.2580    0.8150]);
% $$$ 
% $$$     subplot(1,3,3);
% $$$     bvars = {'TEN_c','ADV_c','FOR_c','RMIX_c'};
% $$$     for gi=1:length(bvars)
% $$$         eval(['var = mean(YvP.' bvars{gi} ',2);']);
% $$$         plot(var/1e15,Pe,'-','color',colors{gi},'linewidth',2);
% $$$         hold on;
% $$$     end
% $$$     plot([0 0],[0 100],'--k');
% $$$     ylim([0 100]);
% $$$     xlabel('Southward Meridional heat transport (PW)');
% $$$     legend('$\partial\mathcal{H}_\phi/\partial t$','$\mathcal{A}_\phi^{advective}$','$\mathcal{F}_\phi$',...
% $$$            '$\mathcal{A}_\phi^{diffusive}$');
% $$$     set(gca,'Position',[0.7    0.1400    0.2580    0.8150]);

    % Time-integrate budget terms:
    for ti =1:length(typs)
        for vi=1:length(bvars)
            eval([typs{ti} 'vP.' bvars{vi} ' = cumsum(' typs{ti} ...
                  'vP.' bvars{vi} '.*repmat(DT_A'',[PL+1 1]),2);']);
        end
    end    

end

    zofP_mean = -interp1(mean(Zv.P,2),Ze,P,'linear'); % Depth axis to
                                                    % go with depth-volume
    yofP_mean = interp1(mean(Yv.P,2)+(1:yL+1)'/1e10,latv_edges,P,'linear'); % Depth axis to
                                                    % go with depth-volume
    
    % Construct climatology and make years even:
    if (PI_or_his == 1)
        yr1 =50; % initial years to throw out
% $$$         yr1 =0;
    else
        yr1 = 0;
    end
    nyrs = floor(tL/12)-yr1; % number of years
    tL = nyrs*12; % number of months
    ti = yr1*12+1; % starting month
    
    for ty = 1:length(typs)
        eval(['fnames = fieldnames(' typs{ty} 'vP);']);
        for vi =1:length(fnames)
            eval([typs{ty} 'vP.' fnames{vi} ' = ' typs{ty} 'vP.' fnames{vi} '(:,ti:(ti+tL-1));']);
        end
        eval(['fnames = fieldnames(' typs{ty} 'v);']);
        for vi =1:length(fnames)
            eval([typs{ty} 'v.' fnames{vi} ' = ' typs{ty} 'v.' fnames{vi} '(:,ti:(ti+tL-1));']);
        end
    end
    TIMESERIES = {'CIN.N34','CIN.TPIr1','CIN.TPIr2','CIN.TPIr3','OHC','CIN.AMOC','CIN.WPOW'};
    for vi=1:length(TIMESERIES)
        eval([TIMESERIES{vi} ' = ' TIMESERIES{vi} '(ti:(ti+tL-1));']);
    end
    time = time(ti:(ti+tL-1));
    DT_A = DT_A(ti:(ti+tL-1));
    Vtot = Vtot(ti:(ti+tL-1));
    CIN.AMOCfull = CIN.AMOC;
    OHCfull = OHC;
 
    %%% De-drifing:
    
    % Calculate cubic drift trend and subtract:
    fnames = fieldnames(TvP);
    if (isfield(TvP,'Tap'))
        ZvP.Tap = TvP.Tap;
        YvP.Tap = TvP.Tap;
    end
    if (isfield(TvP,'Hap'))
        ZvP.Hap = TvP.Hap;
        YvP.Hap = TvP.Hap; % for convience -> delete after...
    end
    if (PI_or_his == 1) % Calculate cubic drift from PI-control
        
    for ti = 1:length(typs)
        for vi = 1:length(fnames)
            tn = typs{ti};
            vn = fnames{vi};

            % Mean:
            eval([tn 'vP.' vn '_mean = mean(' tn 'vP.' vn ',2);']);
            
            % Calculate cubic trend:
            eval(['lg = length(' tn 'vP.' vn '(:,1));']);
            eval([tn 'vP.' vn '_cub = zeros(tL,lg);']);
            eval([tn 'vP.' vn '_cubtr = zeros(4,lg);']);
            for pi=1:lg
                eval(['[' tn 'vP.' vn '_cub(:,pi),' tn 'vP.' vn ...
                      '_cubtr(:,pi)] = cubfit(time,' tn 'vP.' vn ...
                      '(pi,:)'');']);
            end
            eval([tn 'vP.' vn '_cub = ' tn 'vP.' vn '_cub'';']);
            
            % subtract cubic trend:
            eval([tn 'vP.' vn ' = ' tn 'vP.' vn '-' tn 'vP.' vn '_cub;']);
            
        end
    end
    [OHC_cub,OHC_cubtr] = cubfit(time,OHC);
    OHC = OHC-OHC_cub;

    else 

        % Use PI-control cubic drift for de-drifting
        
        Itime = -100;
        t = time-Itime;
        t2 = (time-Itime).^2;
        t3 = (time-Itime).^3;
        
        PIcont = load([baseMAT 'PIcontrolPP_Tb05.mat'],'ZvP','TvP','YvP','OHC_cubtr','time');
        
        for ti = 1:length(typs)
            for vi = 1:length(fnames)
                tn = typs{ti};
                vn = fnames{vi};

                % Mean (from PI control):
                eval([tn 'vP.' vn '_mean = PIcont.' tn 'vP.' vn '_mean;']);
            
                % Calculate cubic trend:
                eval(['lg = length(' tn 'vP.' vn '(:,1));']);
                eval([tn 'vP.' vn '_cub = zeros(tL,lg);']);
                for pi=1:lg
                    eval([tn 'vP.' vn '_cub(:,pi) = PIcont.' ...
                          tn 'vP.' vn '_cubtr(1,pi) + PIcont.' ...
                          tn 'vP.' vn '_cubtr(2,pi).*t + PIcont.' ...
                          tn 'vP.' vn '_cubtr(3,pi).*t2 + PIcont.' ...
                          tn 'vP.' vn '_cubtr(4,pi).*t3;']);
                end
                eval([tn 'vP.' vn '_cub = ' tn 'vP.' vn '_cub'';']);
            
                % subtract cubic trend:
                eval([tn 'vP.' vn ' = ' tn 'vP.' vn '-' tn 'vP.' vn '_cub;']);
            
            end
        end
        OHC_cub = PIcont.OHC_cubtr(1)+PIcont.OHC_cubtr(2).*t + ...
                  PIcont.OHC_cubtr(3).*t2 + PIcont.OHC_cubtr(4).*t3;
        OHC = OHC-OHC_cub;
        OHCcubdd = OHC;
        TIMESERIES = {TIMESERIES{:},'OHCcubdd'};

        % Set date:
        time = time+1850;
        yr1 = 1850;
        
        % Dedrift with linear fit to historical e1:
        if (PI_or_his == 2)
            for ti = 1:length(typs)
                for vi = 1:length(fnames)
                    tn = typs{ti};
                    vn = fnames{vi};

                    % Calculate linear trend:
                    eval(['lg = length(' tn 'vP.' vn '(:,1));']);
                    eval([tn 'vP.' vn '_lin = zeros(tL,lg);']);
                    eval([tn 'vP.' vn '_lintr = zeros(2,lg);']);
                    for pi=1:lg
                        eval(['[' tn 'vP.' vn '_lin(:,pi),' tn 'vP.' vn ...
                              '_lintr(:,pi)] = linfit(time,' tn 'vP.' vn ...
                              '(pi,:)'');']);
                    end
                    eval([tn 'vP.' vn '_lin = ' tn 'vP.' vn '_lin'';']);
            
                    % subtract linear trend:
                    eval([tn 'vP.' vn ' = ' tn 'vP.' vn '-' tn 'vP.' vn '_lin;']);
            
                end
            end
            [OHC_lin,OHC_lintr] = linfit(time,OHC);
            OHC = OHC-OHC_lin;
        elseif (PI_or_his == 0)
            hisNATe1 = load([baseMAT 'hisNATe1PP_Tb05.mat'],'ZvP','TvP','YvP','OHC_lintr','time');
        
            for ti = 1:length(typs)
                for vi = 1:length(fnames)
                    tn = typs{ti};
                    vn = fnames{vi};

                    % Calculate linear trend:
                    eval(['lg = length(' tn 'vP.' vn '(:,1));']);
                    eval([tn 'vP.' vn '_lin = zeros(tL,lg);']);
                    for pi=1:lg
                        eval([tn 'vP.' vn '_lin(:,pi) = hisNATe1.' ...
                              tn 'vP.' vn '_lintr(1,pi) + hisNATe1.' ...
                              tn 'vP.' vn '_lintr(2,pi).*time;']);
                    end
                    eval([tn 'vP.' vn '_lin = ' tn 'vP.' vn '_lin'';']);
            
                    % subtract cubic trend:
                    eval([tn 'vP.' vn ' = ' tn 'vP.' vn '-' tn 'vP.' vn '_lin;']);
            
                end
            end
            OHC_lin = hisNATe1.OHC_lintr(1)+hisNATe1.OHC_lintr(2).*time;
            OHC = OHC-OHC_lin;
        end
    end
        
    for vi=1:length(TIMESERIES)
        eval([TIMESERIES{vi} '_clim = zeros(12,1);']);
        for mi=1:12
            eval([TIMESERIES{vi} '_clim(mi) = mean(' TIMESERIES{vi} ...
                  '(mi:12:end));']);
        end
        eval([TIMESERIES{vi} '_clim = ' TIMESERIES{vi} '_clim - mean(' ...
              TIMESERIES{vi} '_clim);']);
        
        eval([TIMESERIES{vi} ' = ' TIMESERIES{vi} ...
              ' - repmat(' TIMESERIES{vi} '_clim,[nyrs 1]);']);
    end

    for ti = 1:length(typs)
        for vi = 1:length(fnames)
            tn = typs{ti};
            vn = fnames{vi};

            % Calculate Seasonal cycle:
            eval(['lg = length(' tn 'vP.' vn '(:,1));']);
            eval([tn 'vP.' vn '_clim = zeros(lg,12);']);
            for mi=1:12
                eval([tn 'vP.' vn '_clim(:,mi) = mean(' tn 'vP.' vn ...
                      '(:,mi:12:end),2);']);
            end
            
            % Subtract mean of season to avoid changing mean:
            eval([tn 'vP.' vn '_clim = ' tn 'vP.' vn '_clim-repmat(mean(' ...
                  tn 'vP.' vn '_clim,2),[1 12]);']);
            
            % Subtract seasonal cycle:
            eval([tn 'vP.' vn ' = ' tn 'vP.' vn ' - repmat(' tn 'vP.' ...
                  vn '_clim,[1 nyrs]);']);
        end
    end
    
    Vtot_clim = zeros(12,1);
    for mi=1:12
        Vtot_clim(mi) = mean(Vtot(mi:12:end));
    end
    
    % Calculate cumulative integral OHC (overwrite total - just anomalies now):
    TvP.Hp = rho0*Cp*cumsum(TvP.Tp,1)*dP/100.*repmat(Vtot',[PL 1]);
    ZvP.Hp = rho0*Cp*cumsum(ZvP.Tp,1)*dP/100.*repmat(Vtot',[PL 1]);
    YvP.Hp = rho0*Cp*cumsum(YvP.Tp,1)*dP/100.*repmat(Vtot',[PL 1]);
    
    % Calculate tendency term contributions to temperature:
    if (doBUDGET)
    for ti =1:length(typs)
        for vi = 1:length(bvars)
            tn = typs{ti};
            vn = bvars{vi};
            eval([tn 'vP.' vn(1:end-2) 'p = diff(' tn 'vP.' vn ',[],1)/rho0/' ...
                                'Cp/dP*100./repmat(Vtot'',[PL ' ...
                  '1]);']);
            eval([tn 'vP.' vn(1:end-2) 'p_clim = diff(' tn 'vP.' vn '_clim,[],1)/rho0/' ...
                                'Cp/dP*100./repmat(Vtot_clim'',[PL ' ...
                  '1]);']);
        end
    end
    end
    
    % Calc GMSST:
    CIN.GMSST = ZvP.Tp(1,:)';
    
    % Calc TPI:
    CIN.TPIraw = CIN.TPIr2-(CIN.TPIr1+CIN.TPIr3)/2;
    [b,a] = cheby1(6,10,2/13/12);
    CIN.TPI = filter(b,a,CIN.TPIraw);
    % Shift time (not sure why I have to do this???):
    CIN.TPI(1:(end-13*12)) = CIN.TPI(13*12:end-1);
    CIN.TPI((end-13*12+1):end) = 0;
    
    % Add in mean profiles for %-transformations:
    ZvP.mean = zofP_mean;
    TvP.mean = TvP.Tp_mean;
    YvP.mean = yofP_mean;

    % Tyz processing:
    if (doTyz)
        
% $$$         % Pre-process Tyz and Hyz:
% $$$         load([baseMAT 'CM2_PIcontrol__Tyz_ALL.mat']);
% $$$         Tyz = TyzS.H/TyzS.V/rho0/Cp;
% $$$         save([baseMAT 'CM2_PIcontrol__Tyz_ALL_TyzONLY.mat'],'Tyz');
% $$$         % Time-mean for plotting:
% $$$         Tyz_mean = mean(Tyz,3);
% $$$ 
% $$$         load([baseMAT 'CM2_PIcontrol__Tyz_ALL.mat']);
% $$$         Ayz = ndgrid(diff(latv_edges),diff(Ze));
% $$$         Tyz = TyzS.H./repmat(Ayz,[1 1 length(TyzS.H(1,1,:))]);
% $$$         save([baseMAT 'CM2_PIcontrol__Tyz_ALL_HyzONLY.mat'],'Tyz');
                
        % Load:
        load([baseMAT 'CM2_PIcontrol__Tyz_ALL_HyzONLY.mat']);
        ti = yr1*12+1; % starting month
        Tyz = Tyz(:,:,ti:(ti+tL-1));     
        
        % Detrend:
        [Tyz_cub,Tyz_cubtr] = cubfit(time,reshape(Tyz,[yL*zL tL])');
        Tyz_cub = reshape(Tyz_cub',[yL zL tL]);
        Tyz_cubtr = reshape(Tyz_cubtr',[yL zL 4]);
        Tyz = Tyz-Tyz_cub;
        
        % Deseason:
        Tyz_clim = zeros(yL,zL,12);
        for mi=1:12
            Tyz_clim(:,:,mi) = mean(Tyz(:,:,mi:12:end),3);
        end
        Tyz = Tyz - repmat(Tyz_clim,[1 1 nyrs]);
        
    end
    
    
if (saveMAT)
    save([baseMAT saveNAME]);

end

