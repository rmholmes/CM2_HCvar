% Post-processing for HC variability project

clear all;
PI_or_his = 1; % 1 = PI-control, 0 = historical simualtion
if (PI_or_his)
    load('CM2_PIcontrol_ALL.mat');
else
    load('CM2_historical_ALL.mat');
end

% Define a new percentile grid:
dP = 0.25;
Pe = 0:dP:100;
P = (Pe(2:end)+Pe(1:end-1))/2;
PL = length(P);

% Time vector:
tL = length(time);
time = time/365.25; % time in years

% Fix advection budget term by residual:
typs = {'T','Z','Y'};
vars = {'temp_submeso', 'temp_vdiffuse_diff_cbt', 'temp_nonlocal_KPP', ...
        'temp_vdiffuse_sbc','frazil_3d','sw_heat','temp_rivermix', ...
        'neutral_diffusion_temp','neutral_gm_temp', 'temp_vdiffuse_k33', 'mixdownslope_temp', ...
        'temp_sigma_diff','sfc_hflux_pme','temp_eta_smooth'};    
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
    eval([typs{vi} 'v.FOR = ' typs{vi} 'v.temp_vdiffuse_sbc+' typs{vi} 'v.frazil_3d+' typs{vi} 'v.sw_heat' ...
          '+' typs{vi} 'v.temp_rivermix+' typs{vi} 'v.sfc_hflux_pme' ...
          '+' typs{vi} 'v.temp_eta_smooth;']);
    eval([typs{vi} 'v.RMIX = ' typs{vi} 'v.neutral_diffusion_temp+' typs{vi} 'v.temp_vdiffuse_k33+' typs{vi} 'v.mixdownslope_temp' ...
          '+' typs{vi} 'v.temp_sigma_diff;']);
    eval([typs{vi} 'v.VMIX = ' typs{vi} 'v.temp_vdiffuse_diff_cbt+' typs{vi} 'v.temp_nonlocal_KPP;']);
end
vars = {vars{:},'temp_tendency','temp_advection'};
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

% $$$ %%% New approach: Interpolate HC
% $$$     % Interpolate only high-resolution grids for smoothness:
% $$$     zfac = 1;
% $$$ % $$$     Zv.H_c = interp1(Zv.H_c,1:1/zfac:(zL+1));Zv.V_c = interp1(Zv.V_c,1:1/zfac:(zL+1));
% $$$     Tfac = 1;
% $$$ % $$$     Tv.H_c = interp1(Tv.H_c,1:1/Tfac:(TL+1));Tv.V_c = interp1(Tv.V_c,1:1/Tfac:(TL+1));
% $$$     yfac = 1;
% $$$ % $$$     Yv.H_c = interp1(Yv.H_c,1:1/yfac:(yL+1));Yv.V_c = interp1(Yv.V_c,1:1/yfac:(yL+1));
% $$$ 
% $$$     % Remap to ocean percentile
% $$$     Vtot = Zv.V_c(end,:)';
% $$$     Tv.P   = 100*Tv.V_c./repmat(Vtot',[TL*Tfac+1 1]);
% $$$     Zv.P   = 100*Zv.V_c./repmat(Vtot',[zL*zfac+1 1]);
% $$$     Yv.P   = 100*Yv.V_c./repmat(Vtot',[yL*yfac+1 1]);
% $$$ 
% $$$     ZvP.Hp = zeros(PL+1,tL);
% $$$     YvP.Hp = zeros(PL+1,tL);
% $$$     TvP.Hp = zeros(PL+1,tL);
% $$$     for ti = 1:tL
% $$$         ZvP.Hp(:,ti) = interp1((Zv.P(:,ti)+(1:zL*zfac+1)'/1e10),Zv.H_c(:,ti),Pe,'linear');
% $$$         YvP.Hp(:,ti) = interp1((Yv.P(:,ti)+(1:yL*yfac+1)'/1e10),Yv.H_c(:,ti),Pe,'linear');
% $$$         TvP.Hp(:,ti) = interp1((Tv.P(:,ti)+(1:TL*Tfac+1)'/1e10),Tv.H_c(:,ti),Pe,'linear');
% $$$     end
% $$$     ZvP.Hp(1,:) = 0;
% $$$     YvP.Hp(1,:) = 0;
% $$$     TvP.Hp(1,:) = 0;
% $$$ 
% $$$     % Calculate temperatures from HC:
% $$$     ZvP.Tp = (ZvP.Hp(2:end,:)-ZvP.Hp(1:(end-1),:))/rho0/Cp*100/dP./repmat(Vtot',[PL 1]);
% $$$     YvP.Tp = (YvP.Hp(2:end,:)-YvP.Hp(1:(end-1),:))/rho0/Cp*100/dP./repmat(Vtot',[PL 1]);
% $$$     TvP.Tp = (TvP.Hp(2:end,:)-TvP.Hp(1:(end-1),:))/rho0/Cp*100/dP./repmat(Vtot',[PL 1]);
% $$$ 
% $$$     zofP_mean = -interp1(mean(Zv.P(1:zfac:end,:),2),Ze,P,'linear'); % Depth axis to
% $$$                                                     % go with depth-volume
% $$$     yofP_mean = interp1(mean(Yv.P(1:yfac:end,:),2)+(1:yL+1)'/1e10,latv_edges,P,'linear'); % Depth axis to
% $$$                                                     % go with depth-volume

    %%% Old approach: Interpolate temperatures
    % Remap to ocean percentile
    Vtot = Zv.V_c(end,:)';
    Atot = Tv.A_c(1,:)';
    Tv.P   = 100*Tv.V_c./repmat(Vtot',[TL+1 1]);
    Zv.P   = 100*Zv.V_c./repmat(Vtot',[zL+1 1]);
    Yv.P   = 100*Yv.V_c./repmat(Vtot',[yL+1 1]);
    Tv.Pa  = 100*Tv.A_c./repmat(Atot',[TL+1 1]);

    % Calculate z-temperature:
    Zv.T = Zv.H/rho0/Cp./Zv.V;
    Zv.Te = cat(1,Zv.T(1,:),(Zv.T(2:end,:)+Zv.T(1:end-1,:))/2,Zv.T(end,:));

    % Calculate y-temperature:
    Yv.T = Yv.H/rho0/Cp./Yv.V;
    Yv.Te = cat(1,Yv.T(1,:),(Yv.T(2:end,:)+Yv.T(1:end-1,:))/2,Yv.T(end,:));

    % Interpolate variables to P levels:
    TvP.Tp = zeros(PL,tL);
    TvP.Tap = zeros(PL,tL);
    ZvP.Tp = zeros(PL,tL);
    YvP.Tp = zeros(PL,tL);
    bvars = {'TEN_c','ADV_c','FOR_c','RMIX_c','VMIX_c'};
    for ti = 1:tL
        TvP.Tp(:,ti) = interp1(Tv.P(:,ti)+(1:TL+1)'/1e10,Te,P,'linear');
        TvP.Tp(1,ti) = Te(end);
        TvP.Tap(:,ti) = interp1(Tv.Pa(:,ti)+(1:TL+1)'/1e10,Te,P,'linear');
        TvP.Tap(1,ti) = Te(end);
        ZvP.Tp(:,ti) = interp1(Zv.P(:,ti)+(1:zL+1)'/1e10,Zv.Te(:,ti),P,'linear');
        YvP.Tp(:,ti) = interp1(Yv.P(:,ti)+(1:yL+1)'/1e10,Yv.Te(:,ti),P,'linear');
        for vi=1:length(bvars)
            eval(['ZvP.' bvars{vi} '(:,ti) = interp1(Zv.P(:,ti)+(1:zL+1)''/1e10,Zv.' bvars{vi} '(:,ti),Pe,''linear'');']);
            eval(['YvP.' bvars{vi} '(:,ti) = interp1(Yv.P(:,ti)+(1:yL+1)''/1e10,Yv.' bvars{vi} '(:,ti),Pe,''linear'');']);
            eval(['TvP.' bvars{vi} '(:,ti) = interp1(Tv.P(:,ti)+(1:TL+1)''/1e10,Tv.' bvars{vi} '(:,ti),Pe,''linear'');']);
            eval(['ZvP.' bvars{vi} '(1,ti) = 0;']);eval(['ZvP.' bvars{vi} '(end,ti) = Zv.' bvars{vi} '(end,ti);']);
            eval(['TvP.' bvars{vi} '(1,ti) = 0;']);eval(['TvP.' bvars{vi} '(end,ti) = Tv.' bvars{vi} '(1,ti);']);
            eval(['YvP.' bvars{vi} '(1,ti) = 0;']);eval(['YvP.' bvars{vi} '(end,ti) = Yv.' bvars{vi} '(end,ti);']);
        end
    end
    
% $$$     %%%% Check budgets:
% $$$     colors = {'m','b','k','r',[0 0.5 0]};     
% $$$ % $$$     figure;
% $$$ % $$$     set(gcf,'Position',[1921           1        1920        1005]);
% $$$ % $$$     subplot(1,3,1);
% $$$ % $$$     for gi=1:length(bvars)
% $$$ % $$$         eval(['var = mean(ZvP.' bvars{gi} ',2);']);
% $$$ % $$$         plot(var/1e15,Pe,'-','color',colors{gi},'linewidth',2);
% $$$ % $$$         hold on;
% $$$ % $$$     end
% $$$ % $$$     plot([0 0],[0 100],'--k');
% $$$ % $$$     ylabel('Depth percentile');
% $$$ % $$$     set(gca,'ydir','reverse');
% $$$ % $$$     ylim([0 100]);
% $$$ % $$$     xlim([-2 2]);
% $$$ % $$$     xlabel('Vertical heat transport (PW)');
% $$$ % $$$     
% $$$     subplot(1,3,2);
% $$$     for gi=1:length(bvars)
% $$$         eval(['var = mean(TvP.' bvars{gi} ',2);']);
% $$$         plot(var/1e15,Pe,'-','color',colors{gi},'linewidth',2);
% $$$         hold on;
% $$$     end
% $$$     legend('Tendency','Advection','Surface Forcing','Neutral Mixing','Vertical Mixing');
% $$$     plot([0 0],[0 100],'--k');
% $$$     set(gca,'ydir','reverse');
% $$$     ylim([0 100]);
% $$$     ylabel('Temperature percentile');
% $$$     xlabel('Diathermal heat transport (PW)');
% $$$ 
% $$$     subplot(1,3,3);
% $$$     for gi=1:length(bvars)
% $$$         eval(['var = mean(YvP.' bvars{gi} ',2);']);
% $$$         plot(var/1e15,Pe,'-','color',colors{gi},'linewidth',2);
% $$$         hold on;
% $$$     end
% $$$     plot([0 0],[0 100],'--k');
% $$$     ylim([0 100]);
% $$$     ylabel('Latitude percentile');
% $$$     xlabel('Meridional heat transport (PW)');


    zofP_mean = -interp1(mean(Zv.P,2),Ze,P,'linear'); % Depth axis to
                                                    % go with depth-volume
    yofP_mean = interp1(mean(Yv.P,2)+(1:yL+1)'/1e10,latv_edges,P,'linear'); % Depth axis to
                                                    % go with depth-volume

    % Time-integrate budget terms:
    for ti =1:length(typs)
        for vi=1:length(bvars)
            eval([typs{ti} 'vP.' bvars{vi} ' = cumsum(' typs{ti} ...
                  'vP.' bvars{vi} '.*repmat(DT_A'',[PL+1 1]),2);']);
        end
    end    
    
    % Construct climatology and make years even:
    if (PI_or_his)
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
    Vtot = Vtot(ti:(ti+tL-1));
    CIN.AMOCfull = CIN.AMOC;
    OHCfull = OHC;
        
    % Calculate cubic drift trend and subtract:
    fnames = fieldnames(TvP);
    ZvP.Tap = TvP.Tap;
    YvP.Tap = TvP.Tap; % for convience -> delete after...
    if (PI_or_his)
        
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
                eval(['[' tn 'vP.' vn '_cub(:,pi),' tn 'vP.' tn ...
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
    
% $$$     PI.ZvP.Tp_cubtr = ZvP.Tp_cubtr;
% $$$     PI.TvP.Tp_cubtr = TvP.Tp_cubtr;
% $$$     PI.YvP.Tp_cubtr = YvP.Tp_cubtr;
% $$$     PI.time = time;
% $$$     PI.OHCfull = OHCfull;
% $$$     PI.OHC = OHC;
% $$$     PI.OHC_cubtr = OHC_cubtr;
% $$$     PI.TvP.Tp_mean = TvP.Tp_mean;
% $$$     PI.YvP.Tp_mean = YvP.Tp_mean;
% $$$     PI.ZvP.Tp_mean = ZvP.Tp_mean;
% $$$     save('PIcubicFit.mat','PI');

    else
        ' NOT ADJUSTED YET!'
% $$$         Itime = -100;
% $$$         load('PIcubicFit.mat');
% $$$         TvP.Tp_mean = PI.TvP.Tp_mean;
% $$$         YvP.Tp_mean = PI.YvP.Tp_mean;
% $$$         ZvP.Tp_mean = PI.ZvP.Tp_mean;
% $$$ 
% $$$         TvP.Tp_cub = zeros(tL,PL);
% $$$         TvP.Tap_cub = zeros(tL,PL);
% $$$         ZvP.Tp_cub = zeros(tL,PL);
% $$$         YvP.Tp_cub = zeros(tL,PL);
% $$$         
% $$$         t = time-Itime;
% $$$         t2 = (time-Itime).^2;
% $$$         t3 = (time-Itime).^3;
% $$$         
% $$$         for pi = 1:PL
% $$$             TvP.Tp_cub(:,pi) = PI.TvP.Tp_cubtr(1,pi)+PI.TvP.Tp_cubtr(2,pi).*t ...
% $$$                           + PI.TvP.Tp_cubtr(3,pi).*t2 + PI.TvP.Tp_cubtr(4,pi).*t3;
% $$$             YvP.Tp_cub(:,pi) = PI.YvP.Tp_cubtr(1,pi)+PI.YvP.Tp_cubtr(2,pi).*t ...
% $$$                           + PI.YvP.Tp_cubtr(3,pi).*t2 + PI.YvP.Tp_cubtr(4,pi).*t3;
% $$$             ZvP.Tp_cub(:,pi) = PI.ZvP.Tp_cubtr(1,pi)+PI.ZvP.Tp_cubtr(2,pi).*t ...
% $$$                           + PI.ZvP.Tp_cubtr(3,pi).*t2 + PI.ZvP.Tp_cubtr(4,pi).*t3;
% $$$         end
% $$$         OHC_cub = PI.OHC_cubtr(1)+PI.OHC_cubtr(2).*t + ...
% $$$                   PI.OHC_cubtr(3).*t2 + PI.OHC_cubtr(4).*t3;
% $$$     
% $$$         % subtract cubear trend:
% $$$         TvP.Tp = TvP.Tp-TvP.Tp_cub';
% $$$         ZvP.Tp = ZvP.Tp-ZvP.Tp_cub';
% $$$         YvP.Tp = YvP.Tp-YvP.Tp_cub';
% $$$         OHC = OHC-OHC_cub;
% $$$         
% $$$         % Set date:
% $$$         time = time+1850;
% $$$         yr1 = 1850;
    end
        
    for vi=1:length(TIMESERIES)
        eval([TIMESERIES{vi} '_clim = zeros(12,1);']);
        for mi=1:12
            eval([TIMESERIES{vi} '_clim(mi) = mean(' TIMESERIES{vi} ...
                  '(mi:12:end));']);
        end
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
        
% $$$     if (PI_or_his)
% $$$         PI.OHC = OHC;
% $$$         save('PIcubicFit.mat','PI');
% $$$     end
    
    % Calculate cumulative integral OHC (overwrite total - just anomalies now):
    TvP.Hp = rho0*Cp*cumsum(TvP.Tp,1)*dP/100.*repmat(Vtot',[PL 1]);
    ZvP.Hp = rho0*Cp*cumsum(ZvP.Tp,1)*dP/100.*repmat(Vtot',[PL 1]);
    YvP.Hp = rho0*Cp*cumsum(YvP.Tp,1)*dP/100.*repmat(Vtot',[PL 1]);
    
    % Calculate tendency term contributions to temperature:
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
    
% $$$     % Calculate quantities for Time Of Emergence:
% $$$     if (PI_or_his)
% $$$     PIstd.stdZvP.Tp = std(ZvP.Tp,[],2);
% $$$     PIstd.stdTvP.Tp = std(TvP.Tp,[],2);
% $$$     PIstd.stdYvP.Tp = std(YvP.Tp,[],2);
% $$$     PIstd.stdOHC = std(OHC,[],1);
% $$$     save('PIstd.mat','PIstd');
% $$$     else
% $$$         load('PIstd.mat');
% $$$     end

