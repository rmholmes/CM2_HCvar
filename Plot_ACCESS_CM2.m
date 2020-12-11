% Plotting for HC variability project

clear all;
PI_or_his = 1; % 1 = PI-control, 0 = historical simualtion
if (PI_or_his)
    load('CM2_PIcontrol_ALL.mat');
else
    load('CM2_historical_ALL.mat');
end

%%% POST PROCESSING

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
            eval(['ZvP.' bvars{vi} '(:,ti) = interp1(Zv.P(:,ti)+(1:zL+1)''/1e10,Zv.' bvars{vi} '(:,ti),P,''linear'');']);
            eval(['YvP.' bvars{vi} '(:,ti) = interp1(Yv.P(:,ti)+(1:yL+1)''/1e10,Yv.' bvars{vi} '(:,ti),P,''linear'');']);
            eval(['TvP.' bvars{vi} '(:,ti) = interp1(Tv.P(:,ti)+(1:TL+1)''/1e10,Tv.' bvars{vi} '(:,ti),P,''linear'');']);
        end
    end
    
% $$$     %% Check budgets:
% $$$     colors = {'m','b','k','r',[0 0.5 0]};     
% $$$     figure;
% $$$     subplot(1,3,1);
% $$$     for gi=1:length(bvars)
% $$$         eval(['var = mean(ZvP.' bvars{gi} ',2);']);
% $$$         plot(var/1e15,P,'-','color',colors{gi});
% $$$         hold on;
% $$$     end
% $$$     ylabel('Depth percentile');
% $$$     set(gca,'ydir','reverse');
% $$$     ylim([0 100]);
% $$$     xlim([-2 2]);
% $$$     xlabel('Vertical heat transport (PW)');
% $$$     
% $$$     subplot(1,3,2);
% $$$     for gi=1:length(bvars)
% $$$         eval(['var = mean(TvP.' bvars{gi} ',2);']);
% $$$         plot(var/1e15,P,'-','color',colors{gi});
% $$$         hold on;
% $$$     end
% $$$     legend('Tendency','Advection','Surface Forcing','Neutral Mixing','Vertical Mixing');
% $$$     set(gca,'ydir','reverse');
% $$$     ylim([0 100]);
% $$$     ylabel('Temperature percentile');
% $$$     xlabel('Diathermal heat transport (PW)');
% $$$ 
% $$$     subplot(1,3,3);
% $$$     for gi=1:length(bvars)
% $$$         eval(['var = mean(YvP.' bvars{gi} ',2);']);
% $$$         plot(var/1e15,P,'-','color',colors{gi});
% $$$         hold on;
% $$$     end
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
                  'vP.' bvars{vi} '.*repmat(DT_A'',[PL 1]),2);']);
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
% $$$     TvP.Tp = TvP.Tp(:,ti:(ti+tL-1));
% $$$     TvP.Tap = TvP.Tap(:,ti:(ti+tL-1));
% $$$     ZvP.Tp = ZvP.Tp(:,ti:(ti+tL-1)); 
% $$$     YvP.Tp = YvP.Tp(:,ti:(ti+tL-1)); 
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
            eval([tn 'vP.' vn '_cub = zeros(tL,PL);']);
            eval([tn 'vP.' vn '_cubtr = zeros(4,PL);']);
            for pi=1:PL
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
            eval([tn 'vP.' vn '_clim = zeros(PL,12);']);
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
                                'Cp/dP*100./repmat(Vtot'',[PL-1 ' ...
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

% $$$     %%%% Plot total heat content and fits time series:
% $$$     figure;
% $$$     set(gcf,'Position',[1 41 2560 1327.3]);
% $$$     subplot(5,1,1);
% $$$     plot(time,filter_field(CIN.N34,5,'-t'));
% $$$     hold on;
% $$$     plot(time,CIN.GMSST*10,'-r');
% $$$     legend('Nino 3.4','Global Mean SST Anomalies * 10');
% $$$     xlabel('Year');
% $$$     ylabel('SST ($^\circ$C)');
% $$$     xlim([yr1 yr1+nyrs]);
% $$$ 
% $$$     subplot(5,1,2);
% $$$     plot(time,OHC,'-k');
% $$$     hold on;
% $$$     plot(time,OHCfull,'-b');
% $$$     plot(time,OHC_cub,'-r');
% $$$     legend('OHC dedrifted','OHC deseasoned','OHC cubic fix');
% $$$     xlabel('Year');
% $$$     ylabel('Ocean Heat Content (J)');
% $$$     xlim([yr1 yr1+nyrs]);
% $$$ 
% $$$     subplot(5,1,3);
% $$$     plot(time,OHC,'-k');
% $$$     hold on;
% $$$     plot(time,TvP.Hp(end,:),'--g');
% $$$     plot(time,ZvP.Hp(end,:),'--c');
% $$$     plot(time,YvP.Hp(end,:),'--y');
% $$$     legend('OHC dedrifted',...
% $$$            'OHC from cumsum(TvP.Tp)','OHC from cumsum(ZvP.Tp)','OHC from cumsum(YvP.Tp)');
% $$$     xlabel('Year');
% $$$     ylabel('Ocean Heat Content (J)');
% $$$     xlim([yr1 yr1+nyrs]);
% $$$ 
% $$$     subplot(5,1,4);
% $$$     plot(time,CIN.AMOCfull,'-k');
% $$$     hold on;
% $$$     plot(time,filter_field(CIN.AMOCfull,5*12+1,'-t'),'-k','linewidth',2);
% $$$     xlabel('Year');
% $$$     ylabel('AMOC at $26^\circ$N (Sv)');
% $$$     xlim([yr1 yr1+nyrs]);
% $$$ 
% $$$     subplot(5,1,5);
% $$$     plot(time,CIN.WPOW,'-k');
% $$$     hold on;
% $$$     plot(time,filter_field(CIN.WPOW,5*12+1,'-t'),'-k','linewidth',2);
% $$$     xlabel('Year');
% $$$     ylabel('Total Wind Power (W)');
% $$$     xlim([yr1 yr1+nyrs]);
    
% $$$     %%% Time of emergence/historical OHC:
% $$$     
% $$$     figure;
% $$$     plot(time,OHC,'-k');
% $$$     hold on;
% $$$ % $$$     plot([1850 2020],[1 1]*2*PIstd.stdOHC,'--k');
% $$$ % $$$     plot([1850 2020],[1 1]*-2*PIstd.stdOHC,'--k');
% $$$     plot([50 600],[1 1]*2*PIstd.stdOHC,'--k');
% $$$     plot([50 600],[1 1]*-2*PIstd.stdOHC,'--k');
% $$$     xlabel('Year');
% $$$     ylabel('OHC (J)');
% $$$ 
% $$$     figure;
% $$$     plot(PIstd.stdZvP.Tp*2,P,'-k');
% $$$     hold on;
% $$$     plot(PIstd.stdTvP.Tp*2,P,'-r');
% $$$     plot(PIstd.stdYvP.Tp*2,P,'-b');
% $$$     xlabel('2$\sigma$ variability amplitude ($^\circ$C)');
% $$$     ylabel('Percentile');
% $$$     legend('$\Theta(p_z)$','$\Theta(p_\Theta)$',['$\Theta(p_\' ...
% $$$                         'phi)$']);
% $$$     
% $$$     % Example time series:
% $$$     pii = 20;
% $$$     [tmp piii] = min(abs(P-pii));
% $$$     figure;
% $$$     plot(time,ZvP.Tp(piii,:),'-k','linewidth',2);
% $$$     hold on;
% $$$     plot(time,TvP.Tp(piii,:),'-r','linewidth',2);
% $$$ % $$$     plot(time,YvP.Tp(piii,:),'-b','linewidth',2);
% $$$     plot([time(1) time(end)],[1 1]*2*PIstd.stdZvP.Tp(piii),'--k');
% $$$     plot([time(1) time(end)],[1 1]*-2*PIstd.stdZvP.Tp(piii),'--k');
% $$$     plot([time(1) time(end)],[1 1]*2*PIstd.stdTvP.Tp(piii),'--r');
% $$$     plot([time(1) time(end)],[1 1]*-2*PIstd.stdTvP.Tp(piii),'--r');
% $$$ % $$$     plot([time(1) time(end)],[1 1]*2*PIstd.stdYvP.Tp(piii),'--b');
% $$$ % $$$     plot([time(1) time(end)],[1 1]*-2*PIstd.stdYvP.Tp(piii),'--b');

    % Choose whether to remap to T for y-axis:
    remap_to_T = 0;
    if (remap_to_T)
        
        Z_YvPar = ZvP.Tp_mean;
        zlim = [-2 34];
        zlab = 'Horizontally-averaged Temperature ($^\circ$C)';

        T_YvPar = TvP.Tp_mean;
        tlim = [-2 34];
        tlab = 'Percentile-averaged Temperature ($^\circ$C)';
    
        y_YvPar = YvP.Tp_mean;
        yylim = [-80 80];
        ylab = 'Vertical/zonal-averaged Temperature ($^\circ$C)';
    else
        Z_YvPar = P;
        zlim = [0 100];
        zlab = 'Depth Percentile $p_z$';

        T_YvPar = P;
        tlim = zlim;
        tlab = 'Temperature Percentile $p_\Theta$';
        
        y_YvPar = P;
        yylim = zlim;
        ylab = 'Latitude Percentile $p_\phi$';
    end        
    
    %%%%% Plot mean, climatology, std and trends:
    figure;
    set(gcf,'Position',[1 41 2560 1330]);

    subplot(3,4,1);
    plot(ZvP.Tp_mean,P,'linewidth',2);
    xlabel('Temperature ($^\circ$C)');
    xlim([-2 34]);
    ylabel('Depth Percentile');
    ylim([0 100]);
    title('ACCESS-CM2 PI-control mean $\Theta(z)$');
    set(gca,'ydir','reverse')

    subplot(3,4,2);
    [X,Y] = ndgrid(1:12,Z_YvPar);
    contourf(X,Y,ZvP.Tp_clim',[-10 -0.2:0.002:0.2 10],'linestyle','none');
    xlabel('Month');
    ylabel(zlab);
    ylim(zlim);
    title('ACCESS-CM2 PI-control $\Theta(p_z)$ seasonal cycle');
    caxis([-0.2 0.2]);
    colorbar;
    if (~remap_to_T);set(gca,'ydir','reverse');end;

    subplot(3,4,3);
    plot(std(ZvP.Tp_clim,[],2),Z_YvPar,'linewidth',2);
    hold on;
    plot(std(ZvP.Tp,[],2),Z_YvPar,'-r','linewidth',2);
    xlabel('std(Temperature) ($^\circ$C)');
    legend('Seasonal Climatology','Anomalies');
    ylabel(zlab);
    ylim(zlim);
    xlim([0 0.3]);
    title('ACCESS-CM2 PI-control variability $\Theta(p_z)$');
    if (~remap_to_T);set(gca,'ydir','reverse');end;

    subplot(3,4,4);
    plot(ZvP.Tp_cubtr(2,:)',P,'linewidth',2);
    xlabel('Temperature trend ($^\circ$C/year)');
    ylabel('Depth Percentile');
    ylim(zlim);
    xlim([-1.5e-3 1.5e-3]);
    title('ACCESS-CM2 PI-control linear trend $\Theta(p_z)$');
    if (~remap_to_T);set(gca,'ydir','reverse');end;

    subplot(3,4,5);
    plot(TvP.Tp_mean,P,'linewidth',2);
    xlabel('Temperature ($^\circ$C)');
    xlim([-2 34]);
    ylim([0 100]);
    ylabel('Temperature Percentile');
    title('ACCESS-CM2 PI-control mean $\Theta(p)$');
    set(gca,'ydir','reverse')

    subplot(3,4,6);
    [X,Y] = ndgrid(1:12,T_YvPar);
    contourf(X,Y,TvP.Tp_clim',[-10 -0.2:0.002:0.2 10],'linestyle','none');
    xlabel('Month');
    ylabel(tlab);
    title('ACCESS-CM2 PI-control $\Theta(p_\Theta)$ seasonal cycle');
    colorbar;
    colormap(redblue);
    caxis([-0.2 0.2]);
    ylim(tlim);
    if (~remap_to_T);set(gca,'ydir','reverse');end;

    subplot(3,4,7);
    plot(std(TvP.Tp_clim,[],2),T_YvPar,'linewidth',2);
    hold on;
    plot(std(TvP.Tp,[],2),T_YvPar,'-r','linewidth',2);
    xlabel('std(Temperature) ($^\circ$C)');
    legend('Seasonal Climatology','Anomalies');
    ylim(yylim);
    xlim([0 0.3]);
    ylabel(ylab);
    title('ACCESS-CM2 PI-control variability $\Theta(p_\Theta)$');
    if (~remap_to_T);set(gca,'ydir','reverse');end;
% $$$     set(gca,'Position',[0.1300    0.5838    0.2504    0.3412]);

    subplot(3,4,8);
    plot(TvP.Tp_cubtr(2,:)',T_YvPar,'linewidth',2);
    xlabel('Temperature trend ($^\circ$C/year)');
    xlim([-1.5e-3 1.5e-3]);
    ylim(yylim);
    ylabel(ylab);
    title('ACCESS-CM2 PI-control linear trend $\Theta(p_\Theta)$');
    if (~remap_to_T);set(gca,'ydir','reverse');end;

    subplot(3,4,9);
    plot(YvP.Tp_mean,P,'linewidth',2);
    xlabel('Temperature ($^\circ$C)');
    xlim([-1 6.5]);
    ylim([0 100]);
    ylabel('Latitude Percentile');
    title('ACCESS-CM2 PI-control mean $\phi(p)$');

    subplot(3,4,10);
    [X,Y] = ndgrid(1:12,y_YvPar);
    contourf(X,Y,YvP.Tp_clim',[-10 -0.2:0.002:0.2 10],'linestyle','none');
    xlabel('Month');
    ylabel(ylab);
    title('ACCESS-CM2 PI-control $\Theta(p_\phi)$ seasonal cycle');
    colorbar;
    colormap(redblue);
    caxis([-0.2 0.2]);
    ylim(yylim);

    subplot(3,4,11);
    plot(std(YvP.Tp_clim,[],2),y_YvPar,'linewidth',2);
    hold on;
    plot(std(YvP.Tp,[],2),y_YvPar,'-r','linewidth',2);
    xlabel('std(Temperature) ($^\circ$C)');
    legend('Seasonal Climatology','Anomalies');
    ylim(yylim);
    xlim([0 0.3]);
    ylabel(ylab);
    title('ACCESS-CM2 PI-control variability $\Theta(p_\phi)$');
% $$$     set(gca,'Position',[0.1300    0.5838    0.2504    0.3412]);

    subplot(3,4,12);
    plot(YvP.Tp_cubtr(2,:)',y_YvPar,'linewidth',2);
    xlabel('Temperature trend ($^\circ$C/year)');
    xlim([-1.5e-3 1.5e-3]);
    ylim(yylim);
    ylabel(ylab);
    title('ACCESS-CM2 PI-control linear trend $\Theta(p_\phi)$');
% $$$ 
% $$$ 
% $$$     %%% Straight comparison of amplitude of anomalies:
% $$$     figure;
% $$$     set(gcf,'Position',[1 41 2560 1330]);
% $$$     subplot(1,2,1);
% $$$     plot(std(ZvP.Tp,[],2),P,'-k','linewidth',2);
% $$$     hold on;
% $$$     plot(std(TvP.Tp,[],2),P,'-r','linewidth',2);
% $$$     plot(std(YvP.Tp,[],2),P,'-b','linewidth',2);
% $$$     xlabel('Standard Deviation of $\Theta$ ($^\circ$C)');
% $$$     legend('$\Theta(p_z)$','$\Theta(p_\Theta)$','$\Theta(p_\phi)$');
% $$$     ylabel('Percentile');
% $$$     ylim([0 100]);
% $$$     xlim([0 0.3]);
% $$$     title('ACCESS-CM2 PI-control variability');
% $$$     set(gca,'ydir','reverse');
% $$$     
% $$$     % Averages:
% $$$     text(0.1,30,['$\Theta(p_z)$ average = ' sprintf('%3.4f',mean(std(ZvP.Tp,[],2))) '$^\circ$C'],'color','k');
% $$$     text(0.1,35,['$\Theta(p_\Theta)$ average = ' sprintf('%3.4f',mean(std(TvP.Tp,[],2))) '$^\circ$C'],'color','r');
% $$$     text(0.1,40,['$\Theta(p_\phi)$ average = ' sprintf('%3.4f',mean(std(YvP.Tp,[],2))) '$^\circ$C'],'color','b');
% $$$ 
% $$$     ax1=gca;
% $$$     ax1_pos = ax1.Position; % position of first axes
% $$$     ax1_pos(1) = ax1_pos(1)*0.8;
% $$$     ax1_pos(3) = 0.00001;
% $$$     ax2 = axes('Position',ax1_pos,...
% $$$                'Color','none');
% $$$     set(ax2,'FontSize',15);
% $$$     Zticks = [30 18 12:-2:-2];
% $$$     Pticks = zeros(size(Zticks));
% $$$     for ii=1:length(Zticks)
% $$$         Pticks(ii) = interp1(TvP.Tp_mean,P,Zticks(ii),'linear');
% $$$     end
% $$$     set(ax2,'ytick',Pticks);
% $$$     set(ax2,'yticklabel',Zticks);
% $$$     ylabel(ax2,'Temperature ($^\circ$C)');
% $$$     ylim(ax2,zlim);
% $$$     set(ax2,'ydir','reverse');
% $$$     ax1=gca;
% $$$     ax1_pos = ax1.Position; % position of first axes
% $$$     ax1_pos(1) = ax1_pos(1)*0.6;
% $$$     ax1_pos(3) = 0.00001;
% $$$     ax2 = axes('Position',ax1_pos,...
% $$$                    'Color','r');
% $$$     set(ax2,'color','r');
% $$$     set(ax2,'FontSize',15);
% $$$     Zticks = [0:500:4000 5000];
% $$$     Pticks = zeros(size(Zticks));
% $$$     for ii=1:length(Zticks)
% $$$         Pticks(ii) = interp1(zofP_mean,P,-Zticks(ii),'linear');
% $$$     end
% $$$     Pticks(1) = 0;
% $$$     set(ax2,'ytick',Pticks);
% $$$     set(ax2,'yticklabel',Zticks);
% $$$     ylabel(ax2,'Depth (m)');
% $$$     ylim(ax2,zlim);
% $$$     set(ax2,'ydir','reverse');
% $$$     ax1=gca;
% $$$     ax1_pos = ax1.Position; % position of first axes
% $$$     ax1_pos(1) = ax1_pos(1)*0.4;
% $$$     ax1_pos(3) = 0.00001;
% $$$     ax2 = axes('Position',ax1_pos,...
% $$$                    'Color','r');
% $$$     set(ax2,'color','r');
% $$$     set(ax2,'FontSize',15);
% $$$     Zticks = [-75:15:75];
% $$$     Pticks = zeros(size(Zticks));
% $$$     for ii=1:length(Zticks)
% $$$         Pticks(ii) = interp1(yofP_mean,P,Zticks(ii),'linear');
% $$$     end
% $$$     Pticks(1) = 0;
% $$$     set(ax2,'ytick',Pticks);
% $$$     set(ax2,'yticklabel',Zticks);
% $$$     ylabel(ax2,'Latitude ($^\circ$N)');
% $$$     ylim(ax2,zlim);
% $$$     set(gca,'ydir','reverse');
% $$$     
% $$$ 
% $$$     subplot(1,2,2);
% $$$     plot(std(ZvP.Hp,[],2),P,'-k','linewidth',2);
% $$$     hold on;
% $$$     plot(std(TvP.Hp,[],2),P,'-r','linewidth',2);
% $$$     plot(std(YvP.Hp,[],2),P,'-b','linewidth',2);
% $$$     xlabel('Standard Deviation of heat content (J)');
% $$$     legend('$H(p_z)$','$H(p_\Theta)$','$H(p_\phi)$');
% $$$     ylabel('Percentile');
% $$$     ylim([0 100]);
% $$$     title('ACCESS-CM2 PI-control variability');
% $$$     set(gca,'ydir','reverse');
        
    % Choose whether to plot T or H:
    plot_H = 0;
    if (plot_H)
        ZvPar = ZvP.Hp;
        Zcxs = [-0.5e23 1e21 0.5e23];
        Zlab = '$H(p_z)$';

        TvPar = TvP.Hp;
        Tcxs = [-0.5e23 1e21 0.5e23];
        Tlab = '$H(p_\Theta)$';
    
        yvar = YvP.Hp;
        ycxs = [-0.5e23 1e21 0.5e23];
        Ylab = '$H(p_\phi)$';
    else
        ZvPar = ZvP.Tp;
        Zcxs = [-0.2 0.02 0.2];
        Zlab = '$\Theta(p_z)$';
        TvPar = TvP.Tp;
        Tcxs = Zcxs;%[-0.1 0.01 0.1];
        Tlab = '$\Theta(p_\Theta)$';
        yvar = YvP.Tp;
        ycxs = Zcxs;%[-0.1 0.01 0.1];
        Ylab = '$\Theta(p_\phi)$';
    end        
        
    %%%%%% Plot all anomalies time series:

    xlims = [yr1 yr1+nyrs];

    % ENSO Focus settings:
    zlim = [0 10];
    tlim = [0 10];
    yylim = [40 85];
    xlims = [300 350];
    Zcxs = [-0.2 0.01 0.2];
    Tcxs = [-0.2 0.01 0.2];
    ycxs = [-0.2 0.01 0.2];

    [tmp t1] = min(abs(time-xlims(1)));
    [tmp t2] = min(abs(time-xlims(2)));
    
    figure;
    set(gcf,'Position',[1 41 2560 1330]);
    set(gcf,'defaulttextfontsize',15);
    set(gcf,'defaultaxesfontsize',15);
    [X,Y] = ndgrid(time(t1:t2),Z_YvPar);
    subplot(3,1,1);
    contourf(X,Y,ZvPar(:,t1:t2)',[-1e50 Zcxs(1):Zcxs(2):Zcxs(3) 1e50],'linestyle','none');
    caxis([Zcxs(1) Zcxs(3)]);
    ylim(zlim);
    xlim(xlims);
% $$$     xlabel('Year');
    set(gca,'xticklabel',[]);
    ylabel(zlab);
% $$$     cb = colorbar;
% $$$     ylabel(cb,'Temperature Anomalies ($^\circ$C)','Interpreter','latex');
    text(xlims(1)+(xlims(2)-xlims(1))*0.01,0.93*zlim(2),['ACCESS-CM2 PI-control ' Zlab ' anomalies'],'Backgroundcolor','w');
    set(gca,'Position',[0.13 0.69 0.75 0.29]);
    if (~remap_to_T)
        set(gca,'ydir','reverse');
        ax1=gca;
        ax1_pos = ax1.Position; % position of first axes
        ax1_pos(1) = ax1_pos(1)*0.7;
        ax1_pos(3) = 0.00001;
        ax2 = axes('Position',ax1_pos,...
                   'Color','none');
        set(ax2,'FontSize',15);
        Zticks = [0:500:4000 5000];
% $$$         Zticks = [0:50:800];
% $$$         Zticks = [3000:250:5500];
        Pticks = zeros(size(Zticks));
        for ii=1:length(Zticks)
            Pticks(ii) = interp1(zofP_mean,P,-Zticks(ii),'linear');
        end
        Pticks(1) = 0;
        set(ax2,'ytick',Pticks);
        set(ax2,'yticklabel',Zticks);
        ylabel(ax2,'Depth (m)');
        ylim(ax2,zlim);
        set(ax2,'ydir','reverse');
    end        
    
    [X,Y] = ndgrid(time(t1:t2),T_YvPar);
    subplot(3,1,2);
    contourf(X,Y,TvPar(:,t1:t2)',[-1e50 Tcxs(1):Tcxs(2):Tcxs(3) 1e50],'linestyle','none');
    caxis([Tcxs(1) Tcxs(3)]);
    ylim(tlim);
    xlim(xlims);
    colormap(redblue);
% $$$     xlabel('Year');
    set(gca,'xticklabel',[]);
    ylabel(tlab);
    cb = colorbar;
    ylabel(cb,'Temperature Anomalies ($^\circ$C)','Interpreter','latex');
    text(xlims(1)+(xlims(2)-xlims(1))*0.01,0.73*tlim(2),['ACCESS-CM2 PI-control ' Tlab ' anomalies'],'Backgroundcolor','w');
    set(gca,'Position',[0.13 0.37 0.75 0.29]);
    if (~remap_to_T)
        set(gca,'ydir','reverse');
        ax1=gca;
        ax1_pos = ax1.Position; % position of first axes
        ax1_pos(1) = ax1_pos(1)*0.7;
        ax1_pos(3) = 0.00001;
        ax2 = axes('Position',ax1_pos,...
                   'Color','none');
        set(ax2,'FontSize',15);
% $$$         Zticks = [30 18 12:-2:-2];
        Zticks = [30 15 10:-2:0];
% $$$         Zticks = [30 22 18 16:-2:8];
% $$$         Zticks = [2:-0.25:-0.5 -2];
        Pticks = zeros(size(Zticks));
        for ii=1:length(Zticks)
            Pticks(ii) = interp1(TvP.Tp_mean,P,Zticks(ii),'linear');
        end
        set(ax2,'ytick',Pticks);
        set(ax2,'yticklabel',Zticks);
        ylabel(ax2,'Temperature ($^\circ$C)');
        ylim(ax2,zlim);
        set(ax2,'ydir','reverse');
    end

    [X,Y] = ndgrid(time(t1:t2),y_YvPar);
    subplot(3,1,3);
    contourf(X,Y,yvar(:,t1:t2)',[-1e50 ycxs(1):ycxs(2):ycxs(3) 1e50],'linestyle','none');
    caxis([ycxs(1) ycxs(3)]);
    ylim(yylim);
    xlim(xlims);
    xlabel('Year');
    ylabel(ylab);
% $$$     cb = colorbar;
% $$$     ylabel(cb,'Temperature Anomalies ($^\circ$C)','Interpreter','latex');
    text(xlims(1)+(xlims(2)-xlims(1))*0.01,83,['ACCESS-CM2 PI-control ' Ylab ' anomalies'],'Backgroundcolor','w');
    set(gca,'Position',[0.13 0.05 0.75 0.29]);
    if (~remap_to_T)
        ax1=gca;
        ax1_pos = ax1.Position; % position of first axes
        ax1_pos(1) = ax1_pos(1)*0.7;
        ax1_pos(3) = 0.00001;
        ax2 = axes('Position',ax1_pos,...
                   'Color','none');
        set(ax2,'FontSize',15);
        Zticks = [-75:15:60];
        Pticks = zeros(size(Zticks));
        for ii=1:length(Zticks)
            Pticks(ii) = interp1(yofP_mean,P,Zticks(ii),'linear');
        end
        Pticks(1) = 0;
        set(ax2,'ytick',Pticks);
        set(ax2,'yticklabel',Zticks);
        ylabel(ax2,'Latitude ($^\circ$N)');
        ylim(ax2,yylim);
    end

    %%%%%% Plot all anomalies time series:
    xlims = [yr1 yr1+nyrs];

    % ENSO Focus settings:
    zlim = [0 10];
    tlim = [0 10];
    yylim = [40 85];
    xlims = [300 350];
    Zcxs = [-0.2 0.01 0.2];
    Tcxs = [-0.2 0.01 0.2];
    ycxs = [-0.2 0.01 0.2];

    [tmp t1] = min(abs(time-xlims(1)));
    [tmp t2] = min(abs(time-xlims(2)));
    
    bvars = {'TENp','ADVp','FORp','RMIXp','VMIXp'};

    for ti =1:length(typs)
        figure;
        set(gcf,'Position',[1 41 2560 1330]);
        set(gcf,'defaulttextfontsize',15);
        set(gcf,'defaultaxesfontsize',15);
        [X,Y] = ndgrid(time(t1:t2),Z_YvPar);

        subplot(3,2,1);
        contourf(X,Y,ZvPar(:,t1:t2)',[-1e50 Zcxs(1):Zcxs(2):Zcxs(3) 1e50],'linestyle','none');
        caxis([Zcxs(1) Zcxs(3)]);
        ylim(zlim);
        xlim(xlims);
        set(gca,'xticklabel',[]);
        ylabel(zlab);
        text(xlims(1)+(xlims(2)-xlims(1))*0.01,0.93*zlim(2),['ACCESS-CM2 PI-control ' Zlab ' anomalies'],'Backgroundcolor','w');
        set(gca,'ydir','reverse');
        
        for vi = 1:length(bvars)
    
    [X,Y] = ndgrid(time(t1:t2),T_YvPar);
    subplot(3,1,2);
    contourf(X,Y,TvPar(:,t1:t2)',[-1e50 Tcxs(1):Tcxs(2):Tcxs(3) 1e50],'linestyle','none');
    caxis([Tcxs(1) Tcxs(3)]);
    ylim(tlim);
    xlim(xlims);
    colormap(redblue);
% $$$     xlabel('Year');
    set(gca,'xticklabel',[]);
    ylabel(tlab);
    cb = colorbar;
    ylabel(cb,'Temperature Anomalies ($^\circ$C)','Interpreter','latex');
    text(xlims(1)+(xlims(2)-xlims(1))*0.01,0.73*tlim(2),['ACCESS-CM2 PI-control ' Tlab ' anomalies'],'Backgroundcolor','w');
    set(gca,'Position',[0.13 0.37 0.75 0.29]);
    if (~remap_to_T)
        set(gca,'ydir','reverse');
        ax1=gca;
        ax1_pos = ax1.Position; % position of first axes
        ax1_pos(1) = ax1_pos(1)*0.7;
        ax1_pos(3) = 0.00001;
        ax2 = axes('Position',ax1_pos,...
                   'Color','none');
        set(ax2,'FontSize',15);
% $$$         Zticks = [30 18 12:-2:-2];
        Zticks = [30 15 10:-2:0];
% $$$         Zticks = [30 22 18 16:-2:8];
% $$$         Zticks = [2:-0.25:-0.5 -2];
        Pticks = zeros(size(Zticks));
        for ii=1:length(Zticks)
            Pticks(ii) = interp1(TvP.Tp_mean,P,Zticks(ii),'linear');
        end
        set(ax2,'ytick',Pticks);
        set(ax2,'yticklabel',Zticks);
        ylabel(ax2,'Temperature ($^\circ$C)');
        ylim(ax2,zlim);
        set(ax2,'ydir','reverse');
    end

    [X,Y] = ndgrid(time(t1:t2),y_YvPar);
    subplot(3,1,3);
    contourf(X,Y,yvar(:,t1:t2)',[-1e50 ycxs(1):ycxs(2):ycxs(3) 1e50],'linestyle','none');
    caxis([ycxs(1) ycxs(3)]);
    ylim(yylim);
    xlim(xlims);
    xlabel('Year');
    ylabel(ylab);
% $$$     cb = colorbar;
% $$$     ylabel(cb,'Temperature Anomalies ($^\circ$C)','Interpreter','latex');
    text(xlims(1)+(xlims(2)-xlims(1))*0.01,83,['ACCESS-CM2 PI-control ' Ylab ' anomalies'],'Backgroundcolor','w');
    set(gca,'Position',[0.13 0.05 0.75 0.29]);
    if (~remap_to_T)
        ax1=gca;
        ax1_pos = ax1.Position; % position of first axes
        ax1_pos(1) = ax1_pos(1)*0.7;
        ax1_pos(3) = 0.00001;
        ax2 = axes('Position',ax1_pos,...
                   'Color','none');
        set(ax2,'FontSize',15);
        Zticks = [-75:15:60];
        Pticks = zeros(size(Zticks));
        for ii=1:length(Zticks)
            Pticks(ii) = interp1(yofP_mean,P,Zticks(ii),'linear');
        end
        Pticks(1) = 0;
        set(ax2,'ytick',Pticks);
        set(ax2,'yticklabel',Zticks);
        ylabel(ax2,'Latitude ($^\circ$N)');
        ylim(ax2,yylim);
    end


% $$$ 
% $$$     % Time series to go with that:
% $$$     figure;
% $$$     set(gcf,'Position',[1 41 2560 1330]);
% $$$     set(gcf,'defaulttextfontsize',15);
% $$$     set(gcf,'defaultaxesfontsize',15);
% $$$     subplot(2,1,1);
% $$$     plot(time(t1:t2),N34(t1:t2),'-k');
% $$$     hold on;
% $$$     plot(time(t1:t2),GMSST(t1:t2)'*10,'-r');
% $$$     plot(time(t1:t2),OHC(t1:t2)/1e22,'-','color',[0.5 0.5 0.5]);
% $$$     legend('Nino 3.4','Global SST * 10','Global OHC anomaly ($/10^{22} J$)');
% $$$     ylabel('SST ($^\circ$C)');
% $$$     xlim(xlims);
% $$$     xlabel('Year');
% $$$ 
% $$$     %%%%% Spectral analysis:
% $$$     fs = 12;
% $$$     [pxx,f] = pmtm(ZvPar(15,:)',[],[],fs);
% $$$     fL = length(f);
% $$$     TvPar_spec = zeros(fL,PL);
% $$$     ZvPar_spec = zeros(fL,PL);
% $$$     yvar_spec = zeros(fL,PL);
% $$$     for pi=1:PL
% $$$         [ZvPar_spec(:,pi),~] = pmtm(ZvPar(pi,:)',[],[],fs);
% $$$         [TvPar_spec(:,pi),~] = pmtm(TvPar(pi,:)',[],[],fs);
% $$$         [yvar_spec(:,pi),~] = pmtm(yvar(pi,:)',[],[],fs);
% $$$     end
% $$$     
% $$$     figure;
% $$$     set(gcf,'Position',[1 41 2560 1330]);
% $$$     subplot(3,1,1);
% $$$     [X,Y] = ndgrid(f,Z_YvPar);
% $$$     contourf(X,Y,log10(ZvPar_spec),[-4.5:0.25:-0.25],'linestyle','none');
% $$$     xlim([1/250 1]);    
% $$$     xlab = [300 100 50 10 7 5 3 2 1];
% $$$     ylim(zlim);
% $$$     set(gca,'xscale','log');
% $$$     set(gca,'xtick',1./xlab);
% $$$     set(gca,'xticklabel',xlab);
% $$$     xlabel('Period (years)');
% $$$     caxis([-4.5 -1]);
% $$$     ylabel(zlab);
% $$$     title(['Multi-taper power spectra of ACCESS-CM2 PI-control ' ...
% $$$            '$\Theta(p_z)$ anomalies ($^\circ$C$^2$/year)']);
% $$$     cb = colorbar;
% $$$     set(cb,'ytick',[-4:1:-1]);
% $$$     set(cb,'yticklabel',10.^[-4:1:-1]);
% $$$     if (~remap_to_T)
% $$$         set(gca,'ydir','reverse');
% $$$         ax1=gca;
% $$$         ax1_pos = ax1.Position; % position of first axes
% $$$         ax1_pos(1) = ax1_pos(1)*0.7;
% $$$         ax1_pos(3) = 0.00001;
% $$$         ax2 = axes('Position',ax1_pos,...
% $$$                    'Color','none');
% $$$         set(ax2,'FontSize',15);
% $$$         Zticks = [0:500:4000 5000];
% $$$ % $$$         Zticks = [0:200:800];
% $$$ % $$$         Zticks = [3000:250:5500];
% $$$         Pticks = zeros(size(Zticks));
% $$$         for ii=1:length(Zticks)
% $$$             Pticks(ii) = interp1(zofP_mean,P,-Zticks(ii),'linear');
% $$$         end
% $$$         Pticks(1) = 0;
% $$$         set(ax2,'ytick',Pticks);
% $$$         set(ax2,'yticklabel',Zticks);
% $$$         ylabel(ax2,'Depth (m)');
% $$$         ylim(ax2,zlim);
% $$$         set(ax2,'ydir','reverse');
% $$$     end        
% $$$     
% $$$     subplot(3,1,2);
% $$$     [X,Y] = ndgrid(f,T_YvPar);
% $$$     contourf(X,Y,log10(TvPar_spec),[-4.5:0.25:-0.25],'linestyle','none');
% $$$     xlim([1/250 1]);    
% $$$     xlab = [300 100 50 10 7 5 3 2 1];
% $$$     ylim(tlim);
% $$$     set(gca,'xscale','log');
% $$$     set(gca,'xtick',1./xlab);
% $$$     set(gca,'xticklabel',xlab);
% $$$     xlabel('Period (years)');
% $$$     caxis([-4.5 -1]);
% $$$     ylabel(tlab);
% $$$     title(['Multi-taper power spectra of ACCESS-CM2 PI-control ' ...
% $$$            '$\Theta(p_\Theta)$ anomalies ($^\circ$C$^2$/year)']);
% $$$     colorbar;
% $$$     cb = colorbar;
% $$$     set(cb,'ytick',[-4:1:-1]);
% $$$     set(cb,'yticklabel',10.^[-4:1:-1]);
% $$$     if (~remap_to_T)
% $$$         set(gca,'ydir','reverse');
% $$$         ax1=gca;
% $$$         ax1_pos = ax1.Position; % position of first axes
% $$$         ax1_pos(1) = ax1_pos(1)*0.7;
% $$$         ax1_pos(3) = 0.00001;
% $$$         ax2 = axes('Position',ax1_pos,...
% $$$                    'Color','none');
% $$$         set(ax2,'FontSize',15);
% $$$         Zticks = [30 18 12:-2:-2];
% $$$ % $$$         Zticks = [30 22 18 16:-2:8];
% $$$ % $$$         Zticks = [2:-0.25:-0.5 -2];
% $$$         Pticks = zeros(size(Zticks));
% $$$         for ii=1:length(Zticks)
% $$$             Pticks(ii) = interp1(TvP.Tp_mean,P,Zticks(ii),'linear');
% $$$         end
% $$$         set(ax2,'ytick',Pticks);
% $$$         set(ax2,'yticklabel',Zticks);
% $$$         ylabel(ax2,'Temperature ($^\circ$C)');
% $$$         ylim(ax2,zlim);
% $$$         set(ax2,'ydir','reverse');
% $$$     end
% $$$ 
% $$$     subplot(3,1,3);
% $$$     [X,Y] = ndgrid(f,y_YvPar);
% $$$     contourf(X,Y,log10(yvar_spec),[-4.5:0.25:-0.25],'linestyle','none');
% $$$     xlim([1/250 1]);    
% $$$     xlab = [300 100 50 10 7 5 3 2 1];
% $$$     ylim(yylim);
% $$$     set(gca,'xscale','log');
% $$$     set(gca,'xtick',1./xlab);
% $$$     set(gca,'xticklabel',xlab);
% $$$     xlabel('Period (years)');
% $$$     caxis([-4.5 -1]);
% $$$     ylabel(ylab);
% $$$     title(['Multi-taper power spectra of ACCESS-CM2 PI-control ' ...
% $$$            '$\Theta(p_\phi)$ anomalies ($^\circ$C$^2$/year)']);
% $$$     colorbar;
% $$$     cb = colorbar;
% $$$     set(cb,'ytick',[-4:1:-1]);
% $$$     set(cb,'yticklabel',10.^[-4:1:-1]);
% $$$     if (~remap_to_T)
% $$$         ax1=gca;
% $$$         ax1_pos = ax1.Position; % position of first axes
% $$$         ax1_pos(1) = ax1_pos(1)*0.7;
% $$$         ax1_pos(3) = 0.00001;
% $$$         ax2 = axes('Position',ax1_pos,...
% $$$                    'Color','none');
% $$$         set(ax2,'FontSize',15);
% $$$         Zticks = [-75:15:75];
% $$$         Pticks = zeros(size(Zticks));
% $$$         for ii=1:length(Zticks)
% $$$             Pticks(ii) = interp1(yofP_mean,P,Zticks(ii),'linear');
% $$$         end
% $$$         Pticks(1) = 0;
% $$$         set(ax2,'ytick',Pticks);
% $$$         set(ax2,'yticklabel',Zticks);
% $$$         ylabel(ax2,'Latitude ($^\circ$N)');
% $$$         ylim(ax2,yylim);
% $$$     end
% $$$ 
% $$$     colormap(cmocean('dense'));

% $$$     % Spectral average:
% $$$     figure;
% $$$     set(gcf,'Position',[421.7 537.7 1213.3 614.7]);
% $$$     plot(f,log10(mean(ZvPar_spec,2)),'-k','linewidth',2);
% $$$     hold on;
% $$$     plot(f,log10(mean(TvPar_spec,2)),'-r','linewidth',2);
% $$$     plot(f,log10(mean(yvar_spec,2)),'-b','linewidth',2);
% $$$     xlim([1/250 1]);    
% $$$     xlab = [300 100 50 10 7 5 3 2 1];
% $$$     set(gca,'xscale','log');
% $$$     set(gca,'xtick',1./xlab);
% $$$     set(gca,'xticklabel',xlab);
% $$$     xlabel('Period (years)');
% $$$     legend('$\Theta(p_z)$','$\Theta(p_\Theta)$','$\Theta(p_\phi)$');
% $$$     ylabel('Multi-taper power spectra ($^\circ$C$^2$/year)');
% $$$     title('ACCESS-CM2 PI-control Percentile-Averaged Spectra');
% $$$ 
% $$$     % Normalized index spectra:
% $$$     fs = 12;
% $$$     fL = length(f);
% $$$     figure;
% $$$     set(gcf,'Position',[421.7 537.7 1213.3 614.7]);
% $$$     [pxx,f] = pmtm(N34/std(N34),5,[],fs);
% $$$     plot(f,log10(pxx),'-k','linewidth',2);
% $$$     hold on;
% $$$     [pxx,f] = pmtm(AMOC/std(AMOC),5,[],fs);
% $$$     plot(f,log10(pxx),'-b','linewidth',2);
% $$$     ts = filter_field(AMOC,5*12+1,'-t');
% $$$     ts(isnan(ts))= 0;
% $$$     [pxx,f] = pmtm(ts/std(ts),5,[],fs);
% $$$     plot(f,log10(pxx),'--b','linewidth',2);
% $$$     [pxx,f] = pmtm(WPOW/std(WPOW),5,[],fs);
% $$$     plot(f,log10(pxx),'-r','linewidth',2);
% $$$     [pxx,f] = pmtm(GMSST/std(GMSST),5,[],fs);
% $$$     plot(f,log10(pxx),'-','linewidth',2,'color',[0 0.5 0]);
% $$$     [pxx,f] = pmtm(TPI/std(TPI),5,[],fs);
% $$$     plot(f,log10(pxx),'-m','linewidth',2);
% $$$     xlim([1/250 1]);    
% $$$     xlab = [300 100 50 10 7 5 3 2 1];
% $$$     set(gca,'xscale','log');
% $$$     set(gca,'xtick',1./xlab);
% $$$     set(gca,'xticklabel',xlab);
% $$$     xlabel('Period (years)');
% $$$     ylabel('Multi-taper power spectra ($^\circ$C$^2$/year)');
% $$$     legend('Nino 3.4','AMOC','AMOC 5-year low-pass','Wind Power','GMSST','IPO (TPI)');
% $$$     title('ACCESS-CM2 PI-control Normalized Index Spectra');
% $$$ 
    
    %%% EOF analysis:

    nmod = 5;
    ntot = 30;
    vars = {'z','T','y'};
    
    for vi=1:length(vars)
        eval(['X = T' vars{vi} 'p'';']);

        C = X'*X;
        tmp = eigs(C,ntot);
        varfrac = tmp./sum(tmp);
        [V,D] = eigs(C,nmod);
        PC = X*V;
        scale = std(PC,[],1);
        PC = PC./repmat(scale,[tL 1]);
        EOFs = ZvP.Tp*PC/tL;
        
        eval(['PC' vars{vi} ' = PC;']);
        eval(['EOF' vars{vi} ' = EOFs;']);
    end

    cols = {'k','r','b','m',[0 0.5 0]};
    figure;
    subplot(1,3,1);
    for i=1:nmod
        plot(EOFz(:,i),P,'-','color',cols{i});
        hold on;
    end
    xlabel('EOFs');
    ylabel('Depth Percentile');
    set(gca,'ydir','reverse');
    legend('EOF1','EOF2','EOF3','EOF4','EOF5');
    subplot(1,3,2);
    for i=1:nmod
        plot(EOFT(:,i),P,'-','color',cols{i});
        hold on;
    end
    set(gca,'ydir','reverse');
    xlabel('EOFs');
    ylabel('Temperature Percentile');
    subplot(1,3,3);
    for i=1:nmod
        plot(EOFy(:,i),P,'-','color',cols{i});
        hold on;
    end
    set(gca,'ydir','reverse');
    xlabel('EOFs');
    ylabel('Latitude Percentile');
    
    figure;
    subplot(3,1,1);
    for i=1:nmod
        plot(time,PCz(:,i),'-','color',cols{i});
        hold on;
    end
    legend('EOF1','EOF2','EOF3','EOF4','EOF5');
    xlabel('Time');
    ylabel('PC');
    xlim([300 450]);
    subplot(3,1,2);
    for i=1:nmod
        plot(time,PCT(:,i),'-','color',cols{i});
        hold on;
    end
    xlabel('Time');
    ylabel('PC');
    xlim([300 450]);
    subplot(3,1,3);
    for i=1:nmod
        plot(time,PCy(:,i),'-','color',cols{i});
        hold on;
    end
    xlabel('Time');
    ylabel('PC');
    xlim([300 450]);

    xlims = [300 380];
    [tmp t1] = min(abs(time-xlims(1)));
    [tmp t2] = min(abs(time-xlims(2)));
    Zcxs = [-0.1 0.0025 0.1];
    zlim = [0 30];
    var = 'z';
    neof = [3 4];

    figure;
    set(gcf,'Position',[1 41 2560 1330]);
    set(gcf,'defaulttextfontsize',15);
    set(gcf,'defaultaxesfontsize',15);
    [X,Y] = ndgrid(time(t1:t2),Z_YvPar);
    subplot(3,1,1);
    contourf(X,Y,ZvPar(:,t1:t2)',[-1e50 Zcxs(1):Zcxs(2):Zcxs(3) 1e50],'linestyle','none');
    caxis([Zcxs(1) Zcxs(3)]);
    ylim(zlim);
    xlim(xlims);
% $$$     xlabel('Year');
    set(gca,'xticklabel',[]);
    ylabel(zlab);
    cb = colorbar;
    ylabel(cb,'Temperature Anomalies ($^\circ$C)','Interpreter','latex');
    text(xlims(1)+(xlims(2)-xlims(1))*0.01,0.93*zlim(2),['ACCESS-CM2 PI-control ' Zlab ' anomalies'],'Backgroundcolor','w');
    set(gca,'Position',[0.13 0.69 0.75 0.29]);
    if (~remap_to_T)
        set(gca,'ydir','reverse');
        ax1=gca;
        ax1_pos = ax1.Position; % position of first axes
        ax1_pos(1) = ax1_pos(1)*0.7;
        ax1_pos(3) = 0.00001;
        ax2 = axes('Position',ax1_pos,...
                   'Color','none');
        set(ax2,'FontSize',15);
        Zticks = [0:500:4000 5000];
% $$$         Zticks = [0:50:800];
% $$$         Zticks = [3000:250:5500];
        Pticks = zeros(size(Zticks));
        for ii=1:length(Zticks)
            Pticks(ii) = interp1(zofP_mean,P,-Zticks(ii),'linear');
        end
        Pticks(1) = 0;
        set(ax2,'ytick',Pticks);
        set(ax2,'yticklabel',Zticks);
        ylabel(ax2,'Depth (m)');
        ylim(ax2,zlim);
        set(ax2,'ydir','reverse');
    end        

    subplot(3,1,2);
    eval(['varF = EOF' var '(:,neof(1))*PC' var '(:,neof(1))'';']);
    contourf(X,Y,varF(:,t1:t2)',[-1e50 Zcxs(1):Zcxs(2):Zcxs(3) 1e50],'linestyle','none');
    caxis([Zcxs(1) Zcxs(3)]);
    ylim(zlim);
    xlim(xlims);
% $$$     xlabel('Year');
    set(gca,'xticklabel',[]);
    ylabel(zlab);
% $$$     cb = colorbar;
% $$$     ylabel(cb,'Temperature Anomalies ($^\circ$C)','Interpreter','latex');
    text(xlims(1)+(xlims(2)-xlims(1))*0.01,0.93*zlim(2),['ACCESS-CM2 ' ...
                        'PI-control EOF ' num2str(neof(1))],'Backgroundcolor','w');
    set(gca,'Position',[0.13 0.37 0.75 0.29]);
    if (~remap_to_T)
        set(gca,'ydir','reverse');
        ax1=gca;
        ax1_pos = ax1.Position; % position of first axes
        ax1_pos(1) = ax1_pos(1)*0.7;
        ax1_pos(3) = 0.00001;
        ax2 = axes('Position',ax1_pos,...
                   'Color','none');
        set(ax2,'FontSize',15);
        Zticks = [0:500:4000 5000];
% $$$         Zticks = [0:50:800];
% $$$         Zticks = [3000:250:5500];
        Pticks = zeros(size(Zticks));
        for ii=1:length(Zticks)
            Pticks(ii) = interp1(zofP_mean,P,-Zticks(ii),'linear');
        end
        Pticks(1) = 0;
        set(ax2,'ytick',Pticks);
        set(ax2,'yticklabel',Zticks);
        ylabel(ax2,'Depth (m)');
        ylim(ax2,zlim);
        set(ax2,'ydir','reverse');
    end        

    subplot(3,1,3);
    eval(['varF = EOF' var '(:,neof(2))*PC' var '(:,neof(2))'';']);
    contourf(X,Y,varF(:,t1:t2)',[-1e50 Zcxs(1):Zcxs(2):Zcxs(3) 1e50],'linestyle','none');
    caxis([Zcxs(1) Zcxs(3)]);
    ylim(zlim);
    xlim(xlims);
    xlabel('Year');
% $$$     set(gca,'xticklabel',[]);
    ylabel(zlab);
% $$$     cb = colorbar;
% $$$     ylabel(cb,'Temperature Anomalies ($^\circ$C)','Interpreter','latex');
    text(xlims(1)+(xlims(2)-xlims(1))*0.01,0.93*zlim(2),['ACCESS-CM2 ' ...
                        'PI-control EOF ' num2str(neof(2))],'Backgroundcolor','w');
    set(gca,'Position',[0.13 0.05 0.75 0.29]);
    if (~remap_to_T)
        set(gca,'ydir','reverse');
        ax1=gca;
        ax1_pos = ax1.Position; % position of first axes
        ax1_pos(1) = ax1_pos(1)*0.7;
        ax1_pos(3) = 0.00001;
        ax2 = axes('Position',ax1_pos,...
                   'Color','none');
        set(ax2,'FontSize',15);
        Zticks = [0:500:4000 5000];
% $$$         Zticks = [0:50:800];
% $$$         Zticks = [3000:250:5500];
        Pticks = zeros(size(Zticks));
        for ii=1:length(Zticks)
            Pticks(ii) = interp1(zofP_mean,P,-Zticks(ii),'linear');
        end
        Pticks(1) = 0;
        set(ax2,'ytick',Pticks);
        set(ax2,'yticklabel',Zticks);
        ylabel(ax2,'Depth (m)');
        ylim(ax2,zlim);
        set(ax2,'ydir','reverse');
    end        
% $$$ 
% $$$     var = 'z';
% $$$     nmod = 3;
% $$$     figure;
% $$$     for i=1:nmod
% $$$         subplot(nmod,1,i);
% $$$         [X,Y] = ndgrid(time,P);
% $$$         pcolPlot(X,Y,varF');
% $$$         set(gca,'ydir','reverse');
% $$$         xlabel('Time');
% $$$         ylabel('Percentile');
% $$$         caxis([-0.2 0.2]);
% $$$     end
    colormap('redblue');
    
    
    
    
        

    
    
    
% $$$     %%%% Plot Regressions:
% $$$     
% $$$     % ENSO:
% $$$     % $$$     ts = GMSST/std(GMSST);
% $$$ % $$$     label = 'GMSST';
% $$$     ts = N34/std(N34);
% $$$     label = 'ENSO';
% $$$     lags = [-12*2:1:12*2]; %months lag
% $$$     zlim = [0 10];
% $$$     tlim = [0 7];
% $$$     yylim = [40 85];
% $$$     Zcxs = [-0.1 0.001 0.1]; Tcxs = Zcxs;Ycxs=Tcxs;
% $$$     
% $$$     % OHC:
% $$$     ts = OHC/std(OHC);
% $$$     label = 'OHC';
% $$$     lags = [-12*100:1:12*100]; %months lag
% $$$     zlim = [0 50];
% $$$     tlim = [0 50];
% $$$     yylim = [0 100];
% $$$     Zcxs = [-0.03 0.001 0.03]; Tcxs = Zcxs;Ycxs=Tcxs;
% $$$ 
% $$$     % AMOC:
% $$$     ts = filter_field(AMOC,5*12+1,'-t');
% $$$     ts(isnan(ts)) = 0;
% $$$     ts = ts/std(ts);
% $$$     label = 'AMOC';
% $$$     lags = [-12*30:1:12*30]; %months lag
% $$$     zlim = [0 40];
% $$$     tlim = [0 40];
% $$$     yylim = [0 100];
% $$$     Zcxs = [-0.02 0.001 0.02]; Tcxs = Zcxs;Ycxs=Tcxs;
% $$$ 
% $$$     % IPO:
% $$$     ts = TPI;
% $$$     ts = ts/std(ts);
% $$$     label = 'IPO';
% $$$     lags = [-12*30:1:12*30]; %months lag
% $$$     zlim = [0 40];
% $$$     tlim = [0 40];
% $$$     yylim = [0 100];
% $$$     Zcxs = [-0.02 0.001 0.02]; Tcxs = Zcxs;Ycxs=Tcxs;
% $$$ 
% $$$     % WPOW:
% $$$     ts = WPOW;
% $$$     ts = ts/std(ts);
% $$$     label = 'Wind Power';
% $$$     lags = [-12*2:1:12*2]; %months lag
% $$$     zlim = [0 10];
% $$$     tlim = [0 7];
% $$$     yylim = [40 85];
% $$$     Zcxs = [-0.05 0.001 0.05]; Tcxs = Zcxs;Ycxs=Tcxs;
% $$$ 
% $$$     % WPOW low-pass:
% $$$     ts = filter_field(WPOW,10*12+1,'-t');
% $$$     ts(isnan(ts))= 0;
% $$$     ts = ts/std(ts);
% $$$     label = 'Wind Power 10-year low-pass';
% $$$     lags = [-12*30:1:12*30]; %months lag
% $$$     zlim = [0 40];
% $$$     tlim = [0 40];
% $$$     yylim = [0 100];
% $$$     Zcxs = [-0.02 0.001 0.02]; Tcxs = Zcxs;Ycxs=Tcxs;
% $$$ 
% $$$     
% $$$     ll = length(lags);
% $$$ 
% $$$     ZvPar_lr = zeros(PL,ll);
% $$$     TvPar_lr = zeros(PL,ll);
% $$$     yvar_lr = zeros(PL,ll);
% $$$     for ii=1:ll
% $$$         lag = lags(ii);
% $$$         TS = zeros(size(ts));
% $$$         if (lag<0)
% $$$             TS(1:(end+lag)) = ts((-lag+1):end);
% $$$             TS((end+lag+1):end) = 0;
% $$$         elseif (lag == 0)
% $$$             TS = ts;
% $$$         else
% $$$             TS(1:(lag)) = 0;
% $$$             TS((lag+1):end) = ts(1:(end-lag));
% $$$         end
% $$$         
% $$$         ZvPar_lr(:,ii) = ZvPar*TS/(sum(TS.^2));
% $$$         TvPar_lr(:,ii) = TvPar*TS/(sum(TS.^2));
% $$$         yvar_lr(:,ii) = yvar*TS/(sum(TS.^2));
% $$$     end
% $$$ 
% $$$     figure;
% $$$     set(gcf,'Position',[1 41 2560 1330]);
% $$$ 
% $$$     subplot(3,1,1);
% $$$     [X,Y] = ndgrid(lags/12,Z_YvPar);
% $$$     contourf(X,Y,ZvPar_lr',[-1e50 Zcxs(1):Zcxs(2):Zcxs(3) 1e50],'linestyle','none');
% $$$     caxis([Zcxs(1) Zcxs(3)]);
% $$$     xlabel('Lag (years)');
% $$$     ylabel(zlab);
% $$$     ylim(zlim);
% $$$     title(['ACCESS-CM2 PI-control $\Theta(p_z)$' label ' regression']);
% $$$     cb = colorbar;
% $$$     ylabel(cb,'Temperature Anomalies ($^\circ$C / $\sigma$)','Interpreter','latex');
% $$$     if (~remap_to_T)
% $$$         set(gca,'ydir','reverse');
% $$$         ax1=gca;
% $$$         ax1_pos = ax1.Position; % position of first axes
% $$$         ax1_pos(1) = ax1_pos(1)*0.7;
% $$$         ax1_pos(3) = 0.00001;
% $$$         ax2 = axes('Position',ax1_pos,...
% $$$                    'Color','none');
% $$$         set(ax2,'FontSize',15);
% $$$         Zticks = [0:500:4000 5000];
% $$$ % $$$         Zticks = [0:200:2000];
% $$$ % $$$         Zticks = [0:100:800];
% $$$         Zticks = [0:250:4000 5000];
% $$$ % $$$         Zticks = [3000:250:5500];
% $$$         Pticks = zeros(size(Zticks));
% $$$         for ii=1:length(Zticks)
% $$$             Pticks(ii) = interp1(zofP_mean,P,-Zticks(ii),'linear');
% $$$         end
% $$$         Pticks(1) = 0;
% $$$         set(ax2,'ytick',Pticks);
% $$$         set(ax2,'yticklabel',Zticks);
% $$$         ylabel(ax2,'Depth (m)');
% $$$         ylim(ax2,zlim);
% $$$         set(ax2,'ydir','reverse');
% $$$     end
% $$$ 
% $$$     subplot(3,1,2);
% $$$     [X,Y] = ndgrid(lags/12,T_YvPar);
% $$$     contourf(X,Y,TvPar_lr',[-1e50 Tcxs(1):Tcxs(2):Tcxs(3) 1e50],'linestyle','none');
% $$$     caxis([Tcxs(1) Tcxs(3)]);
% $$$     xlabel('Lag (years)');
% $$$     ylabel(tlab);
% $$$     title(['ACCESS-CM2 PI-control $\Theta(p_\Theta)$' label ' regression']);
% $$$     cb = colorbar;
% $$$     ylabel(cb,'Temperature Anomalies ($^\circ$C / $\sigma$)','Interpreter','latex');
% $$$     colormap(redblue);
% $$$     ylim(tlim);
% $$$     if (~remap_to_T)
% $$$         set(gca,'ydir','reverse');
% $$$         ax1=gca;
% $$$         ax1_pos = ax1.Position; % position of first axes
% $$$         ax1_pos(1) = ax1_pos(1)*0.7;
% $$$         ax1_pos(3) = 0.00001;
% $$$         ax2 = axes('Position',ax1_pos,...
% $$$                    'Color','none');
% $$$         set(ax2,'FontSize',15);
% $$$         Zticks = [30 18 12:-2:-2];
% $$$ % $$$         Zticks = [30 26 22 18:-2:10];
% $$$ % $$$         Zticks = [30 22 18 16:-2:8];
% $$$ % $$$         Zticks = [2:-0.25:-0.5 -2];
% $$$         Pticks = zeros(size(Zticks));
% $$$         for ii=1:length(Zticks)
% $$$             Pticks(ii) = interp1(TvP.Tp_mean,P,Zticks(ii),'linear');
% $$$         end
% $$$         set(ax2,'ytick',Pticks);
% $$$         set(ax2,'yticklabel',Zticks);
% $$$         ylabel(ax2,'Temperature ($^\circ$C)');
% $$$         ylim(ax2,zlim);
% $$$         set(ax2,'ydir','reverse');
% $$$     end
% $$$ 
% $$$     subplot(3,1,3);
% $$$     [X,Y] = ndgrid(lags/12,y_YvPar);
% $$$     contourf(X,Y,yvar_lr',[-1e50 Ycxs(1):Ycxs(2):Ycxs(3) 1e50],'linestyle','none');
% $$$     caxis([Ycxs(1) Ycxs(3)]);
% $$$     xlabel('Lag (years)');
% $$$     ylabel(ylab);
% $$$     ylim(yylim);
% $$$     title(['ACCESS-CM2 PI-control $\Theta(p_\phi)$' label ' regression']);
% $$$     cb = colorbar;
% $$$     ylabel(cb,'Temperature Anomalies ($^\circ$C / $\sigma$)','Interpreter','latex');
% $$$     if (~remap_to_T)
% $$$         ax1=gca;
% $$$         ax1_pos = ax1.Position; % position of first axes
% $$$         ax1_pos(1) = ax1_pos(1)*0.7;
% $$$         ax1_pos(3) = 0.00001;
% $$$         ax2 = axes('Position',ax1_pos,...
% $$$                    'Color','none');
% $$$         set(ax2,'FontSize',15);
% $$$         Zticks = [-75:15:75];
% $$$ % $$$         Zticks = [-70:10:80];
% $$$         Pticks = zeros(size(Zticks));
% $$$         for ii=1:length(Zticks)
% $$$             Pticks(ii) = interp1(yofP_mean,P,Zticks(ii),'linear');
% $$$         end
% $$$         Pticks(1) = 0;
% $$$         set(ax2,'ytick',Pticks);
% $$$         set(ax2,'yticklabel',Zticks);
% $$$         ylabel(ax2,'Latitude ($^\circ$N)');
% $$$         ylim(ax2,yylim);
% $$$     end
    
    
% $$$     %%% Area plots:
% $$$     figure;
% $$$     subplot(2,3,1);
% $$$     plot(TvP.Tap_mean,P,'linewidth',2);
% $$$     xlabel('Temperature ($^\circ$C)');
% $$$     xlim([-2 34]);
% $$$     ylabel('Surface Area Percentile');
% $$$     ylim([0 100]);
% $$$     title('ACCESS-CM2 PI-control mean $\Theta(p_\Theta)$');
% $$$     set(gca,'ydir','reverse')
% $$$ 
% $$$     subplot(2,3,2);
% $$$     [X,Y] = ndgrid(1:12,T_YvPar);
% $$$     contourf(X,Y,TvP.Tap_clim',[-10 -1:0.01:1 10],'linestyle','none');
% $$$     xlabel('Month');
% $$$     ylabel(tlab);
% $$$     ylim(tlim);
% $$$     title('ACCESS-CM2 PI-control $\Theta(p_\Theta)$ seasonal cycle');
% $$$     caxis([-1 1]);
% $$$     colorbar;
% $$$     if (~remap_to_T);set(gca,'ydir','reverse');end;
% $$$     colormap(redblue);
% $$$ 
% $$$     subplot(2,3,3);
% $$$     plot(std(TvP.Tap_clim,[],2),T_YvPar,'linewidth',2);
% $$$     hold on;
% $$$     plot(std(TvP.Tap,[],2),T_YvPar,'-r','linewidth',2);
% $$$     xlabel('std(Temperature) ($^\circ$C)');
% $$$     legend('Seasonal Climatology','Anomalies');
% $$$     ylabel(tlab);
% $$$     ylim(tlim);
% $$$     xlim([0 0.75]);
% $$$     title('ACCESS-CM2 PI-control variability $\Theta(p_\Theta)$');
% $$$     if (~remap_to_T);set(gca,'ydir','reverse');end;
% $$$ 
% $$$     %%%%%% Plot all anomalies time series:
% $$$ 
% $$$     xlims = [yr1 yr1+nyrs];
% $$$     xlims = [300 400];
% $$$ 
% $$$     [tmp t1] = min(abs(time-xlims(1)));
% $$$     [tmp t2] = min(abs(time-xlims(2)));
% $$$     [X,Y] = ndgrid(time(t1:t2),P);
% $$$     subplot(2,3,[4 5 6]);
% $$$     contourf(X,Y,TvP.Tap(:,t1:t2)',[-1e50 -1:0.02:1 1e50],'linestyle','none');
% $$$     caxis([-1 1]);
% $$$     ylim(tlim);
% $$$     xlim(xlims);
% $$$     colormap(redblue);
% $$$     xlabel('Year');
% $$$     ylabel(tlab);
% $$$     cb = colorbar;
% $$$     ylabel(cb,'Temperature Anomalies ($^\circ$C)','Interpreter','latex');
% $$$     title(['ACCESS-CM2 PI-control ' Tlab ' anomalies']);
% $$$     if (~remap_to_T)
% $$$         set(gca,'ydir','reverse');
% $$$         ax1=gca;
% $$$         ax1_pos = ax1.Position; % position of first axes
% $$$         ax1_pos(1) = ax1_pos(1)*0.7;
% $$$         ax1_pos(3) = 0.00001;
% $$$         ax2 = axes('Position',ax1_pos,...
% $$$                    'Color','none');
% $$$         set(ax2,'FontSize',15);
% $$$         Zticks = [30 18 12:-2:-2];
% $$$         Pticks = zeros(size(Zticks));
% $$$         for ii=1:length(Zticks)
% $$$             Pticks(ii) = interp1(TvP.Tp_mean,P,Zticks(ii),'linear');
% $$$         end
% $$$         set(ax2,'ytick',Pticks);
% $$$         set(ax2,'yticklabel',Zticks);
% $$$         ylabel(ax2,'Temperature ($^\circ$C)');
% $$$         ylim(ax2,tlim);
% $$$         set(ax2,'ydir','reverse');
% $$$     end
% $$$ 
% $$$     %%%%% Spectral analysis:
% $$$     fs = 12;
% $$$     [pxx,f] = pmtm(ZvPar(15,:)',[],[],fs);
% $$$     fL = length(f);
% $$$     TvPar_spec = zeros(fL,PL);
% $$$     for pi=1:PL
% $$$         [TvPar_spec(:,pi),~] = pmtm(TvP.Tap(pi,:)',[],[],fs);
% $$$     end
% $$$     
% $$$     figure;
% $$$     set(gcf,'Position',[1 41 2560 1330]);
% $$$     subplot(2,1,1);
% $$$     [X,Y] = ndgrid(f,T_YvPar);
% $$$     contourf(X,Y,log10(TvPar_spec),[-4.5:0.25:-0.25],'linestyle','none');
% $$$     xlim([1/250 1]);    
% $$$     xlab = [300 100 50 10 7 5 3 2 1];
% $$$     ylim(tlim);
% $$$     set(gca,'xscale','log');
% $$$     set(gca,'xtick',1./xlab);
% $$$     set(gca,'xticklabel',xlab);
% $$$     xlabel('Period (years)');
% $$$     caxis([-3 -1]);
% $$$     ylabel(tlab);
% $$$     title(['Multi-taper power spectra of ACCESS-CM2 PI-control ' ...
% $$$            '$\Theta(p_\Theta)$ anomalies ($^\circ$C$^2$/year)']);
% $$$     colorbar;
% $$$     cb = colorbar;
% $$$     set(cb,'ytick',[-4:1:-1]);
% $$$     set(cb,'yticklabel',10.^[-4:1:-1]);
% $$$     if (~remap_to_T)
% $$$         set(gca,'ydir','reverse');
% $$$         ax1=gca;
% $$$         ax1_pos = ax1.Position; % position of first axes
% $$$         ax1_pos(1) = ax1_pos(1)*0.7;
% $$$         ax1_pos(3) = 0.00001;
% $$$         ax2 = axes('Position',ax1_pos,...
% $$$                    'Color','none');
% $$$         set(ax2,'FontSize',15);
% $$$         Zticks = [30 18 12:-2:-2];
% $$$         Pticks = zeros(size(Zticks));
% $$$         for ii=1:length(Zticks)
% $$$             Pticks(ii) = interp1(TvP.Tp_mean,P,Zticks(ii),'linear');
% $$$         end
% $$$         set(ax2,'ytick',Pticks);
% $$$         set(ax2,'yticklabel',Zticks);
% $$$         ylabel(ax2,'Temperature ($^\circ$C)');
% $$$         ylim(ax2,zlim);
% $$$         set(ax2,'ydir','reverse');
% $$$     end
% $$$     colormap(cmocean('dense'));


% $$$ end

% $$$ % term options:
% $$$ haveRedi = 1; % 1 = Redi diffusion is on, 0 = off
% $$$ haveMDS = 1; % 1 = MDS is on, 0 = off;
% $$$ haveSIG = 1; % 1 = SIG is on, 0 = off;
% $$$ 
% $$$ % files:
% $$$ filenames = dir(baseL);
% $$$ 
% $$$ for fi=1:length(filenames)
% $$$     fname = [baseL filenames(fi).name];
% $$$     if (length(strfind(fname,'ocean_month.nc'))>0)
% $$$         dstr = fname(end-8:end);
% $$$ 
% $$$         % grids:
% $$$         area = ncread(fname,'area_t');
% $$$         time = ncread(fname,'time');
% $$$         [xL,yL] = size(area);
% $$$         tL = length(time);
% $$$ 
% $$$ GWBmon.SWH    = zeros(TL+1,tL);GWBmon.VDS    = zeros(TL+1,tL);
% $$$ GWBmon.RMX    = zeros(TL+1,tL);GWBmon.PME    = zeros(TL+1,tL);
% $$$ GWBmon.FRZ    = zeros(TL+1,tL);GWBmon.ETS    = zeros(TL+1,tL);
% $$$ GWBmon.VDF    = zeros(TL+1,tL);GWBmon.KNL    = zeros(TL+1,tL);
% $$$ GWBmon.TEN    = zeros(TL+1,tL);
% $$$ if (haveRedi)
% $$$ GWBmon.K33    = zeros(TL+1,tL);
% $$$ GWBmon.RED    = zeros(TL+1,tL);
% $$$ end
% $$$ if (haveMDS)
% $$$ GWBmon.MDS   = zeros(TL+1,tL);
% $$$ end
% $$$ if (haveSIG)
% $$$ GWBmon.SIG   = zeros(TL+1,tL);
% $$$ end
% $$$ 
% $$$ for zi=1:zL
% $$$     for ti=1:tL
% $$$         sprintf('Calculating MON/AN binned time %03d of %03d, depth %02d of %02d',ti,tL,zi,zL)
% $$$ 
% $$$         temp = ncread(fname,'temp',[1 1 zi ti],[xL yL 1 1]);
% $$$         if (max(max(temp))>120);temp = temp-273.15;end;
% $$$         temp(abs(temp)>100) = NaN;
% $$$         
% $$$         TEN = area.*ncread(fname,'temp_tendency',[1 1 zi ti],[xL yL 1 1]);
% $$$         RMX = area.*ncread(fname,'temp_rivermix',[1 1 zi ti],[xL yL 1 1]);
% $$$         VDS = area.*ncread(fname,'temp_vdiffuse_sbc',[1 1 zi ti],[xL yL 1 1]);
% $$$         SWH = area.*ncread(fname,'sw_heat',[1 1 zi ti],[xL yL 1 1]);
% $$$         VDF = area.*ncread(fname,'temp_vdiffuse_diff_cbt',[1 1 zi ti],[xL yL 1 1]);
% $$$         KNL = area.*ncread(fname,'temp_nonlocal_KPP',[1 1 zi ti],[xL yL 1 1]);
% $$$         FRZ = area.*ncread(fname,'frazil_3d',[1 1 zi ti],[xL yL 1 1]);
% $$$         if (haveRedi)
% $$$             K33 = area.*ncread(fname,'temp_vdiffuse_k33',[1 1 zi ti],[xL yL 1 1]);
% $$$             RED = area.*ncread(fname,'neutral_diffusion_temp',[1 1 zi ti],[xL yL 1 1]);
% $$$         end
% $$$         if (haveMDS)
% $$$             MDS = area.*ncread(fname,'mixdownslope_temp',[1 1 zi ti],[xL yL 1 1]);
% $$$         end
% $$$         if (haveSIG)
% $$$             SIG = area.*ncread(fname,'temp_sigma_diff',[1 1 zi ti],[xL yL 1 1]);
% $$$         end
% $$$ 
% $$$         if (zi == 1)
% $$$             ETS = area.*ncread(fname,'temp_eta_smooth',[1 1 ti],[xL yL 1]);
% $$$             PME = area.*ncread(fname,'sfc_hflux_pme',[1 1 ti],[xL yL 1]);
% $$$         end
% $$$                 
% $$$         for Ti=1:TL
% $$$             %Accumulate sums:
% $$$             inds = find(temp>=Te(Ti) & temp<Te(Ti+1));
% $$$             GWBmon.TEN(Ti,ti) = GWBmon.TEN(Ti,ti)+nansum(TEN(inds));
% $$$             GWBmon.RMX(Ti,ti) = GWBmon.RMX(Ti,ti)+nansum(RMX(inds));
% $$$             GWBmon.VDS(Ti,ti) = GWBmon.VDS(Ti,ti)+nansum(VDS(inds));
% $$$             GWBmon.SWH(Ti,ti) = GWBmon.SWH(Ti,ti)+nansum(SWH(inds));
% $$$             GWBmon.VDF(Ti,ti) = GWBmon.VDF(Ti,ti)+nansum(VDF(inds));
% $$$             GWBmon.KNL(Ti,ti) = GWBmon.KNL(Ti,ti)+nansum(KNL(inds));
% $$$             GWBmon.FRZ(Ti,ti) = GWBmon.FRZ(Ti,ti)+nansum(FRZ(inds));
% $$$             if (haveRedi)
% $$$                 GWBmon.K33(Ti,ti) = GWBmon.K33(Ti,ti)+nansum(K33(inds));
% $$$                 GWBmon.RED(Ti,ti) = GWBmon.RED(Ti,ti)+nansum(RED(inds));
% $$$             end
% $$$             if (haveMDS)
% $$$                 GWBmon.MDS(Ti,ti) = GWBmon.MDS(Ti,ti)+nansum(MDS(inds));
% $$$             end
% $$$             if (haveSIG)
% $$$                 GWBmon.SIG(Ti,ti) = GWBmon.SIG(Ti,ti)+nansum(SIG(inds));
% $$$             end
% $$$             
% $$$             if (zi == 1)
% $$$                 GWBmon.ETS(Ti,ti) = GWBmon.ETS(Ti,ti)+nansum(ETS(inds));
% $$$                 GWBmon.PME(Ti,ti) = GWBmon.PME(Ti,ti)+nansum(PME(inds));
% $$$             end
% $$$         end
% $$$         inds = find(temp>=Te(TL+1));
% $$$         GWBmon.TEN(TL+1,ti) = GWBmon.TEN(TL+1,ti)+nansum(TEN(inds));
% $$$         GWBmon.RMX(TL+1,ti) = GWBmon.RMX(TL+1,ti)+nansum(RMX(inds));
% $$$         GWBmon.VDS(TL+1,ti) = GWBmon.VDS(TL+1,ti)+nansum(VDS(inds));
% $$$         GWBmon.SWH(TL+1,ti) = GWBmon.SWH(TL+1,ti)+nansum(SWH(inds));
% $$$         GWBmon.VDF(TL+1,ti) = GWBmon.VDF(TL+1,ti)+nansum(VDF(inds));
% $$$         GWBmon.KNL(TL+1,ti) = GWBmon.KNL(TL+1,ti)+nansum(KNL(inds));
% $$$         GWBmon.FRZ(TL+1,ti) = GWBmon.FRZ(TL+1,ti)+nansum(FRZ(inds));
% $$$         if (haveRedi)
% $$$             GWBmon.K33(TL+1,ti) = GWBmon.K33(TL+1,ti)+nansum(K33(inds));
% $$$             GWBmon.RED(TL+1,ti) = GWBmon.RED(TL+1,ti)+nansum(RED(inds));
% $$$         end
% $$$         if (haveMDS)
% $$$             GWBmon.MDS(TL+1,ti) = GWBmon.MDS(TL+1,ti)+nansum(MDS(inds));
% $$$         end
% $$$         if (haveSIG)
% $$$             GWBmon.SIG(TL+1,ti) = GWBmon.SIG(TL+1,ti)+nansum(SIG(inds));
% $$$         end
% $$$ 
% $$$         if (zi == 1)
% $$$             GWBmon.ETS(TL+1,ti) = GWBmon.ETS(TL+1,ti)+nansum(ETS(inds));
% $$$             GWBmon.PME(TL+1,ti) = GWBmon.PME(TL+1,ti)+nansum(PME(inds));
% $$$         end
% $$$     end
% $$$ end
% $$$ 
% $$$ % Integrate to get to T'>T:
% $$$ names = fieldnames(GWBmon);
% $$$ for i=1:length(names)
% $$$     for ti=1:tL
% $$$         eval(['GWBmon.' names{i} '(:,ti) = flipud(cumsum(flipud(GWBmon.' ...
% $$$               names{i} '(:,ti))));']);
% $$$     end
% $$$ end
% $$$ 
% $$$ save([outD model dstr '_GlobalHBud.mat'],'GWBmon','T','Te','-v7.3');
% $$$     
% $$$     end
% $$$ end
% $$$
    
    %%% BUDGET TESTING:
    % Fix advection terms by residual:
    vars = {'temp_submeso', ...
             'temp_vdiffuse_diff_cbt', 'temp_nonlocal_KPP', ...
             'temp_vdiffuse_sbc','frazil_3d','sw_heat','temp_rivermix', ...
             'neutral_diffusion_temp','neutral_gm_temp', ...
             'temp_vdiffuse_k33', 'mixdownslope_temp', ...
             'temp_sigma_diff','sfc_hflux_pme','temp_eta_smooth'};    
    temp_advectionT = temp_tendencyT;
    temp_advectionz = temp_tendencyz;
    temp_advectiony = temp_tendencyy;
    for vi=1:length(vars)
        eval(['temp_advectionT = temp_advectionT - ' vars{vi} 'T;']);
        eval(['temp_advectionz = temp_advectionz - ' vars{vi} 'z;']);
        eval(['temp_advectiony = temp_advectiony - ' vars{vi} 'y;']);
    end
    
    % group vars:
    typs = {'T','z','y'};
    for vi = 1:length(typs)
        eval(['TEN' typs{vi} ' = temp_tendency' typs{vi} ';']);
        eval(['ADV' typs{vi} ' = temp_advection' typs{vi} '+temp_submeso' typs{vi} '+neutral_gm_temp' typs{vi} ';']);
        eval(['FOR' typs{vi} ' = temp_vdiffuse_sbc' typs{vi} '+frazil_3d' typs{vi} '+sw_heat' typs{vi} ...
              '+temp_rivermix' typs{vi} '+sfc_hflux_pme' typs{vi} ...
              '+temp_eta_smooth' typs{vi} ';']);
        eval(['RMIX' typs{vi} ' = neutral_diffusion_temp' typs{vi} '+temp_vdiffuse_k33' typs{vi} '+mixdownslope_temp' typs{vi} ...
              '+temp_sigma_diff' typs{vi} ';']);
        eval(['VMIX' typs{vi} ' = temp_vdiffuse_diff_cbt' typs{vi} '+temp_nonlocal_KPP' typs{vi} ';']);
    end
% $$$     ADV = temp_advection  
% $$$     
% $$$     gvars = {{'Tendency',{'temp_tendency'}},...
% $$$              {'Advection',{'temp_advection','temp_submeso','neutral_gm_temp'}}, ...
% $$$              {'Forcing',{'temp_vdiffuse_sbc','frazil_3d','sw_heat','temp_rivermix','sfc_hflux_pme','temp_eta_smooth'}}, ...
% $$$              {'Redi Mixing',{'neutral_diffusion_temp','temp_vdiffuse_k33','mixdownslope_temp','temp_sigma_diff'}}, ...
% $$$              {'Vertical Mixing',{'temp_vdiffuse_diff_cbt','temp_nonlocal_KPP'}}};
    vars = {'TEN','ADV','FOR','RMIX','VMIX'};
    colors = {'m','b','k','r',[0 0.5 0]};     
        
    figure;
    subplot(1,3,1);
    for gi=1:length(vars)
        eval(['var = mean(' vars{gi} 'z,2);']);
        plot(cumsum(var,'reverse')/1e15,Z,'-','color',colors{gi});
        hold on;
    end
    ylabel('Depth (m)');
    set(gca,'ydir','reverse');
    ylim([0 500]);
    xlim([-2 2]);
    xlabel('Vertical heat transport (PW)');
    
    subplot(1,3,2);
    for gi=1:length(vars)
        eval(['var = mean(' vars{gi} 'T,2);']);
        plot(cumsum(var,'reverse')/1e15,T,'-','color',colors{gi});
        hold on;
    end
    legend(vars);
    ylabel('Temperature ($^\circ$C)');
    xlabel('Diathermal heat transport (PW)');

    subplot(1,3,3);
    for gi=1:length(vars)
        eval(['var = mean(' vars{gi} 'y,2);']);
        plot(cumsum(var)/1e15,latv,'-','color',colors{gi});
        hold on;
    end
    ylabel('Latitude ($^\circ$N)');
    xlabel('Meridional heat transport (PW)');
    
    
    



end

%% Budget variables summary:


% $$$ # Surface heat fluxes (not including surface volume flux terms):
% $$$ SFCH = temp_vdiffuse_sbc + sw_heat + frazil_3d + # 3D vars
% $$$       temp_eta_smooth; # 2D vars
% $$$       
% $$$ # Surface heat fluxes from surface volume fluxes
% $$$ SFCV = temp_rivermix + # 3D vars
% $$$        sfc_hflux_pme # 2D vars
% $$$        
% $$$ # Shortwave redistribution
% $$$ SWR = sw_heat # 3D vars
% $$$ 
% $$$ # Vertical mixing
% $$$ VMIX = temp_vdiffuse_diff_cbt + temp_nonlocal_KPP # 3D vars
% $$$ 
% $$$ # Miscellaneous mixing
% $$$ SMIX = mixdownslope_temp + temp_sigma_diff # 3D vars
% $$$ 
% $$$ # Neutral diffusion
% $$$ RMIX = temp_vdiffuse_k33 + neutral_diffusion_temp # 3D vars
% $$$ 
% $$$ # Total tendency
% $$$ TEN = temp_tendency # 3D vars
% $$$ 
% $$$ # Total external surface forcing
% $$$ SFC = SFCH + SFCV
% $$$ 
% $$$ # Total internal surface forcing
% $$$ SFCI = SFC - rho0*Cp*THETA*SVF
% $$$ 
% $$$ # Total explicit mixing
% $$$ MIX = VMIX+SMIX+RMIX
% $$$ 
% $$$ # Numerical mixing (by residual)
% $$$ NMIX = dHI/dt - SFCI - MIX
