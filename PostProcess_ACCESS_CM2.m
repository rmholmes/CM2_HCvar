% Post-processing for HC variability project. This requires processed
% (from Process_ACCESS-CM2.m) and collated (from Collate_ACCESS-CM2)
% output .mat files.

%%%%%% OPTIONS %%%
clear all;
baseMAT = 'D:/DATA/access-cm2/';

% Streamline post-processing:
doBUDGET = 1;
doTyz = 0;
doMOCyz = 0;
doTAU = 0;

load([baseMAT 'CM2_PIcontrolTb05__ALL.mat']);
saveNAME = 'PIcontrolTb05PP_Tint.mat';

saveMAT = 1;

%%%%%%

% Define a new percentile grid:

% linear:
dP = 0.25;
Pe = 0:dP:100;
P = (Pe(2:end)+Pe(1:end-1))/2;
PL = length(P);
dP = repmat(dP,[PL 1]);

% Time vector:
tL = length(time);
time = time/365.25; % time in years

% The temp_advection term in the ACCESS-CM2 output is
% corrupted. Fix this by calculating it by residual:
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
% advection = tendency - submeso - diff_cbt - nonlocal_KPP - sbc -
% frazil - sw_heat - rivermix - neutral_diffusion_temp -
% neutral_gm_temp - vdiffuse_k33 - mixdownslope_temp - sigma_diff -
% hflux_pme - eta_smooth

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

% internal-heat surface volume flux term:
if (doBUDGET)
    % pme_river = kg/m^3 * m/sec * m^2 = kg/sec
    % JSH = pme_river*Theta*Cp = kg/sec * degC * J kg-1 degC-1 = W
    JSH_c = Tv.pme_river_c*Cp.*repmat(Te',[1 tL]);
    Tv.FOR_c = Tv.FOR_c - JSH_c;
    Tv.ADV_c = Tv.ADV_c + JSH_c;
    Tv = rmfield(Tv,'pme_river');Zv = rmfield(Zv,'pme_river');Yv = rmfield(Yv,'pme_river');
    Tv = rmfield(Tv,'pme_river_c');Zv = rmfield(Zv,'pme_river_c');Yv = rmfield(Yv,'pme_river_c');
end

OHC = squeeze(Zv.H_c(end,:))';


%%% Remap to percentile

% Define percentiles using volume:
Vtot = Zv.V_c(end,:)';
Atot = Tv.A_c(1,:)';
Tv.P   = 100*Tv.V_c./repmat(Tv.V_c(1,:),[TL+1 1]);
Zv.P   = 100*Zv.V_c./repmat(Zv.V_c(end,:),[zL+1 1]);
Yv.P   = 100*Yv.V_c./repmat(Yv.V_c(end,:),[yL+1 1]);
Tv.Pa  = 100*Tv.A_c./repmat(Atot',[TL+1 1]);

%%% Interpolate temperatures
    % Calculate z-temperature:
    Zv.T = Zv.H/rho0/Cp./Zv.V;
    Zv.Te = cat(1,Zv.T(1,:),(Zv.T(2:end,:)+Zv.T(1:end-1,:))/2,Zv.T(end,:));

    % Calculate y-temperature:
    Yv.T = Yv.H/rho0/Cp./Yv.V;
    Yv.Te = cat(1,Yv.T(1,:),(Yv.T(2:end,:)+Yv.T(1:end-1,:))/2,Yv.T(end,:));

    % Interpolate variables to P levels:
    for ti = 1:tL
        [Pun,Iun] = unique(Tv.P(:,ti));
        TvP.Tp(:,ti) = interp1(Pun,Te(Iun),P,'linear');
        
        [Pun,Iun] = unique(Tv.Pa(:,ti));
        TvP.Tap(:,ti) = interp1(Pun,Te(Iun),P,'linear');
        
        [Pun,Iun] = unique(Zv.P(:,ti));
        ZvP.Tp(:,ti) = interp1(Pun,Zv.Te(Iun,ti),P,'linear');
        
        [Pun,Iun] = unique(Yv.P(:,ti));
        YvP.Tp(:,ti) = interp1(Pun,Yv.Te(Iun,ti),P,'linear');
        
        [Pun,Iun] = unique(Tv.P(:,ti));
        TvP.A_c(:,ti) =  interp1(Pun,Tv.A_c(Iun,ti),Pe,'linear');
    end
    TvP.A = diff(TvP.A_c,[],1);
end

if (doBUDGET)
    bvars = {'TEN_c','ADV_c','ADVGM_c','FOR_c','RMIX_c','VMIX_c'};
    for ti = 1:tL
        for vi=1:length(bvars)
            [Pun,Iun] = unique(Zv.P(:,ti));
            eval(['ZvP.' bvars{vi} '(:,ti) = interp1(Pun,Zv.' bvars{vi} '(Iun,ti),Pe,''linear'');']);
            [Pun,Iun] = unique(Yv.P(:,ti));
            eval(['YvP.' bvars{vi} '(:,ti) = interp1(Pun,Yv.' bvars{vi} '(Iun,ti),Pe,''linear'');']);
            [Pun,Iun] = unique(Tv.P(:,ti));
            eval(['TvP.' bvars{vi} '(:,ti) = interp1(Pun,Tv.' bvars{vi} '(Iun,ti),Pe,''linear'');']);
        end
    end

    %%%% Pre-anomaly budget plotting:
    plot_budget_ACCESS_CM2;

    %%%% Time-integrate budget terms:
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
    yr1 =50; % initial years to throw out
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
    if (isfield(TvP,'A_c'))
        ZvP.A_c = TvP.A_c;
        YvP.A_c = TvP.A_c; % for convience -> delete after...
    end
    if (isfield(TvP,'A'))
        ZvP.A = TvP.A;
        YvP.A = TvP.A; % for convience -> delete after...
    end
    % Calculate cubic drift from PI-control        
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
    TvP.Hp = rho0*Cp*cumsum(TvP.Tp,1).*repmat(dP,[1 tL])/100.*repmat(Vtot',[PL 1]);
    ZvP.Hp = rho0*Cp*cumsum(ZvP.Tp,1).*repmat(dP,[1 tL])/100.*repmat(Vtot',[PL 1]);
    YvP.Hp = rho0*Cp*cumsum(YvP.Tp,1).*repmat(dP,[1 tL])/100.*repmat(Vtot',[PL 1]);
    
    % Calculate tendency term contributions to temperature:
    if (doBUDGET)
    for ti =1:length(typs)
        for vi = 1:length(bvars)
            tn = typs{ti};
            vn = bvars{vi};
            eval([tn 'vP.' vn(1:end-2) 'p = diff(' tn 'vP.' vn ',[],1)/rho0/' ...
                                'Cp./repmat(dP,[1 tL])*100./repmat(Vtot'',[PL ' ...
                  '1]);']);
            eval([tn 'vP.' vn(1:end-2) 'p_clim = diff(' tn 'vP.' vn '_clim,[],1)/rho0/' ...
                                'Cp./repmat(dP,[1 12])*100./repmat(Vtot_clim'',[PL ' ...
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
        Tyz_mean = mean(Tyz,3);
        
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
        
        saveNAMETyz = [saveNAME(1:end-4) '_Tyz.mat'];
        save([baseMAT saveNAMETyz],'Tyz','Tyz_clim','Tyz_mean');
    end    

    % MOCyz processing:
    if (doMOCyz)
                
        % Load:
        load([baseMAT 'CM2_PIcontrolTb05__MOCyz_ALL.mat']);
        ti = yr1*12+1; % starting month
        MOC   = MOCyzS.MOC;
        MOCgm = MOCyzS.MOCgm;
        clear MOCyzS;
        MOC = cumsum(MOC(:,:,ti:(ti+tL-1)),2)/1e6/rho0;     
        MOCgm = MOCgm(:,:,ti:(ti+tL-1))/1e6/rho0;     
        MOC_mean = mean(MOC,3);
        MOCgm_mean = mean(MOCgm,3);
        
        % Detrend:
        [MOC_cub,MOC_cubtr] = cubfit(time,reshape(MOC,[yL*zL tL])');
        MOC_cub = reshape(MOC_cub',[yL zL tL]);
        MOC_cubtr = reshape(MOC_cubtr',[yL zL 4]);
        MOC = MOC-MOC_cub;
        
        % Deseason:
        MOC_clim = zeros(yL,zL,12);
        for mi=1:12
            MOC_clim(:,:,mi) = mean(MOC(:,:,mi:12:end),3);
        end
        MOC = MOC - repmat(MOC_clim,[1 1 nyrs]);

        % Detrend:
        [MOCgm_cub,MOCgm_cubtr] = cubfit(time,reshape(MOCgm,[yL*zL tL])');
        MOCgm_cub = reshape(MOCgm_cub',[yL zL tL]);
        MOCgm_cubtr = reshape(MOCgm_cubtr',[yL zL 4]);
        MOCgm = MOCgm-MOCgm_cub;
        
        % Deseason:
        MOCgm_clim = zeros(yL,zL,12);
        for mi=1:12
            MOCgm_clim(:,:,mi) = mean(MOCgm(:,:,mi:12:end),3);
        end
        MOCgm = MOCgm - repmat(MOCgm_clim,[1 1 nyrs]);

        saveNAMEMOC = [saveNAME(1:end-4) '_MOC.mat'];
        save([baseMAT saveNAMEMOC],'MOC','MOCgm','MOC_clim','MOCgm_clim','MOC_mean','MOCgm_mean');
    end    

    if (doTAU)
                
        % Load:
        load([baseMAT 'CM2_PIcontrolTb05__taux_ALL.mat']);
        ti = yr1*12+1; % starting month
        taux   = tauxS.tau_x(:,ti:(ti+tL-1));
        taux_mean = mean(taux,2);
        
        % Detrend:
        [taux_cub,taux_cubtr] = cubfit(time,reshape(taux,[yL tL])');
        taux_cub = taux_cub';
        taux_cubtr = taux_cubtr';
        taux = taux-taux_cub;
        
        % Deseason:
        taux_clim = zeros(yL,12);
        for mi=1:12
            taux_clim(:,mi) = mean(taux(:,mi:12:end),2);
        end
        taux = taux - repmat(taux_clim,[1 nyrs]);

        saveNAMEtau = [saveNAME(1:end-4) '_taux.mat'];
        save([baseMAT saveNAMEtau],'taux','taux_mean','taux_clim','taux_cub');
    end    

if (saveMAT)
    save([baseMAT saveNAME]);
end

