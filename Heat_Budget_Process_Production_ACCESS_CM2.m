% Process surface forcing and parameterized mixing terms from
% ACCESS-CM2 runs binned into temperature coordinates.

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

clear all;

plot_only = 1;
PI_or_his = 1; % 1 = PI-control, 0 = historical simualtion
if (PI_or_his)
    mname = 'ACCESS_PIcontrol_TLZvolp.mat';
    mnameBUD = 'ACCESS_PIcontrol_TLZbud.mat';
else
    mname = 'ACCESS_Historical_TLZvolp.mat';
    mnameBUD = 'ACCESS_Historical_TLZbud.mat';
end


%%%%%% DATA PROCESSING

if (~plot_only)

    if (PI_or_his)
        base = '/g/data/p66/cm2704/archive/bi889/history/ocn/';
        name = 'PIcontrol';
        fname = [base 'ocean_month.nc-08500630'];
    else
        base = '/g/data/p66/cm2704/archive/bj594/history/ocn/';
        name = 'historical';
        fname = [base 'ocean_month.nc-18500630'];
    end

    outD = '/scratch/e14/rmh561/access-cm2/';

    % Constants:
    Cp = 3992.10322329649; % J kg-1 degC-1
    rho0 = 1035; % kgm-3

    % Define a temperature grid:
    dT = 0.1;
    Te = -3:dT:34;
    T = (Te(2:end)+Te(1:end-1))/2;
    TL = length(T);
    
    % Constant grid parameters:
    area = ncread(fname,'area_t');
% $$$     lon = ncread(fname,'geolon_t');
% $$$     lat = ncread(fname,'geolat_t');
    lonv = ncread(fname,'xt_ocean');
    latv = ncread(fname,'yt_ocean');
    latu = ncread(fname,'yu_ocean');
    latv_edges = cat(1,latv(1)-(latu(1)-latv(1)),latu);
    
    [xL,yL] = size(area);
    zL = 50;

    % Depth grid:
    Z = ncread(fname,'st_ocean');
    Ze = ncread(fname,'st_edges_ocean');
    zL = length(Z);
    
    % 3D mask:
    temp = ncread(fname,'temp',[1 1 1 1],[xL yL zL 1]);
    mask = ~isnan(temp);
    
    % A(z):
    A = zeros(zL,1);
    for zi=1:zL
        A(zi) = nansum(nansum(area(mask(:,:,zi))));
    end
    
    % Save grid info:
    save(mname,'Cp','rho0','dT','Te','T','TL',...
               'xL','yL','zL','Z','Ze','zL','A', ...
               'latv','latv_edges');
    save(mnameBUD,'Cp','rho0','dT','Te','T','TL',...
               'xL','yL','zL','Z','Ze','zL','A', ...
               'latv','latv_edges');
    
    %%%%%% Standard non-budget variables (mname):
    
    % Initialize variables:
    Hz = []; % Heat content at depth z (J)
    Vz = []; % Volume at depth z (m3)
    HT = []; % Heat content at temperature T (J)
    VT = []; % Volume at temperature T (m3)
    Hy = []; % Heat content at latitude y (J)
    Vy = []; % Volume at latitude y (m3)
    time = []; % time axis
    DT_A = []; % averaging time
    AT = []; % Area at temperature T (m2)

    % SST-based climate indices:
    N34 = []; % Nino 3.4
    TPIr1 = []; % Tripole index (IPO) region 1
    TPIr2 = []; % Tripole index (IPO) region 2
    TPIr3 = []; % Tripole index (IPO) region 3

    [tmp, N34x1] = min(abs(lonv+170));  [tmp, N34x2] = min(abs(lonv+120));
    [tmp, N34y1] = min(abs(latv+5));    [tmp, N34y2] = min(abs(latv-5));
    
    [tmp, TPIr1x1] = min(abs(lonv+220));  [tmp, TPIr1x2] = min(abs(lonv+145));
    [tmp, TPIr1y1] = min(abs(latv-25));    [tmp, TPIr1y2] = min(abs(latv-45));
    [tmp, TPIr2x1] = min(abs(lonv+190));  [tmp, TPIr2x2] = min(abs(lonv+90));
    [tmp, TPIr2y1] = min(abs(latv+10));    [tmp, TPIr2y2] = min(abs(latv-10));
    [tmp, TPIr3x1] = min(abs(lonv+230));  [tmp, TPIr3x2] = min(abs(lonv+160));
    [tmp, TPIr3y1] = min(abs(latv+50));    [tmp, TPIr3y2] = min(abs(latv+15));
    
    % AMOC index:
    AMOC = []; % AMOC index (Sv)
    potrho = ncread(fname,'potrho_edges');
    potrhoL = length(potrho)-1;
    
    [tmp, AMOClt] = min(abs(latv-26)); [tmp, AMOCrho] = min(abs(potrho-1035.5));
    [tmp, AMOCln1] = min(abs(lonv+103));    [tmp, AMOCln2] = min(abs(lonv+5));
    
    % Wind power index:
    WPOW = [];
    
    % Start file loop:
    files = dir(base);

    for fi = 1:length(files)
        if (strfind(files(fi).name,'month'))

            fname = [base files(fi).name];
            sprintf('Doing %03d of %03d',fi,length(files))
            time_t = ncread(fname,'time');
            DT_A_t = ncread(fname,'average_DT')*86400;
        
            time = cat(1,time,time_t);
            DT_A = cat(1,DT_A,DT_A_t);

            tL = length(time_t);
            
            temp = ncread(fname,'temp');
            temp(~mask) = NaN;
            if (max(max(temp))>120); temp=temp-273.15;end;
            V = ncread(fname,'dht').*repmat(area,[1 1 zL tL]);
            V(isnan(V)) = 0;
            H = rho0*Cp*temp.*V;
            SST = squeeze(temp(:,:,1,:));

            % Latitude space H and V:
            Hy_t = squeeze(nansum(nansum(H,1),3));
            Vy_t = squeeze(nansum(nansum(V,1),3));

            Hy = cat(2,Hy,Hy_t);
            Vy = cat(2,Vy,Vy_t);

            % Depth space H and V:
            Hz_t = squeeze(nansum(nansum(H,1),2));
            Vz_t = squeeze(nansum(nansum(V,1),2));
            
            Hz = cat(2,Hz,Hz_t);
            Vz = cat(2,Vz,Vz_t);
            
            % Temp space H and V:
            VT_t = zeros(TL,tL);
            HT_t = zeros(TL,tL);
            AT_t = zeros(TL,tL);
            for Ti=1:TL
                %Accumulate sums:
                inds = temp>=Te(Ti) & temp<Te(Ti+1);
                VT_t(Ti,:) = VT_t(Ti,:) + squeeze(nansum(nansum(nansum(V.*inds,1),2),3))';
                HT_t(Ti,:) = HT_t(Ti,:) + squeeze(nansum(nansum(nansum(H.*inds,1),2),3))';
                indsS = SST>=Te(Ti) & SST<Te(Ti+1);
                AT_t(Ti,:) = AT_t(Ti,:) + squeeze(nansum(nansum(repmat(area,[1 1 tL]).*indsS,1),2))';
            end
            % Account for water warmer than max temperature (possible for CM2):
            inds = temp>Te(end);
            VT_t(end,:) = VT_t(end,:) + squeeze(nansum(nansum(nansum(V.*inds,1),2),3))';
            HT_t(end,:) = HT_t(end,:) + squeeze(nansum(nansum(nansum(H.*inds,1),2),3))';
            indsS = SST>Te(end);
            AT_t(end,:) = AT_t(end,:) + squeeze(nansum(nansum(repmat(area,[1 1 tL]).*indsS,1),2))';
            
            HT = cat(2,HT,HT_t);
            VT = cat(2,VT,VT_t);
            AT = cat(2,AT,AT_t);
            
            % Calculate SST based indices:
            N34_t = squeeze(nanmean(nanmean(temp(N34x1:N34x2,N34y1:N34y2,1,:),1),2));
            N34 = cat(1,N34,N34_t);

            TPIr1_t = squeeze(nanmean(nanmean(temp(TPIr1x1:TPIr1x2,TPIr1y1:TPIr1y2,1,:),1),2));
            TPIr1 = cat(1,TPIr1,TPIr1_t);
            TPIr2_t = squeeze(nanmean(nanmean(temp(TPIr2x1:TPIr2x2,TPIr2y1:TPIr2y2,1,:),1),2));
            TPIr2 = cat(1,TPIr2,TPIr2_t);
            TPIr3_t = squeeze(nanmean(nanmean(temp(TPIr3x1:TPIr3x2,TPIr3y1:TPIr3y2,1,:),1),2));
            TPIr3 = cat(1,TPIr3,TPIr3_t);
            
            % AMOC index:
            AMOC_t = squeeze(nansum( ...
            nansum(ncread(fname,'ty_trans_rho',[AMOCln1 AMOClt AMOCrho 1],[AMOCln2-AMOCln1+1 1 potrhoL-AMOCrho+1 tL]),3) + ...
                   ncread(fname,'ty_trans_rho_gm',[AMOCln1 AMOClt AMOCrho 1],[AMOCln2-AMOCln1+1 1 1 tL]),1));
            AMOC = cat(1,AMOC,-AMOC_t/1e6/rho0);
            
            % WPOW index:
            WPOW_t = squeeze(nansum(nansum(ncread(fname,'wind_power_u')+ncread(fname,'wind_power_v'),1),2));
            WPOW = cat(1,WPOW,WPOW_t);
            
            if (mod(fi,5)==0)
                save(mname,'time','DT_A', ...
                     'Hz','Vz','HT','VT','Hy','Vy','AT', ...
                     'N34','TPIr1','TPIr2','TPIr3','AMOC','WPOW', ...
                     '-append');
            end
        end
    end

    %%%%%% Budget variables (mnameBUD):

else % Plotting

    load(mname);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% POST PROCESSING
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Define a new percentile grid:
    dP = 0.25;
    Pe = 0:dP:100;
    P = (Pe(2:end)+Pe(1:end-1))/2;
    PL = length(P);
    
    tL = length(time);
    time = time/365.25; % time in years
    
    % Do cumulative sums:
    Hz_c = cat(1,zeros(1,tL),cumsum(Hz,1));
    Vz_c = cat(1,zeros(1,tL),cumsum(Vz,1));
    Hy_c = cat(1,zeros(1,tL),cumsum(Hy,1));
    Vy_c = cat(1,zeros(1,tL),cumsum(Vy,1));
    HT_c = cat(1,cumsum(HT,1,'reverse'),zeros(1,tL));
    VT_c = cat(1,cumsum(VT,1,'reverse'),zeros(1,tL));
    AT_c = cat(1,cumsum(AT,1,'reverse'),zeros(1,tL));

    OHC = squeeze(Hz_c(end,:))';

% $$$     %%% New approach: Interpolate HC
% $$$     % Interpolate only high-resolution grids for smoothness:
% $$$     zfac = 20;
% $$$     Hz_c = interp1(Hz_c,1:1/zfac:(zL+1));Vz_c = interp1(Vz_c,1:1/zfac:(zL+1));
% $$$     Tfac = 20;
% $$$     HT_c = interp1(HT_c,1:1/Tfac:(TL+1));VT_c = interp1(VT_c,1:1/Tfac:(TL+1));
% $$$     yfac = 2;
% $$$     Hy_c = interp1(Hy_c,1:1/yfac:(yL+1));Vy_c = interp1(Vy_c,1:1/yfac:(yL+1));
% $$$ 
% $$$ 
% $$$     % Remap to ocean percentile
% $$$     Vtot = Vz_c(end,:)';
% $$$     PT   = 100*VT_c./repmat(Vtot',[TL*Tfac+1 1]);
% $$$     Pz   = 100*Vz_c./repmat(Vtot',[zL*zfac+1 1]);
% $$$     Py   = 100*Vy_c./repmat(Vtot',[yL*yfac+1 1]);
% $$$ 
% $$$     Hzp = zeros(PL+1,tL);
% $$$     Hyp = zeros(PL+1,tL);
% $$$     HTp = zeros(PL+1,tL);
% $$$     for ti = 1:tL
% $$$         Hzp(:,ti) = interp1((Pz(:,ti)+(1:zL*zfac+1)'/1e10),Hz_c(:,ti),Pe,'linear');
% $$$         Hyp(:,ti) = interp1((Py(:,ti)+(1:yL*yfac+1)'/1e10),Hy_c(:,ti),Pe,'linear');
% $$$         HTp(:,ti) = interp1((PT(:,ti)+(1:TL*Tfac+1)'/1e10),HT_c(:,ti),Pe,'linear');
% $$$     end
% $$$     Hzp(1,:) = 0;
% $$$     Hyp(1,:) = 0;
% $$$     HTp(1,:) = 0;
% $$$ 
% $$$     % Calculate temperatures from HC:
% $$$     Tzp = (Hzp(2:end,:)-Hzp(1:(end-1),:))/rho0/Cp*100/dP./repmat(Vtot',[PL 1]);
% $$$     Typ = (Hyp(2:end,:)-Hyp(1:(end-1),:))/rho0/Cp*100/dP./repmat(Vtot',[PL 1]);
% $$$     TTp = (HTp(2:end,:)-HTp(1:(end-1),:))/rho0/Cp*100/dP./repmat(Vtot',[PL 1]);
% $$$ 
% $$$     zofP_mean = -interp1(mean(Pz(1:zfac:end,:),2),Ze,P,'linear'); % Depth axis to
% $$$                                                     % go with depth-volume
% $$$     yofP_mean = interp1(mean(Py(1:yfac:end,:),2)+(1:yL+1)'/1e10,latv_edges,P,'linear'); % Depth axis to
% $$$                                                     % go with depth-volume

    %%% Old approach: Interpolate temperatures
    % Remap to ocean percentile
    Vtot = Vz_c(end,:)';
    Atot = AT_c(1,:)';
    PT   = 100*VT_c./repmat(Vtot',[TL+1 1]);
    Pz   = 100*Vz_c./repmat(Vtot',[zL+1 1]);
    Py   = 100*Vy_c./repmat(Vtot',[yL+1 1]);
    PTa  = 100*AT_c./repmat(Atot',[TL+1 1]);

    % Calculate z-temperature:
    Tz = Hz/rho0/Cp./Vz;
    Tze = cat(1,Tz(1,:),(Tz(2:end,:)+Tz(1:end-1,:))/2,Tz(end,:));

    % Calculate y-temperature:
    Ty = Hy/rho0/Cp./Vy;
    Tye = cat(1,Ty(1,:),(Ty(2:end,:)+Ty(1:end-1,:))/2,Ty(end,:));

    TTp = zeros(PL,tL);
    TTap = zeros(PL,tL);
    Tzp = zeros(PL,tL);
    Typ = zeros(PL,tL);
    for ti = 1:tL
        TTp(:,ti) = interp1(PT(:,ti)+(1:TL+1)'/1e10,Te,P,'linear');
        TTp(1,ti) = Te(end);
        TTap(:,ti) = interp1(PTa(:,ti)+(1:TL+1)'/1e10,Te,P,'linear');
        TTap(1,ti) = Te(end);
        Tzp(:,ti) = interp1(Pz(:,ti)+(1:zL+1)'/1e10,Tze(:,ti),P,'linear');
        Typ(:,ti) = interp1(Py(:,ti)+(1:yL+1)'/1e10,Tye(:,ti),P,'linear');
    end

    zofP_mean = -interp1(mean(Pz,2),Ze,P,'linear'); % Depth axis to
                                                    % go with depth-volume
    yofP_mean = interp1(mean(Py,2)+(1:yL+1)'/1e10,latv_edges,P,'linear'); % Depth axis to
                                                    % go with depth-volume

    % Construct climatology and make years even:
    if (PI_or_his)
        yr1 =50; % initial years to throw out
    else
        yr1 = 0;
    end
    nyrs = floor(tL/12)-yr1; % number of years
    tL = nyrs*12; % number of months
    ti = yr1*12+1; % starting month
    TTp = TTp(:,ti:(ti+tL-1));
    TTap = TTap(:,ti:(ti+tL-1));
    Tzp = Tzp(:,ti:(ti+tL-1)); 
    Typ = Typ(:,ti:(ti+tL-1)); 
    
    TIMESERIES = {'N34','TPIr1','TPIr2','TPIr3','OHC','AMOC','WPOW'};
    for vi=1:length(TIMESERIES)
        eval([TIMESERIES{vi} ' = ' TIMESERIES{vi} '(ti:(ti+tL-1));']);
    end
    time = time(ti:(ti+tL-1));
    Vtot = Vtot(ti:(ti+tL-1));
    AMOCfull = AMOC;
    OHCfull = OHC;
        
    % Calculate cubic drift trend and subtract:
    if (PI_or_his)
        
    % Calculate means:
    TTp_mean = mean(TTp,2);
    TTap_mean = mean(TTap,2);
    Tzp_mean = mean(Tzp,2);
    Typ_mean = mean(Typ,2);    

    TTp_cubtr = zeros(4,PL);
    TTp_cub = zeros(tL,PL);
    TTap_cubtr = zeros(4,PL);
    TTap_cub = zeros(tL,PL);
    Tzp_cubtr = zeros(4,PL);
    Tzp_cub = zeros(tL,PL);
    Typ_cubtr = zeros(4,PL);
    Typ_cub = zeros(tL,PL);
    for pi=1:PL
        [Tzp_cub(:,pi),Tzp_cubtr(:,pi)] = cubfit(time,Tzp(pi,:)');
        [Typ_cub(:,pi),Typ_cubtr(:,pi)] = cubfit(time,Typ(pi,:)');
        [TTp_cub(:,pi),TTp_cubtr(:,pi)] = cubfit(time,TTp(pi,:)');
        [TTap_cub(:,pi),TTap_cubtr(:,pi)] = cubfit(time,TTap(pi,:)');
    end
    TTp_cub = TTp_cub';
    TTap_cub = TTap_cub';
    Tzp_cub = Tzp_cub';
    Typ_cub = Typ_cub';
    [OHC_cub,OHC_cubtr] = cubfit(time,OHC);
        
    % subtract cubear trend:
    TTp = TTp-TTp_cub;
    Tzp = Tzp-Tzp_cub;
    Typ = Typ-Typ_cub;
    OHC = OHC-OHC_cub;
    
    PI.Tzp_cubtr = Tzp_cubtr;
    PI.TTp_cubtr = TTp_cubtr;
    PI.Typ_cubtr = Typ_cubtr;
    PI.time = time;
    PI.OHCfull = OHCfull;
    PI.OHC = OHC;
    PI.OHC_cubtr = OHC_cubtr;
    PI.TTp_mean = TTp_mean;
    PI.Typ_mean = Typ_mean;
    PI.Tzp_mean = Tzp_mean;
    save('PIcubicFit.mat','PI');

    else
        Itime = -100;
        load('PIcubicFit.mat');
        TTp_mean = PI.TTp_mean;
        Typ_mean = PI.Typ_mean;
        Tzp_mean = PI.Tzp_mean;

        TTp_cub = zeros(tL,PL);
        TTap_cub = zeros(tL,PL);
        Tzp_cub = zeros(tL,PL);
        Typ_cub = zeros(tL,PL);
        
        t = time-Itime;
        t2 = (time-Itime).^2;
        t3 = (time-Itime).^3;
        
        for pi = 1:PL
            TTp_cub(:,pi) = PI.TTp_cubtr(1,pi)+PI.TTp_cubtr(2,pi).*t ...
                          + PI.TTp_cubtr(3,pi).*t2 + PI.TTp_cubtr(4,pi).*t3;
            Typ_cub(:,pi) = PI.Typ_cubtr(1,pi)+PI.Typ_cubtr(2,pi).*t ...
                          + PI.Typ_cubtr(3,pi).*t2 + PI.Typ_cubtr(4,pi).*t3;
            Tzp_cub(:,pi) = PI.Tzp_cubtr(1,pi)+PI.Tzp_cubtr(2,pi).*t ...
                          + PI.Tzp_cubtr(3,pi).*t2 + PI.Tzp_cubtr(4,pi).*t3;
        end
        OHC_cub = PI.OHC_cubtr(1)+PI.OHC_cubtr(2).*t + ...
                  PI.OHC_cubtr(3).*t2 + PI.OHC_cubtr(4).*t3;
    
        % subtract cubear trend:
        TTp = TTp-TTp_cub';
        Tzp = Tzp-Tzp_cub';
        Typ = Typ-Typ_cub';
        OHC = OHC-OHC_cub;
        
        % Set date:
        time = time+1850;
        yr1 = 1850;
    end
        
    % Seasonal cycle:
    TTp_clim = zeros(PL,12);
    TTap_clim = zeros(PL,12);
    Tzp_clim = zeros(PL,12);
    Typ_clim = zeros(PL,12);
    for vi=1:length(TIMESERIES)
        eval([TIMESERIES{vi} '_clim = zeros(12,1);']);
    end
    for mi=1:12
        TTp_clim(:,mi) = mean(TTp(:,mi:12:end),2);
        TTap_clim(:,mi) = mean(TTap(:,mi:12:end),2);
        Tzp_clim(:,mi) = mean(Tzp(:,mi:12:end),2);
        Typ_clim(:,mi) = mean(Typ(:,mi:12:end),2);
        for vi=1:length(TIMESERIES)
            eval([TIMESERIES{vi} '_clim(mi) = mean(' TIMESERIES{vi} ...
                  '(mi:12:end));']);
        end
    end
    
    % Subtract mean of season to avoid changing mean:
    TTp_clim = TTp_clim-repmat(mean(TTp_clim,2),[1 12]);
    Tzp_clim = Tzp_clim-repmat(mean(Tzp_clim,2),[1 12]);
    Typ_clim = Typ_clim-repmat(mean(Typ_clim,2),[1 12]);
    for vi=1:length(TIMESERIES)
        eval([TIMESERIES{vi} '_clim = ' TIMESERIES{vi} '_clim - mean(' TIMESERIES{vi} ');']);
    end
    
    % Subtract seasonal cycle:
    TTp = TTp - repmat(TTp_clim,[1 nyrs]);
    TTap = TTap - repmat(TTap_clim,[1 nyrs]);
    Tzp = Tzp - repmat(Tzp_clim,[1 nyrs]);
    Typ = Typ - repmat(Typ_clim,[1 nyrs]);
    for vi=1:length(TIMESERIES)
        eval([TIMESERIES{vi} ' = ' TIMESERIES{vi} ...
              ' - repmat(' TIMESERIES{vi} '_clim,[nyrs 1]);']);
    end
    
    if (PI_or_his)
        PI.OHC = OHC;
        save('PIcubicFit.mat','PI');
    end
    
    % Calculate cumulative integral OHC (overwrite total - just anomalies now):
    HTp = rho0*Cp*cumsum(TTp,1)*dP/100.*repmat(Vtot',[PL 1]);
    Hzp = rho0*Cp*cumsum(Tzp,1)*dP/100.*repmat(Vtot',[PL 1]);
    Hyp = rho0*Cp*cumsum(Typ,1)*dP/100.*repmat(Vtot',[PL 1]);
    
    % Calc GMSST:
    GMSST = Tzp(1,:)';
    
    % Calc TPI:
    TPIraw = TPIr2-(TPIr1+TPIr3)/2;
    [b,a] = cheby1(6,10,2/13/12);
    TPI = filter(b,a,TPIraw);
    % Shift time (not sure why I have to do this???):
    TPI(1:(end-13*12)) = TPI(13*12:end-1);
    TPI((end-13*12+1):end) = 0;
    
    % Calculate quantities for Time Of Emergence:
    if (PI_or_his)
    PIstd.stdTzp = std(Tzp,[],2);
    PIstd.stdTTp = std(TTp,[],2);
    PIstd.stdTyp = std(Typ,[],2);
    PIstd.stdOHC = std(OHC,[],1);
    save('PIstd.mat','PIstd');
    else
        load('PIstd.mat');
    end
    
    
% $$$     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$     %%% PLOTTING
% $$$     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$     
% $$$     %%%% Plot total heat content and fits time series:
% $$$     figure;
% $$$     set(gcf,'Position',[1 41 2560 1327.3]);
% $$$     subplot(5,1,1);
% $$$     plot(time,filter_field(N34,5,'-t'));
% $$$     hold on;
% $$$     plot(time,GMSST*10,'-r');
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
% $$$     plot(time,HTp(end,:),'--g');
% $$$     plot(time,Hzp(end,:),'--c');
% $$$     plot(time,Hyp(end,:),'--y');
% $$$     legend('OHC dedrifted',...
% $$$            'OHC from cumsum(TTp)','OHC from cumsum(Tzp)','OHC from cumsum(Typ)');
% $$$     xlabel('Year');
% $$$     ylabel('Ocean Heat Content (J)');
% $$$     xlim([yr1 yr1+nyrs]);
% $$$ 
% $$$     subplot(5,1,4);
% $$$     plot(time,AMOCfull,'-k');
% $$$     hold on;
% $$$     plot(time,filter_field(AMOCfull,5*12+1,'-t'),'-k','linewidth',2);
% $$$     xlabel('Year');
% $$$     ylabel('AMOC at $26^\circ$N (Sv)');
% $$$     xlim([yr1 yr1+nyrs]);
% $$$ 
% $$$     subplot(5,1,5);
% $$$     plot(time,WPOW,'-k');
% $$$     hold on;
% $$$     plot(time,filter_field(WPOW,5*12+1,'-t'),'-k','linewidth',2);
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
% $$$     plot(PIstd.stdTzp*2,P,'-k');
% $$$     hold on;
% $$$     plot(PIstd.stdTTp*2,P,'-r');
% $$$     plot(PIstd.stdTyp*2,P,'-b');
% $$$     xlabel('2$\sigma$ variability amplitude ($^\circ$C)');
% $$$     ylabel('Percentile');
% $$$     legend('$\Theta(p_z)$','$\Theta(p_\Theta)$',['$\Theta(p_\' ...
% $$$                         'phi)$']);
% $$$     
% $$$     % Example time series:
% $$$     pii = 20;
% $$$     [tmp piii] = min(abs(P-pii));
% $$$     figure;
% $$$     plot(time,Tzp(piii,:),'-k','linewidth',2);
% $$$     hold on;
% $$$     plot(time,TTp(piii,:),'-r','linewidth',2);
% $$$ % $$$     plot(time,Typ(piii,:),'-b','linewidth',2);
% $$$     plot([time(1) time(end)],[1 1]*2*PIstd.stdTzp(piii),'--k');
% $$$     plot([time(1) time(end)],[1 1]*-2*PIstd.stdTzp(piii),'--k');
% $$$     plot([time(1) time(end)],[1 1]*2*PIstd.stdTTp(piii),'--r');
% $$$     plot([time(1) time(end)],[1 1]*-2*PIstd.stdTTp(piii),'--r');
% $$$ % $$$     plot([time(1) time(end)],[1 1]*2*PIstd.stdTyp(piii),'--b');
% $$$ % $$$     plot([time(1) time(end)],[1 1]*-2*PIstd.stdTyp(piii),'--b');

    % Choose whether to remap to T for y-axis:
    remap_to_T = 0;
    if (remap_to_T)
        
        Z_Yvar = Tzp_mean;
        zlim = [-2 34];
        zlab = 'Horizontally-averaged Temperature ($^\circ$C)';

        T_Yvar = TTp_mean;
        tlim = [-2 34];
        tlab = 'Percentile-averaged Temperature ($^\circ$C)';
    
        y_Yvar = Typ_mean;
        yylim = [-80 80];
        ylab = 'Vertical/zonal-averaged Temperature ($^\circ$C)';
    else
        Z_Yvar = P;
        zlim = [0 100];
        zlab = 'Depth Percentile $p_z$';

        T_Yvar = P;
        tlim = zlim;
        tlab = 'Temperature Percentile $p_\Theta$';
        
        y_Yvar = P;
        yylim = zlim;
        ylab = 'Latitude Percentile $p_\phi$';
    end        
% $$$     
% $$$     %%%%% Plot mean, climatology, std and trends:
% $$$     figure;
% $$$     set(gcf,'Position',[1 41 2560 1330]);
% $$$ 
% $$$     subplot(3,4,1);
% $$$     plot(Tzp_mean,P,'linewidth',2);
% $$$     xlabel('Temperature ($^\circ$C)');
% $$$     xlim([-2 34]);
% $$$     ylabel('Depth Percentile');
% $$$     ylim([0 100]);
% $$$     title('ACCESS-CM2 PI-control mean $\Theta(z)$');
% $$$     set(gca,'ydir','reverse')
% $$$ 
% $$$     subplot(3,4,2);
% $$$     [X,Y] = ndgrid(1:12,Z_Yvar);
% $$$     contourf(X,Y,Tzp_clim',[-10 -0.2:0.002:0.2 10],'linestyle','none');
% $$$     xlabel('Month');
% $$$     ylabel(zlab);
% $$$     ylim(zlim);
% $$$     title('ACCESS-CM2 PI-control $\Theta(p_z)$ seasonal cycle');
% $$$     caxis([-0.2 0.2]);
% $$$     colorbar;
% $$$     if (~remap_to_T);set(gca,'ydir','reverse');end;
% $$$ 
% $$$     subplot(3,4,3);
% $$$     plot(std(Tzp_clim,[],2),Z_Yvar,'linewidth',2);
% $$$     hold on;
% $$$     plot(std(Tzp,[],2),Z_Yvar,'-r','linewidth',2);
% $$$     xlabel('std(Temperature) ($^\circ$C)');
% $$$     legend('Seasonal Climatology','Anomalies');
% $$$     ylabel(zlab);
% $$$     ylim(zlim);
% $$$     xlim([0 0.3]);
% $$$     title('ACCESS-CM2 PI-control variability $\Theta(p_z)$');
% $$$     if (~remap_to_T);set(gca,'ydir','reverse');end;
% $$$ 
% $$$     subplot(3,4,4);
% $$$     plot(Tzp_cubtr(2,:)',P,'linewidth',2);
% $$$     xlabel('Temperature trend ($^\circ$C/year)');
% $$$     ylabel('Depth Percentile');
% $$$     ylim(zlim);
% $$$     xlim([-1.5e-3 1.5e-3]);
% $$$     title('ACCESS-CM2 PI-control linear trend $\Theta(p_z)$');
% $$$     if (~remap_to_T);set(gca,'ydir','reverse');end;
% $$$ 
% $$$     subplot(3,4,5);
% $$$     plot(TTp_mean,P,'linewidth',2);
% $$$     xlabel('Temperature ($^\circ$C)');
% $$$     xlim([-2 34]);
% $$$     ylim([0 100]);
% $$$     ylabel('Temperature Percentile');
% $$$     title('ACCESS-CM2 PI-control mean $\Theta(p)$');
% $$$     set(gca,'ydir','reverse')
% $$$ 
% $$$     subplot(3,4,6);
% $$$     [X,Y] = ndgrid(1:12,T_Yvar);
% $$$     contourf(X,Y,TTp_clim',[-10 -0.2:0.002:0.2 10],'linestyle','none');
% $$$     xlabel('Month');
% $$$     ylabel(tlab);
% $$$     title('ACCESS-CM2 PI-control $\Theta(p_\Theta)$ seasonal cycle');
% $$$     colorbar;
% $$$     colormap(redblue);
% $$$     caxis([-0.2 0.2]);
% $$$     ylim(tlim);
% $$$     if (~remap_to_T);set(gca,'ydir','reverse');end;
% $$$ 
% $$$     subplot(3,4,7);
% $$$     plot(std(TTp_clim,[],2),T_Yvar,'linewidth',2);
% $$$     hold on;
% $$$     plot(std(TTp,[],2),T_Yvar,'-r','linewidth',2);
% $$$     xlabel('std(Temperature) ($^\circ$C)');
% $$$     legend('Seasonal Climatology','Anomalies');
% $$$     ylim(yylim);
% $$$     xlim([0 0.3]);
% $$$     ylabel(ylab);
% $$$     title('ACCESS-CM2 PI-control variability $\Theta(p_\Theta)$');
% $$$     if (~remap_to_T);set(gca,'ydir','reverse');end;
% $$$ % $$$     set(gca,'Position',[0.1300    0.5838    0.2504    0.3412]);
% $$$ 
% $$$     subplot(3,4,8);
% $$$     plot(TTp_cubtr(2,:)',T_Yvar,'linewidth',2);
% $$$     xlabel('Temperature trend ($^\circ$C/year)');
% $$$     xlim([-1.5e-3 1.5e-3]);
% $$$     ylim(yylim);
% $$$     ylabel(ylab);
% $$$     title('ACCESS-CM2 PI-control linear trend $\Theta(p_\Theta)$');
% $$$     if (~remap_to_T);set(gca,'ydir','reverse');end;
% $$$ 
% $$$     subplot(3,4,9);
% $$$     plot(Typ_mean,P,'linewidth',2);
% $$$     xlabel('Temperature ($^\circ$C)');
% $$$     xlim([-1 6.5]);
% $$$     ylim([0 100]);
% $$$     ylabel('Latitude Percentile');
% $$$     title('ACCESS-CM2 PI-control mean $\phi(p)$');
% $$$ 
% $$$     subplot(3,4,10);
% $$$     [X,Y] = ndgrid(1:12,y_Yvar);
% $$$     contourf(X,Y,Typ_clim',[-10 -0.2:0.002:0.2 10],'linestyle','none');
% $$$     xlabel('Month');
% $$$     ylabel(ylab);
% $$$     title('ACCESS-CM2 PI-control $\Theta(p_\phi)$ seasonal cycle');
% $$$     colorbar;
% $$$     colormap(redblue);
% $$$     caxis([-0.2 0.2]);
% $$$     ylim(yylim);
% $$$ 
% $$$     subplot(3,4,11);
% $$$     plot(std(Typ_clim,[],2),y_Yvar,'linewidth',2);
% $$$     hold on;
% $$$     plot(std(Typ,[],2),y_Yvar,'-r','linewidth',2);
% $$$     xlabel('std(Temperature) ($^\circ$C)');
% $$$     legend('Seasonal Climatology','Anomalies');
% $$$     ylim(yylim);
% $$$     xlim([0 0.3]);
% $$$     ylabel(ylab);
% $$$     title('ACCESS-CM2 PI-control variability $\Theta(p_\phi)$');
% $$$ % $$$     set(gca,'Position',[0.1300    0.5838    0.2504    0.3412]);
% $$$ 
% $$$     subplot(3,4,12);
% $$$     plot(Typ_cubtr(2,:)',y_Yvar,'linewidth',2);
% $$$     xlabel('Temperature trend ($^\circ$C/year)');
% $$$     xlim([-1.5e-3 1.5e-3]);
% $$$     ylim(yylim);
% $$$     ylabel(ylab);
% $$$     title('ACCESS-CM2 PI-control linear trend $\Theta(p_\phi)$');
% $$$ 
% $$$ 
% $$$     %%% Straight comparison of amplitude of anomalies:
% $$$     figure;
% $$$     set(gcf,'Position',[1 41 2560 1330]);
% $$$     subplot(1,2,1);
% $$$     plot(std(Tzp,[],2),P,'-k','linewidth',2);
% $$$     hold on;
% $$$     plot(std(TTp,[],2),P,'-r','linewidth',2);
% $$$     plot(std(Typ,[],2),P,'-b','linewidth',2);
% $$$     xlabel('Standard Deviation of $\Theta$ ($^\circ$C)');
% $$$     legend('$\Theta(p_z)$','$\Theta(p_\Theta)$','$\Theta(p_\phi)$');
% $$$     ylabel('Percentile');
% $$$     ylim([0 100]);
% $$$     xlim([0 0.3]);
% $$$     title('ACCESS-CM2 PI-control variability');
% $$$     set(gca,'ydir','reverse');
% $$$     
% $$$     % Averages:
% $$$     text(0.1,30,['$\Theta(p_z)$ average = ' sprintf('%3.4f',mean(std(Tzp,[],2))) '$^\circ$C'],'color','k');
% $$$     text(0.1,35,['$\Theta(p_\Theta)$ average = ' sprintf('%3.4f',mean(std(TTp,[],2))) '$^\circ$C'],'color','r');
% $$$     text(0.1,40,['$\Theta(p_\phi)$ average = ' sprintf('%3.4f',mean(std(Typ,[],2))) '$^\circ$C'],'color','b');
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
% $$$         Pticks(ii) = interp1(TTp_mean,P,Zticks(ii),'linear');
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
% $$$     plot(std(Hzp,[],2),P,'-k','linewidth',2);
% $$$     hold on;
% $$$     plot(std(HTp,[],2),P,'-r','linewidth',2);
% $$$     plot(std(Hyp,[],2),P,'-b','linewidth',2);
% $$$     xlabel('Standard Deviation of heat content (J)');
% $$$     legend('$H(p_z)$','$H(p_\Theta)$','$H(p_\phi)$');
% $$$     ylabel('Percentile');
% $$$     ylim([0 100]);
% $$$     title('ACCESS-CM2 PI-control variability');
% $$$     set(gca,'ydir','reverse');
        
    % Choose whether to plot T or H:
    plot_H = 0;
    if (plot_H)
        Zvar = Hzp;
        Zcxs = [-0.5e23 1e21 0.5e23];
        Zlab = '$H(p_z)$';

        Tvar = HTp;
        Tcxs = [-0.5e23 1e21 0.5e23];
        Tlab = '$H(p_\Theta)$';
    
        yvar = Hyp;
        ycxs = [-0.5e23 1e21 0.5e23];
        Ylab = '$H(p_\phi)$';
else
        Zvar = Tzp;
        Zcxs = [-0.2 0.02 0.2];
        Zlab = '$\Theta(p_z)$';
        Tvar = TTp;
        Tcxs = Zcxs;%[-0.1 0.01 0.1];
        Tlab = '$\Theta(p_\Theta)$';
        yvar = Typ;
        ycxs = Zcxs;%[-0.1 0.01 0.1];
        Ylab = '$\Theta(p_\phi)$';
    end        
        
    %%%%%% Plot all anomalies time series:

    xlims = [yr1 yr1+nyrs];

% $$$     % ENSO Focus settings:
% $$$     zlim = [0 10];
% $$$     tlim = [0 10];
% $$$     yylim = [40 85];
% $$$     xlims = [300 350];
% $$$     Zcxs = [-0.2 0.01 0.2];
% $$$     Tcxs = [-0.2 0.01 0.2];
% $$$     ycxs = [-0.2 0.01 0.2];
% $$$ 
    [tmp t1] = min(abs(time-xlims(1)));
    [tmp t2] = min(abs(time-xlims(2)));

    figure;
    set(gcf,'Position',[1 41 2560 1330]);
    set(gcf,'defaulttextfontsize',15);
    set(gcf,'defaultaxesfontsize',15);
    [X,Y] = ndgrid(time(t1:t2),Z_Yvar);
    subplot(3,1,1);
    contourf(X,Y,Zvar(:,t1:t2)',[-1e50 Zcxs(1):Zcxs(2):Zcxs(3) 1e50],'linestyle','none');
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
    
    [X,Y] = ndgrid(time(t1:t2),T_Yvar);
    subplot(3,1,2);
    contourf(X,Y,Tvar(:,t1:t2)',[-1e50 Tcxs(1):Tcxs(2):Tcxs(3) 1e50],'linestyle','none');
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
            Pticks(ii) = interp1(TTp_mean,P,Zticks(ii),'linear');
        end
        set(ax2,'ytick',Pticks);
        set(ax2,'yticklabel',Zticks);
        ylabel(ax2,'Temperature ($^\circ$C)');
        ylim(ax2,zlim);
        set(ax2,'ydir','reverse');
    end

    [X,Y] = ndgrid(time(t1:t2),y_Yvar);
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
% $$$     [pxx,f] = pmtm(Zvar(15,:)',[],[],fs);
% $$$     fL = length(f);
% $$$     Tvar_spec = zeros(fL,PL);
% $$$     Zvar_spec = zeros(fL,PL);
% $$$     yvar_spec = zeros(fL,PL);
% $$$     for pi=1:PL
% $$$         [Zvar_spec(:,pi),~] = pmtm(Zvar(pi,:)',[],[],fs);
% $$$         [Tvar_spec(:,pi),~] = pmtm(Tvar(pi,:)',[],[],fs);
% $$$         [yvar_spec(:,pi),~] = pmtm(yvar(pi,:)',[],[],fs);
% $$$     end
% $$$     
% $$$     figure;
% $$$     set(gcf,'Position',[1 41 2560 1330]);
% $$$     subplot(3,1,1);
% $$$     [X,Y] = ndgrid(f,Z_Yvar);
% $$$     contourf(X,Y,log10(Zvar_spec),[-4.5:0.25:-0.25],'linestyle','none');
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
% $$$     [X,Y] = ndgrid(f,T_Yvar);
% $$$     contourf(X,Y,log10(Tvar_spec),[-4.5:0.25:-0.25],'linestyle','none');
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
% $$$             Pticks(ii) = interp1(TTp_mean,P,Zticks(ii),'linear');
% $$$         end
% $$$         set(ax2,'ytick',Pticks);
% $$$         set(ax2,'yticklabel',Zticks);
% $$$         ylabel(ax2,'Temperature ($^\circ$C)');
% $$$         ylim(ax2,zlim);
% $$$         set(ax2,'ydir','reverse');
% $$$     end
% $$$ 
% $$$     subplot(3,1,3);
% $$$     [X,Y] = ndgrid(f,y_Yvar);
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
% $$$     plot(f,log10(mean(Zvar_spec,2)),'-k','linewidth',2);
% $$$     hold on;
% $$$     plot(f,log10(mean(Tvar_spec,2)),'-r','linewidth',2);
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
        EOFs = Tzp*PC/tL;
        
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
    [X,Y] = ndgrid(time(t1:t2),Z_Yvar);
    subplot(3,1,1);
    contourf(X,Y,Zvar(:,t1:t2)',[-1e50 Zcxs(1):Zcxs(2):Zcxs(3) 1e50],'linestyle','none');
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
% $$$     Zvar_lr = zeros(PL,ll);
% $$$     Tvar_lr = zeros(PL,ll);
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
% $$$         Zvar_lr(:,ii) = Zvar*TS/(sum(TS.^2));
% $$$         Tvar_lr(:,ii) = Tvar*TS/(sum(TS.^2));
% $$$         yvar_lr(:,ii) = yvar*TS/(sum(TS.^2));
% $$$     end
% $$$ 
% $$$     figure;
% $$$     set(gcf,'Position',[1 41 2560 1330]);
% $$$ 
% $$$     subplot(3,1,1);
% $$$     [X,Y] = ndgrid(lags/12,Z_Yvar);
% $$$     contourf(X,Y,Zvar_lr',[-1e50 Zcxs(1):Zcxs(2):Zcxs(3) 1e50],'linestyle','none');
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
% $$$     [X,Y] = ndgrid(lags/12,T_Yvar);
% $$$     contourf(X,Y,Tvar_lr',[-1e50 Tcxs(1):Tcxs(2):Tcxs(3) 1e50],'linestyle','none');
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
% $$$             Pticks(ii) = interp1(TTp_mean,P,Zticks(ii),'linear');
% $$$         end
% $$$         set(ax2,'ytick',Pticks);
% $$$         set(ax2,'yticklabel',Zticks);
% $$$         ylabel(ax2,'Temperature ($^\circ$C)');
% $$$         ylim(ax2,zlim);
% $$$         set(ax2,'ydir','reverse');
% $$$     end
% $$$ 
% $$$     subplot(3,1,3);
% $$$     [X,Y] = ndgrid(lags/12,y_Yvar);
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
% $$$     plot(TTap_mean,P,'linewidth',2);
% $$$     xlabel('Temperature ($^\circ$C)');
% $$$     xlim([-2 34]);
% $$$     ylabel('Surface Area Percentile');
% $$$     ylim([0 100]);
% $$$     title('ACCESS-CM2 PI-control mean $\Theta(p_\Theta)$');
% $$$     set(gca,'ydir','reverse')
% $$$ 
% $$$     subplot(2,3,2);
% $$$     [X,Y] = ndgrid(1:12,T_Yvar);
% $$$     contourf(X,Y,TTap_clim',[-10 -1:0.01:1 10],'linestyle','none');
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
% $$$     plot(std(TTap_clim,[],2),T_Yvar,'linewidth',2);
% $$$     hold on;
% $$$     plot(std(TTap,[],2),T_Yvar,'-r','linewidth',2);
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
% $$$     contourf(X,Y,TTap(:,t1:t2)',[-1e50 -1:0.02:1 1e50],'linestyle','none');
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
% $$$             Pticks(ii) = interp1(TTp_mean,P,Zticks(ii),'linear');
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
% $$$     [pxx,f] = pmtm(Zvar(15,:)',[],[],fs);
% $$$     fL = length(f);
% $$$     Tvar_spec = zeros(fL,PL);
% $$$     for pi=1:PL
% $$$         [Tvar_spec(:,pi),~] = pmtm(TTap(pi,:)',[],[],fs);
% $$$     end
% $$$     
% $$$     figure;
% $$$     set(gcf,'Position',[1 41 2560 1330]);
% $$$     subplot(2,1,1);
% $$$     [X,Y] = ndgrid(f,T_Yvar);
% $$$     contourf(X,Y,log10(Tvar_spec),[-4.5:0.25:-0.25],'linestyle','none');
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
% $$$             Pticks(ii) = interp1(TTp_mean,P,Zticks(ii),'linear');
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
