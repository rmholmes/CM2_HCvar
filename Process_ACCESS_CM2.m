% Process heat budget terms and outputs for ACCESS-CM2 run
%
% INPUTS:
% fname = six-month ocean_month.nc file name to process.
% oname = output .mat file-name.
% msk   = mask key

% Choices:
dT = 0.2; % temperature grid size.

%%%% Grid (time-constant) and time info:

% Constants:
Cp = 3992.10322329649; % J kg-1 degC-1
rho0 = 1035; % kgm-3

% Define a temperature grid:
Te = -3:dT:34;
T = (Te(2:end)+Te(1:end-1))/2;
TL = length(T);
    
% Constant grid parameters:
area = ncread(fname,'area_t');
lon = ncread(fname,'geolon_t');
lat = ncread(fname,'geolat_t');
lonv = ncread(fname,'xt_ocean');
latv = ncread(fname,'yt_ocean');
latu = ncread(fname,'yu_ocean');
latv_edges = cat(1,latv(1)-(latu(1)-latv(1)),latu);
[xL,yL] = size(area);    

% Depth grid:
Z = ncread(fname,'st_ocean');
Ze = ncread(fname,'st_edges_ocean');
zL = length(Z);

% 3D mask:
temp = ncread(fname,'temp',[1 1 1 1],[xL yL zL 1]);
mask = ~isnan(temp);

% region mask:
if (msk == '')
    'Global'
    mask2D = ones(xL,yL);
elseif (msk = 'NH')
    'Northern Hemisphere'
    mask2D = 1*(lat>0);
elseif (msk = 'SH')
    'Southern Hemisphere'
    mask2D = 1*(lat<=0);
end

% A(z):
A = zeros(zL,1);
for zi=1:zL
    A(zi) = nansum(nansum(area(mask(:,:,zi).*mask2D)));
end

% Time info:
time = ncread(fname,'time');
DT_A = ncread(fname,'average_DT')*86400;
tL = length(time);

% Save grid info (doesn't vary with time):
save(oname,'Cp','rho0','dT','Te','T','TL',...
     'xL','yL','zL','Z','Ze','zL','A', ...
     'latv','latv_edges','time','DT_A','tL');

%%%% Load processing variables:
temp = ncread(fname,'temp');
temp(~mask) = NaN;
if (max(max(temp))>120); temp=temp-273.15;end;
V = ncread(fname,'dht').*repmat(area.*mask2D,[1 1 zL tL]);
V(isnan(V)) = 0;
H = rho0*Cp*temp.*V;
SST = squeeze(temp(:,:,1,:));

%%%% Heat and volume:

% Latitude space:
Yv.H = squeeze(nansum(nansum(H,1),3));
Yv.V = squeeze(nansum(nansum(V,1),3));

% Depth space H and V:
Zv.H = squeeze(nansum(nansum(H,1),2));
Zv.V = squeeze(nansum(nansum(V,1),2));

% Temp space H and V:
Tv.V = zeros(TL,tL);
Tv.H = zeros(TL,tL);
Tv.A = zeros(TL,tL);
for Ti=1:TL
    %Accumulate sums:
    inds = temp>=Te(Ti) & temp<Te(Ti+1);
    Tv.V(Ti,:) = Tv.V(Ti,:) + squeeze(nansum(nansum(nansum(V.*inds,1),2),3))';
    Tv.H(Ti,:) = Tv.H(Ti,:) + squeeze(nansum(nansum(nansum(H.*inds,1),2),3))';
    indsS = SST>=Te(Ti) & SST<Te(Ti+1);
    Tv.A(Ti,:) = Tv.A(Ti,:) + squeeze(nansum(nansum(repmat(area.*mask2D,[1 1 tL]).*indsS,1),2))';
end
% Account for water warmer than max temperature (possible for CM2):
inds = temp>Te(end);
Tv.V(end,:) = Tv.V(end,:) + squeeze(nansum(nansum(nansum(V.*inds,1),2),3))';
Tv.H(end,:) = Tv.H(end,:) + squeeze(nansum(nansum(nansum(H.*inds,1),2),3))';
indsS = SST>Te(end);
Tv.A(end,:) = Tv.A(end,:) + squeeze(nansum(nansum(repmat(area.*mask2D,[1 1 tL]).*indsS,1),2))';

% $$$ Yv.time = time; Tv.DT_A = DT_A;
% $$$ Zv.time = time; Zv.DT_A = DT_A;
% $$$ Tv.time = time; Tv.DT_A = DT_A;
save(oname,'Yv','Zv','Tv','-append');

%%%% Climate indices:

% grid locations:
[tmp, N34x1] = min(abs(lonv+170));  [tmp, N34x2] = min(abs(lonv+120));
[tmp, N34y1] = min(abs(latv+5));    [tmp, N34y2] = min(abs(latv-5));

[tmp, TPIr1x1] = min(abs(lonv+220));  [tmp, TPIr1x2] = min(abs(lonv+145));
[tmp, TPIr1y1] = min(abs(latv-25));    [tmp, TPIr1y2] = min(abs(latv-45));
[tmp, TPIr2x1] = min(abs(lonv+190));  [tmp, TPIr2x2] = min(abs(lonv+90));
[tmp, TPIr2y1] = min(abs(latv+10));    [tmp, TPIr2y2] = min(abs(latv-10));
[tmp, TPIr3x1] = min(abs(lonv+230));  [tmp, TPIr3x2] = min(abs(lonv+160));
[tmp, TPIr3y1] = min(abs(latv+50));    [tmp, TPIr3y2] = min(abs(latv+15));

potrho = ncread(fname,'potrho_edges');
potrhoL = length(potrho)-1;

[tmp, AMOClt] = min(abs(latv-26)); [tmp, AMOCrho] = min(abs(potrho-1035.5));
[tmp, AMOCln1] = min(abs(lonv+103));    [tmp, AMOCln2] = min(abs(lonv+5));

% Calculate SST based indices:
CIN.N34 = squeeze(nanmean(nanmean(temp(N34x1:N34x2,N34y1:N34y2,1,:),1),2));

CIN.TPIr1 = squeeze(nanmean(nanmean(temp(TPIr1x1:TPIr1x2,TPIr1y1:TPIr1y2,1,:),1),2));
CIN.TPIr2 = squeeze(nanmean(nanmean(temp(TPIr2x1:TPIr2x2,TPIr2y1:TPIr2y2,1,:),1),2));
CIN.TPIr3 = squeeze(nanmean(nanmean(temp(TPIr3x1:TPIr3x2,TPIr3y1:TPIr3y2,1,:),1),2));

% AMOC index:
CIN.AMOC = (-1/1e6/rho0)*squeeze(nansum( ...
    nansum(ncread(fname,'ty_trans_rho',[AMOCln1 AMOClt AMOCrho 1],[AMOCln2-AMOCln1+1 1 potrhoL-AMOCrho+1 tL]),3) + ...
    ncread(fname,'ty_trans_rho_gm',[AMOCln1 AMOClt AMOCrho 1],[AMOCln2-AMOCln1+1 1 1 tL]),1));

% WPOW index:
CIN.WPOW = squeeze(nansum(nansum(ncread(fname,'wind_power_u')+ncread(fname,'wind_power_v'),1),2));

% $$$ CIN.time = time; CIN.DT_A = DT_A;
save(oname,'CIN','-append');

%%% BUDGET variables:

% Variable lists:
bvars3D = {'temp_tendency','temp_submeso', ...
           'temp_vdiffuse_diff_cbt', 'temp_nonlocal_KPP', ...
           'temp_vdiffuse_sbc','frazil_3d','sw_heat','temp_rivermix', ...
           'neutral_diffusion_temp','neutral_gm_temp','temp_vdiffuse_k33',...
           'mixdownslope_temp','temp_sigma_diff'};
bvars2D = {'sfc_hflux_pme','temp_eta_smooth','pme_river'}; % pme_river is mass_pmepr_on_nrho not binned (see ocean_sbc).
bvars = {bvars3D{:},bvars2D{:}};
% NOTE: temp_advection is corrupted in the CM2 output - reconstruct by residual.

% Initialize T vars:
for vi =1:length(bvars)
    eval(['Tv.' bvars{vi} ' = zeros(TL,tL);']);
end

% Loop through tendency variables:
for vi =1:length(bvars3D)
    var = ncread(fname,bvars3D{vi}).*repmat(area.*mask2D,[1 1 zL tL]);
    
    % Depth binning:
    eval(['Zv.' bvars3D{vi} ' = squeeze(nansum(nansum(var,1),2));']);
    
    % Latitude binning:
    eval(['Yv.' bvars3D{vi} ' = squeeze(nansum(nansum(var,1),3));']);

    % Temperature binning:
    for Ti=1:TL
        %Accumulate sums:
        inds = temp>=Te(Ti) & temp<Te(Ti+1);
        eval(['Tv.' bvars3D{vi} '(Ti,:) = Tv.' bvars3D{vi} ...
              '(Ti,:) + squeeze(nansum(nansum(nansum(var.*inds,1),2),3))'';']);
    end
    inds = temp>Te(end);
    eval(['Tv.' bvars3D{vi} '(end,:) = Tv.' bvars3D{vi} ...
          '(end,:) + squeeze(nansum(nansum(nansum(var.*inds,1),2),3))'';']);
end
for vi =1:length(bvars2D)
    var = ncread(fname,bvars2D{vi}).*repmat(area.*mask2D,[1 1 tL]);

    % Depth binning:
    varz = squeeze(nansum(nansum(var,1),2));
    eval(['Zv.' bvars2D{vi} ' = cat(1,varz'',zeros(zL-1,tL));']);
    
    % Latitude binning:
    eval(['Yv.' bvars2D{vi} ' = squeeze(nansum(var,1));']);
    
    % Temperature binning:
    for Ti=1:TL
        indsS = SST>=Te(Ti) & SST<Te(Ti+1);
        eval(['Tv.' bvars2D{vi} '(Ti,:) = Tv.' bvars2D{vi} ...
              '(Ti,:) + squeeze(nansum(nansum(var.*indsS,1),2))'';']);
    end
    indsS = SST>Te(end);
    eval(['Tv.' bvars2D{vi} '(end,:) = Tv.' bvars2D{vi} ...
          '(end,:) + squeeze(nansum(nansum(var.*indsS,1),2))'';']);
end

save(oname,'Tv','Zv','Yv','-append');


