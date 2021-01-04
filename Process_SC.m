% Calculate a seasonal cycle from the CM2 output.
%

base = '/g/data/p66/cm2704/archive/bi889/history/ocn/';
name = 'PIcontrol'
obase = '/scratch/e14/rmh561/access-cm2/HCvar/'

Tvars3D = {'temp_tendency','temp_submeso', ...
           'temp_vdiffuse_diff_cbt', 'temp_nonlocal_KPP', ...
           'temp_vdiffuse_sbc','frazil_3d','sw_heat','temp_rivermix', ...
           'neutral_diffusion_temp','neutral_gm_temp','temp_vdiffuse_k33',...
           'mixdownslope_temp','temp_sigma_diff'};
Tvars2D = {'sfc_hflux_pme','temp_eta_smooth','pme_river'}; % pme_river is mass_pmepr_on_nrho not binned (see ocean_sbc).
Tvars = {Tvars3D{:},Tvars2D{:}};

Svars3D = {'salt_tendency','salt_advection','salt_submeso', ...
           'salt_vdiffuse_diff_cbt', 'salt_nonlocal_KPP', ...
           'salt_vdiffuse_sbc','salt_rivermix', ...
           'neutral_diffusion_salt','neutral_gm_salt','salt_vdiffuse_k33',...
           'mixdownslope_salt','salt_sigma_diff'};
% $$$ sfc_salt_flux_coupler
% $$$ sfc_salt_flux_pme
% $$$ sfc_salt_flux_prec
% $$$ sfc_salt_flux_evap
Svars2D = {'salt_eta_smooth'}; % pme_river is mass_pmepr_on_nrho not binned (see ocean_sbc).
Svars = {Svars3D{:},Svars2D{:}};

Nvars = {'temp','salt','dht','area_t','geolon_t','geolat_t','st_ocean','st_edges_ocean'};
vars = {Nvars{:},Tvars{:},Svars{:}};




fname = [base 'ocean_month.nc-

         
% Constants:
Cp = 3992.10322329649; % J kg-1 degC-1
rho0 = 1035; % kgm-3

% Constant grid parameters:
area = ncread(fname,'area_t');
lon = ncread(fname,'geolon_t');
lat = ncread(fname,'geolat_t');
Z = ncread(fname,'st_ocean');
Ze = ncread(fname,'st_edges_ocean');

[xL,yL] = size(area);    
zL = length(Z);

% 3D mask:
temp = ncread(fname,'temp',[1 1 1 1],[xL yL zL 1]);
mask = ~isnan(temp);

temp = zeros(xL,yL,zL,12);
salt = zeros(xL,yL,zL,12);


% $$$ % Save grid info (doesn't vary with time):
% $$$ save(oname,'Cp','rho0','dT','Te','T','TL',...
% $$$      'xL','yL','zL','Z','Ze','zL','A', ...
% $$$      'latv','latv_edges','time','DT_A','tL');

temp = ncread(fname,'temp'); % Conservative Temperatuere
salt = ncread(fname,'salt'); % Practical Salinity
dzt = ncread(fname,'dzt');   % Time-varying Layer thickness

% Variable lists:
% NOTE: temp_advection is corrupted in the CM2 output - reconstruct by residual.


