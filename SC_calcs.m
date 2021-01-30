%% Calculations of F, N, SST':

fname = '/srv/ccrc/data03/z3500785/access-cm2/ocean_month_SC.nc';

area = ncread(fname,'area_t');
lat = ncread(fname,'geolat_t');
[xL,yL] = size(area);

DT = ncread(fname,'time_bounds');
DT = diff(DT,[],1);
DTc = cat(2,0,cumsum(DT));
time = ncread(fname,'time');

% Surface flux:
Qnet = squeeze(nansum(ncread(fname,'temp_vdiffuse_sbc'),3))+ncread(fname,'sfc_hflux_pme')+...
       squeeze(nansum(ncread(fname,'temp_rivermix'),3))+squeeze(nansum(ncread(fname,'frazil_3d'),3));
Qnet = Qnet.*repmat(area,[1 1 12]); % multiply by area
Qnet = Qnet - repmat(monmean(Qnet,3,DT),[1 1 12]);

% SST anomalies:
SST = ncread(fname,'temp');
SST = squeeze(SST(:,:,1,:));
SST = SST - repmat(monmean(SST,3,DT),[1 1 12]);

% Solar irradiance:
days = 0.5:1:(sum(DT)-0.5);
Qsr = zeros(xL,yL,12);
for xi = 1:xL
    [PHI,DAYS] = ndgrid(lat(xi,:),days);
    Qsi = Qday(PHI,DAYS);
    for ti = 1:12
        Qsr(xi,:,ti) = mean(Qsi(:,(DTc(ti)+1):DTc(ti+1)),2);
    end
end 
Qsr = Qsr.*repmat(area,[1 1 12]); % multiply by area
Qsr = Qsr - repmat(monmean(Qsr,3,DT),[1 1 12]);

% Masking/NH/SH
mask = ~isnan(SST(:,:,1));
maskNH = mask & lat>0;
maskSH = mask & lat<0;

areaNH = sum(area(maskNH));
areaSH = sum(area(maskSH));
QnetNH = zeros(12,1);
QnetSH = zeros(12,1);
QsrNH = zeros(12,1);
QsrSH = zeros(12,1);

SSTNH = zeros(12,1);
SSTSH = zeros(12,1);

for ti = 1:12
    Qneti = Qnet(:,:,ti);
    SSTi = SST(:,:,ti);
    Qsri = Qsr(:,:,ti);

    QnetNH(ti) = sum(Qneti(maskNH));
    QnetSH(ti) = sum(Qneti(maskSH));
    QsrNH(ti) = sum(Qsri(maskNH));
    QsrSH(ti) = sum(Qsri(maskSH));

    SSTNH(ti) = sum(SSTi(maskNH).*area(maskNH))/areaNH;
    SSTSH(ti) = sum(SSTi(maskSH).*area(maskSH))/areaSH;
end

figure;
subplot(2,1,1);
plot(1:12,QnetNH/areaNH,'-k');
hold on;
plot(1:12,QnetSH/areaSH,'-r');
plot(1:12,QsrNH/areaNH,'--k');
plot(1:12,QsrSH/areaSH,'--r');
xlim([1 12]);
ylabel('N'' and F'' (Wm$^{-2}$)');
legend('NH N''','SH N''','NH F''','SH F''');
subplot(2,1,2);
plot(1:12,SSTNH,'-k');
hold on;
plot(1:12,SSTSH,'-r');
xlim([1 12]);
ylabel('SST anomaly ($^\circ$C)');
xlabel('Month');
legend('NH','SH');

% response parameters:
alphaNH = ((QsrNH - QnetNH)./SSTNH)/areaNH;
alphaSH = ((QsrSH - QnetSH)./SSTSH)/areaSH;

figure;
plot(1:12,alphaNH);
hold on;
plot(1:12,alphaSH);

% Just straight amplitudes:
QnetNHamp = sqrt(var(QnetNH))/areaNH;
QnetSHamp = sqrt(var(QnetSH))/areaSH;
QsrNHamp = sqrt(var(QsrNH))/areaNH;
QsrSHamp = sqrt(var(QsrSH))/areaSH;
SSTNHamp = sqrt(var(SSTNH));
SSTSHamp = sqrt(var(SSTSH));

alphaNH = ((QsrNHamp-QnetNHamp)/SSTNHamp);
alphaSH = ((QsrSHamp-QnetSHamp)/SSTSHamp);

figure;
subplot(1,4,2);
bar([1 2],[QnetNHamp QnetSHamp]);
ylabel('N'' (Wm$^{-2}$)');
set(gca,'xtick',[1 2]);
set(gca,'xticklabel',{'NH','SH'});
ylim([0 120]);
xlim([0 3]);
subplot(1,4,1);
bar([1 2],[QsrNHamp QsrSHamp]);
ylabel('F'' (Wm$^{-2}$)');
set(gca,'xtick',[1 2]);
set(gca,'xticklabel',{'NH','SH'});
ylim([0 120]);
xlim([0 3]);
subplot(1,4,3);
bar([1 2],[SSTNHamp SSTSHamp]);
ylabel('SST'' ($^\circ$C)');
set(gca,'xtick',[1 2]);
set(gca,'xticklabel',{'NH','SH'});
ylim([0 2]);
xlim([0 3]);
subplot(1,4,4);
bar([1 2],[alphaNH alphaSH]);
ylabel('$\alpha$ (Wm$^{-2}$degC$^{-1}$)');
set(gca,'xtick',[1 2]);
set(gca,'xticklabel',{'NH','SH'});
ylim([0 25]);
xlim([0 3]);

%% T-z overturning using tendencies:

fname = '/srv/ccrc/data03/z3500785/access-cm2/ocean_month_SC.nc';
area = ncread(fname,'area_t');
[xL,yL] = size(area);
DT = ncread(fname,'time_bounds');
DT = diff(DT,[],1);

% Constants:
Cp = 3992.10322329649; % J kg-1 degC-1
rho0 = 1035; % kgm-3

% Define a temperature grid:
dT = 1;
Te = -3:dT:34;
T = (Te(2:end)+Te(1:end-1))/2;
TL = length(T);

% Depth grid:
Z = ncread(fname,'st_ocean');
Ze = ncread(fname,'st_edges_ocean');
zL = length(Z);

bvars3D = {'temp_submeso', ...
           'temp_vdiffuse_diff_cbt', 'temp_nonlocal_KPP', ...
           'temp_vdiffuse_sbc','frazil_3d','sw_heat','temp_rivermix', ...
           'neutral_diffusion_temp','neutral_gm_temp','temp_vdiffuse_k33',...
           'mixdownslope_temp','temp_sigma_diff'};
bvars2D = {'sfc_hflux_pme','temp_eta_smooth'};
bvars = {bvars3D{:},bvars2D{:}};

adv = 0;

if (adv)

temp_advection = ncread(fname,'temp_tendency');
for vi = 1:length(bvars3D)
    var = ncread(fname,bvars3D{vi});
    temp_advection = temp_advection-var;
end
for vi = 1:length(bvars2D)
    var = ncread(fname,bvars2D{vi});
    temp_advection(:,:,1,:) = temp_advection(:,:,1,:)-permute(var,[1 ...
                       2 4 3]);
end

tend = temp_advection + ncread(fname,'neutral_gm_temp')+ncread(fname,'temp_submeso');
else
tend = ncread(fname,'temp_vdiffuse_diff_cbt') + ...
       ncread(fname,'temp_nonlocal_KPP') + ...
       ncread(fname,'temp_vdiffuse_sbc') + ...
       ncread(fname,'frazil_3d') + ...
       ncread(fname,'sw_heat') + ...
       ncread(fname,'temp_rivermix') + ...
       ncread(fname,'neutral_diffusion_temp') + ...
       ncread(fname,'temp_vdiffuse_k33') + ...
       ncread(fname,'mixdownslope_temp') + ...
       ncread(fname,'temp_sigma_diff');% + ...
tend(:,:,1,:) = tend(:,:,1,:) + permute(ncread(fname,'sfc_hflux_pme'),[1 2 4 3]);
end

temp = ncread(fname,'temp')-273.15;
% $$$ V = ncread(fname,'dht').*repmat(area,[1 1 zL 12]);


% Temp space calculations:
PSI = zeros(TL,zL,12);
for ti=1:12
    ti
    for Ti=1:TL
        for zi=1:zL
            tempn = temp(:,:,zi,ti);
            tendn = tend(:,:,zi,ti);
% $$$             Vn = V(:,:,zi,ti);
            inds = tempn>=Te(Ti) & tempn<Te(Ti+1);
            
            
% $$$             PSI(Ti,zi,ti) = nansum(nansum(tendn(inds).*Vn(inds),1),2);
            PSI(Ti,zi,ti) = nansum(nansum(tendn(inds).*area(inds),1),2)/rho0/Cp/dT;
        end
    end
end

% Cumsum in z:
PSIc = cat(2,zeros(TL,1,12),cumsum(PSI,2,'reverse'));

figure;
set(gcf,'Position',[1         316        1280         690]);
[X,Y] = ndgrid(T,-Ze);
contourf(X,Y,monmean(PSIc,3,DT)/1e6,[-1e10 -50:1:50 1e10], ...
         'linestyle','none');
caxis([-25 25]);
cb = colorbar;
ylabel(cb,'Sv');
colormap('redblue');
ylim([-3000 0]);
xlabel('Temperature ($^\circ$C)');
ylabel('Depth (m)');
if (adv)
    title('T-z overturning from advection tendency');
else
    title('T-z overturning from diabatic tendencies');
end

figure;
set(gcf,'Position',[1         316        1280         690]);
set(gcf,'defaulttextfontsize',10);
set(gcf,'defaultaxesfontsize',10);
labs = {'Jan','Feb','Mar','Apr','May','Jun',...
        'Jul','Aug','Sep','Oct','Nov','Dec'};
for i=1:12
    subplot(3,4,i);
    contourf(X,Y,PSIc(:,:,i)/1e6,[-1e10 -25:1:25 1e10], ...
             'linestyle','none');
    caxis([-25 25]);
    colormap('redblue');
    ylim([-3000 0]);
    if (i == 1 | i == 5 | i == 9)
        ylabel('Depth (m)');
    else
        set(gca,'yticklabel',[]);
    end
    if (i >= 9)
        xlabel('Temperature ($^\circ$C)');
    else
        set(gca,'xticklabel',[]);
    end    
    title(labs{i});
end

    

