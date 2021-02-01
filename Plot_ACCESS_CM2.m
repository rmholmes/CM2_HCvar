%%%% Plotting scripts for HC variability project

% Load specific experiments:
clear all;
baseMAT = '/srv/ccrc/data03/z3500785/access-cm2/';
% $$$ load([baseMAT 'PIcontrolPP.mat']);
% $$$ load([baseMAT 'PIcontrolPP_hpi.mat']);
load([baseMAT 'PIcontrolPP_Tb05.mat']);

%%% Percentile converter:
Ps = [5 10 12 20 25 40 50 75 95];
for pi=1:length(Ps)
    [tmp pii] = min(abs(P-Ps(pi)));
    sprintf(['Percentile %3.1f is %5.1fm, %3.1fC, %3.1f degrees ' ...
             'latitude'],Ps(pi),-ZvP.mean(pii),TvP.mean(pii), ...
            YvP.mean(pii))
end

%%%% Plot total heat content and fits time series:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
set(gcf,'Position',[1 41 2560 1327.3]);
subplot(5,1,1);
plot(time,filter_field(CIN.N34,5,'-t'));
hold on;
plot(time,CIN.GMSST*10,'-r');
legend('Nino 3.4','Global Mean SST Anomalies * 10');
xlabel('Year');
ylabel('SST ($^\circ$C)');
xlim([yr1 yr1+nyrs]);

subplot(5,1,2);
plot(time,OHC,'-k');
hold on;
plot(time,OHCfull,'-b');
plot(time,OHC_cub,'-r');
legend('OHC dedrifted','OHC deseasoned','OHC cubic fix');
xlabel('Year');
ylabel('Ocean Heat Content (J)');
xlim([yr1 yr1+nyrs]);

subplot(5,1,3);
plot(time,OHC,'-k');
hold on;
plot(time,TvP.Hp(end,:),'--g');
plot(time,ZvP.Hp(end,:),'--c');
plot(time,YvP.Hp(end,:),'--y');
legend('OHC dedrifted',...
       'OHC from cumsum(TvP.Tp)','OHC from cumsum(ZvP.Tp)','OHC from cumsum(YvP.Tp)');
xlabel('Year');
ylabel('Ocean Heat Content (J)');
xlim([yr1 yr1+nyrs]);

subplot(5,1,4);
plot(time,CIN.AMOCfull,'-k');
hold on;
plot(time,filter_field(CIN.AMOCfull,5*12+1,'-t'),'-k','linewidth',2);
xlabel('Year');
ylabel('AMOC at $26^\circ$N (Sv)');
xlim([yr1 yr1+nyrs]);

subplot(5,1,5);
plot(time,CIN.WPOW,'-k');
hold on;
plot(time,filter_field(CIN.WPOW,5*12+1,'-t'),'-k','linewidth',2);
xlabel('Year');
ylabel('Total Wind Power (W)');
xlim([yr1 yr1+nyrs]);

%%% Historical OHC:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

his = load([baseMAT 'hisPP.mat']);
PI = load([baseMAT 'PIcontrolPP.mat']);
hisNATe1 = load([baseMAT 'hisNATe1PP.mat']);
hisAERe1 = load([baseMAT 'hisAERe1PP.mat']);
hisGHGe1 = load([baseMAT 'hisGHGe1PP.mat']);
% $$$ hisNATe2 = load([baseMAT 'hisNATe2PP.mat']);
% $$$ hisNATe3 = load([baseMAT 'hisNATe3PP.mat']);

% Raw OHC:
figure;
set(gcf,'defaulttextfontsize',15);
set(gcf,'defaultaxesfontsize',15);
subplot(3,1,1);
plot(PI.time+1850-100,PI.OHCfull,'-k');
hold on;
plot(PI.time+1850-100,PI.OHC_cub,'-g');
plot(his.time,his.OHCfull,'-b');
plot(hisNATe1.time,hisNATe1.OHCfull,'-r');
plot(hisAERe1.time,hisAERe1.OHCfull,'-m');
plot(hisGHGe1.time,hisGHGe1.OHCfull,'-c');
% $$$ plot(hisNATe2.time,hisNATe2.OHCfull,'-','color',[0 0.5 0]);
% $$$ plot(hisNATe3.time,hisNATe3.OHCfull,'-m');
legend('PI-control Raw OHC','PI-control Cubic Fit',...
       'Historial Raw OHC','Historical-NAT Raw OHC','Historical-AER Raw OHC');
xlim([1800 2100]);
xlabel('Year');
ylabel('OHC (J)');

subplot(3,1,2);
plot(PI.time+1850-100,PI.OHC,'-k');
hold on;
plot(his.time,his.OHCcubdd,'-b');
plot(hisNATe1.time,hisNATe1.OHCcubdd,'-r');
plot(hisAERe1.time,hisAERe1.OHCcubdd,'-m');
plot(hisGHGe1.time,hisGHGe1.OHCcubdd,'-c');
% $$$ plot(hisNATe2.time,hisNATe2.OHC,'-','color',[0 0.5 0]);
% $$$ plot(hisNATe3.time,hisNATe3.OHC,'-m');
% $$$ plot([1800 2350],[1 1]*2*std(PI.OHC),'--k');
% $$$ plot([1800 2350],[1 1]*-2*std(PI.OHC),'--k');
xlabel('Year');
ylabel('OHC Anomaly (J)');
legend('PI-control OHC cubic-detrend','Historial OHC cubic-detrend', ...
       'Historical-NAT OHC cubic-detrend','Historical-AER OHC cubic-detrend', ...
       'Historical-GHG OHC cubic-detrend');
xlim([1800 2100]);

subplot(3,1,3);
plot(PI.time+1850-100,PI.OHC,'-k');
hold on;
plot(his.time,his.OHC,'-b');
plot(hisNATe1.time,hisNATe1.OHC,'-r');
plot(hisAERe1.time,hisAERe1.OHC,'-r');
% $$$ plot(hisNATe2.time,hisNATe2.OHC,'-','color',[0 0.5 0]);
% $$$ plot(hisNATe3.time,hisNATe3.OHC,'-m');
% $$$ plot([1800 2350],[1 1]*2*std(PI.OHC),'--k');
% $$$ plot([1800 2350],[1 1]*-2*std(PI.OHC),'--k');
xlabel('Year');
ylabel('OHC (J)');
legend('PI-control OHC cubic-detrend','Historial OHC cubic-and-linear-detrend', ...
       'Historical-NAT OHC cubic-and-linear-detrend','Historical-AER OHC cubic-and-linear-detrend');
xlim([1800 2100]);

% $$$ figure;
% $$$ plot(PIstd.stdZvP.Tp*2,P,'-k');
% $$$ hold on;
% $$$ plot(PIstd.stdTvP.Tp*2,P,'-r');
% $$$ plot(PIstd.stdYvP.Tp*2,P,'-b');
% $$$ xlabel('2$\sigma$ variability amplitude ($^\circ$C)');
% $$$ ylabel('Percentile');
% $$$ legend('$\Theta(p_z)$','$\Theta(p_\Theta)$',['$\Theta(p_\' ...
% $$$                     'phi)$']);
% $$$ 
% Example time series:
pii = 20;
[tmp piii] = min(abs(P-pii));
figure;
plot(his.time,his.ZvP.Tp(piii,:),'-k','linewidth',2);
hold on;
plot(his.time,his.TvP.Tp(piii,:),'-r','linewidth',2);
plot(hisNATe1.time,hisNATe1.ZvP.Tp(piii,:),'--k','linewidth',2);
plot(hisNATe1.time,hisNATe1.TvP.Tp(piii,:),'--r','linewidth',2);
plot([his.time(1) his.time(end)],[1 1]*2*std(hisNATe1.ZvP.Tp(piii,:)),'--k');
plot([his.time(1) his.time(end)],[1 1]*-2*std(hisNATe1.ZvP.Tp(piii,:)),'--k');
plot([his.time(1) his.time(end)],[1 1]*2*std(hisNATe1.TvP.Tp(piii,:)),'--r');
plot([his.time(1) his.time(end)],[1 1]*-2*std(hisNATe1.TvP.Tp(piii,:)),'--r');
    
%%%%% Plot mean, climatology, std and trends:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%% Plot seasonal cycle and budget terms: %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% $$$ bvars = {'TENp','ADVp','FORp','RMIXp','VMIXp'};
% $$$ names = {'Tendency','Advection','Forcing', ...
% $$$          'Redi Mixing','Vertical Mixing'};
% $$$ bvars = {'ADVp','ADVGMp'};
% $$$ names = {'Total Advection','GM Advection'};

% Depth percentile:
Zcxs = [-0.2 0.01 0.2];
% $$$ Zcxs = [-1 0.05 1];
Zlab = '$\Theta(p_z)$';
zlim = [0 10];zlab = 'Depth Percentile $p_z$';

figure;
set(gcf,'Position',[1          36        1920         970]);
set(gcf,'defaulttextfontsize',15);
set(gcf,'defaultaxesfontsize',15);

[X,Y] = ndgrid(1:12,P);
subplot(3,2,1);
contourf(X,Y,ZvP.Tp_clim',[-1e50 Zcxs(1):Zcxs(2):Zcxs(3) 1e50],'linestyle','none');
caxis([Zcxs(1) Zcxs(3)]);
ylim(zlim);
ylabel(zlab);
set(gca,'ydir','reverse');
set(gca,'xticklabel',[]);
title(['Seasonal ' Zlab ' Anomalies']);
cb = colorbar;
ylabel(cb,'$^\circ$C');

for vi = 1:length(bvars)
    [X,Y] = ndgrid(1:12,Pe);
    eval(['var = ZvP.' bvars{vi} '_clim;']);
    subplot(3,2,1+vi);
    contourf(avg(X,2),avg(Y,2),var',[-1e50 Zcxs(1):Zcxs(2):Zcxs(3) 1e50],'linestyle','none');
    caxis([Zcxs(1) Zcxs(3)]);
    ylim(zlim);
    if (vi<4)
        set(gca,'xticklabel',[]);
    else
        xlabel('Month');
    end
    ylabel(zlab);
    cb = colorbar;
    ylabel(cb,'$^\circ$C');
    set(gca,'ydir','reverse');
    title(names{vi});
end
colormap(redblue);

% Temperature Percentile:
Tcxs = [-0.2 0.01 0.2];
tlim = [0 10];tlim = [95 100];
Tcxs = [-0.05 0.001 0.05];
tlim = [0 100];
% $$$ Tcxs = [-0.1 0.001 0.1];
% $$$ tlim = [0 100];
Tlab = '$\Theta(p_\Theta)$';
tlab = 'Temperature Percentile $p_\Theta$';

% $$$ figure;
% $$$ set(gcf,'Position',[1          36        1920         970]);
% $$$ set(gcf,'defaulttextfontsize',15);
% $$$ set(gcf,'defaultaxesfontsize',15);

[X,Y] = ndgrid(1:12,P);
subplot(3,2,3);
contourf(X,Y,TvP.Tp_clim',[-1e50 Tcxs(1):Tcxs(2):Tcxs(3) 1e50],'linestyle','none');
caxis([Tcxs(1) Tcxs(3)]);
ylim(tlim);
ylabel(tlab);
set(gca,'ydir','reverse');
set(gca,'xticklabel',[]);
title(['Seasonal ' Tlab ' Anomalies']);
cb = colorbar;
ylabel(cb,'$^\circ$C');

for vi = 1:length(bvars)
    [X,Y] = ndgrid(1:12,Pe);
    eval(['var = TvP.' bvars{vi} '_clim;']);
    subplot(3,2,1+vi);
    contourf(avg(X,2),avg(Y,2),var',[-1e50 Tcxs(1):Tcxs(2):Tcxs(3) 1e50],'linestyle','none');
    caxis([Tcxs(1) Tcxs(3)]);
    ylim(tlim);
    if (vi<4)
        set(gca,'xticklabel',[]);
    else
        xlabel('Month');
    end
    ylabel(tlab);
    cb = colorbar;
    ylabel(cb,'$^\circ$C');
    set(gca,'ydir','reverse');
    title(names{vi});
end
colormap(redblue);

% Latitude percentile:
Ycxs = [-0.1 0.005 0.1];Ylab = '$\Theta(p_\phi)$';
yylim = [0 100];
ylab = 'Latitude Percentile $p_\phi$';

figure;
set(gcf,'Position',[1          36        1920         970]);
set(gcf,'defaulttextfontsize',15);
set(gcf,'defaultaxesfontsize',15);

[X,Y] = ndgrid(1:12,P);
subplot(3,2,1);
contourf(X,Y,YvP.Tp_clim',[-1e50 Ycxs(1):Ycxs(2):Ycxs(3) 1e50],'linestyle','none');
caxis([Ycxs(1) Ycxs(3)]);
ylim(yylim);
ylabel(ylab);
set(gca,'xticklabel',[]);
title(['Seasonal ' Ylab ' Anomalies']);
cb = colorbar;
ylabel(cb,'$^\circ$C');

for vi = 1:length(bvars)
    [X,Y] = ndgrid(1:12,Pe);
    eval(['var = YvP.' bvars{vi} '_clim;']);
    subplot(3,2,1+vi);
    contourf(avg(X,2),avg(Y,2),var',[-1e50 Ycxs(1):Ycxs(2):Ycxs(3) 1e50],'linestyle','none');
    caxis([Ycxs(1) Ycxs(3)]);
    ylim(yylim);
    if (vi<4)
        set(gca,'xticklabel',[]);
    else
        xlabel('Month');
    end
    ylabel(ylab);
    cb = colorbar;
    ylabel(cb,'$^\circ$C');
    title(names{vi});
end
colormap(redblue);

%%% Straight comparison of amplitude of anomalies, along with 
%%% Mean curves:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
set(gcf,'Position',[23         208        1807         716]);
set(gcf,'defaulttextfontsize',15);
set(gcf,'defaultaxesfontsize',15);
subplot(1,3,1);
plot(ZvP.Tp_mean,P,'-k','linewidth',2);
hold on;
plot(TvP.Tp_mean,P,'-r','linewidth',2);
plot(YvP.Tp_mean,P,'-b','linewidth',2);
xlabel('Mean $\Theta$ ($^\circ$C)');
xlim([-2 30]);
ylabel('Percentile');
ylim([0 100]);
legend('$\Theta_z(p_z)$','$\Theta_\Theta(p_\Theta)$','$\Theta_\phi(p_\phi)$','Location','SouthEast');
set(gca,'ydir','reverse')
pos1 = [0.1759    0.1030    0.2134    0.8150];
set(gca,'Position',pos1);
text(-1.9,2,'(a)');

ax2_pos = pos1;
ax2_pos(1) = pos1(1)-pos1(1)/4;
ax2_pos(3) = 0.00001;
ax3_pos = ax2_pos;
ax3_pos(1) = pos1(1)-1.9*pos1(1)/4;
ax4_pos = ax2_pos;
ax4_pos(1) = pos1(1)-3.1*pos1(1)/4;

ax2 = axes('Position',ax2_pos,'XColor','r','YColor','r','FontSize',15);
set(ax2,'FontSize',15);
Zticks = [30 18 12:-2:4 3:-1:-1];
Pticks = zeros(size(Zticks));
for ii=1:length(Zticks)
    Pticks(ii) = interp1(TvP.Tp_mean,P,Zticks(ii),'linear');
end
set(ax2,'ytick',Pticks);
set(ax2,'yticklabel',Zticks);
ylabel(ax2,'Temperature ($^\circ$C)');
ylim(ax2,[0 100]);
set(ax2,'ydir','reverse');

ax3 = axes('Position',ax3_pos,'XColor','k','YColor','k','FontSize',15);
Zticks = [0:500:4000 5000];
Pticks = zeros(size(Zticks));
for ii=1:length(Zticks)
    Pticks(ii) = interp1(zofP_mean,P,-Zticks(ii),'linear');
end
Pticks(1) = 0;
set(ax3,'ytick',Pticks);
set(ax3,'yticklabel',Zticks);
ylabel(ax3,'Depth (m)');
ylim(ax3,[0 100]);
set(ax3,'ydir','reverse');

ax4 = axes('Position',ax4_pos,'XColor','b','YColor','b','FontSize',15);
Zticks = [-75:15:60];
Pticks = zeros(size(Zticks));
for ii=1:length(Zticks)
    Pticks(ii) = interp1(yofP_mean,P,Zticks(ii),'linear');
end
Pticks(1) = 0;
set(ax4,'ytick',Pticks);
set(ax4,'yticklabel',Zticks);
ylabel(ax4,'Latitude ($^\circ$N)');
ylim(ax4,[0 100]);
set(gca,'ydir','reverse');

subplot(1,3,2);
plot(std(ZvP.Tp,[],2),P,'-k','linewidth',2);
hold on;
plot(std(TvP.Tp,[],2),P,'-r','linewidth',2);
plot(std(YvP.Tp,[],2),P,'-b','linewidth',2);
xlabel('Standard Deviation of $\Theta$ ($^\circ$C)');
ylim([0 100]);
xlim([0 0.1]);
set(gca,'ydir','reverse');
set(gca,'yticklabel',[]);
set(gca,'Position',[0.4108    0.1030    0.2134    0.8150]);
text(0.001,2,'(b)');

% Averages:
text(0.04,30,['$\overline{\sigma_z}$ = ' sprintf('%3.3f',mean(std(ZvP.Tp,[],2))) '$^\circ$C'],'color','k');
text(0.04,35,['$\overline{\sigma_\Theta}$ = ' sprintf('%3.3f',mean(std(TvP.Tp,[],2))) '$^\circ$C'],'color','r');
text(0.04,40,['$\overline{\sigma_\phi}$ = ' sprintf('%3.3f',mean(std(YvP.Tp,[],2))) '$^\circ$C'],'color','b');

subplot(1,3,3);
plot(std(ZvP.Hp/1e22,[],2),P,'-k','linewidth',2);
hold on;
plot(std(TvP.Hp/1e22,[],2),P,'-r','linewidth',2);
plot(std(YvP.Hp/1e22,[],2),P,'-b','linewidth',2);
xlabel('Standard Deviation of $\mathcal{H}$ ($10^{22}$ J)');
ylim([0 100]);
xlim([0 1.6]);
set(gca,'ydir','reverse');
set(gca,'yticklabel',[]);
set(gca,'Position',[0.6504    0.1030    0.2134    0.8150]);
text(0.001,2,'(c)');

text(0.1,30,['$\overline{\sigma_z}$ = ' sprintf('%3.2f',mean(std(ZvP.Hp,[],2))/1e22) ' $\times10^{22}J$'],'color','k');
text(0.1,35,['$\overline{\sigma_\Theta}$ = ' sprintf('%3.2f',mean(std(TvP.Hp,[],2))/1e22) ' $\times10^{22}J$'],'color','r');
text(0.1,40,['$\overline{\sigma_\phi}$ = ' sprintf('%3.2f',mean(std(YvP.Hp,[],2))/1e22) ' $\times10^{22}J$'],'color','b');

%%%%%% Plot all anomalies P-t:  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% Choose whether to plot T or H:
% $$$ plot_H = 0;
% $$$ if (plot_H)
% $$$     ZvPar = ZvP.Hp;
% $$$     Zcxs = [-0.5e23 1e21 0.5e23];
% $$$     Zlab = '$H(p_z)$';
% $$$ 
% $$$     TvPar = TvP.Hp;
% $$$     Tcxs = [-0.5e23 1e21 0.5e23];
% $$$     Tlab = '$H(p_\Theta)$';
% $$$     
% $$$     yvar = YvP.Hp;
% $$$     ycxs = [-0.5e23 1e21 0.5e23];
% $$$     Ylab = '$H(p_\phi)$';
% $$$ else
ZvPar = ZvP.Tp;TvPar = TvP.Tp;yvar = YvP.Tp;
Zcxs = [-0.06 0.002 0.06];
% $$$ Zcxs = [-0.2 0.02 0.2];
Zlab = '$\Theta_z(p_z)$';
Tcxs = Zcxs;
Tlab = '$\Theta_\Theta(p_\Theta)$';
ycxs = Zcxs;
Ylab = '$\Theta_\phi(p_\phi)$';

zlim = [0 100];
tlim = [0 100];
yylim = [0 100];
xlims = [yr1 yr1+nyrs];

% $$$ % ENSO Focus settings:
% $$$ zlim = [0 10];
% $$$ tlim = [0 10];
% $$$ yylim = [40 85];
% $$$ xlims = [300 350];
% $$$ Zcxs = [-0.2 0.01 0.2];
% $$$ Tcxs = [-0.2 0.01 0.2];
% $$$ ycxs = [-0.2 0.01 0.2];

[tmp t1] = min(abs(time-xlims(1)));
[tmp t2] = min(abs(time-xlims(2)));

figure;
set(gcf,'Position',[1          40        1890         963]);
set(gcf,'defaulttextfontsize',15);
set(gcf,'defaultaxesfontsize',15);

x = mean(reshape(time(t1:t2),[12 (t2-t1+1)/12]),1)';
y = Z_YvPar;
v = squeeze(mean(reshape(ZvPar(:,t1:t2)',[12 (t2-t1+1)/12 PL]),1));
[X,Y] = ndgrid(x,y);
axes1 = subplot(3,1,1);
contourf(X,Y,v,[-1e50 Zcxs(1):Zcxs(2):Zcxs(3) 1e50],'linestyle','none');
caxis([Zcxs(1) Zcxs(3)]);
ylim(zlim);
xlim(xlims);
% $$$     xlabel('Year');
set(gca,'xticklabel',[]);
ylabel(zlab);
% $$$     cb = colorbar;
% $$$     ylabel(cb,'Temperature Anomalies ($^\circ$C)','Interpreter','latex');
text(xlims(1)+(xlims(2)-xlims(1))*0.01,0.93*zlim(2),['(a) ' Zlab ' anomalies'],'Backgroundcolor','w');
set(gca,'Position',[0.13 0.7 0.75 0.28]);
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

x = mean(reshape(time(t1:t2),[12 (t2-t1+1)/12]),1)';
y = T_YvPar;
v = squeeze(mean(reshape(TvPar(:,t1:t2)',[12 (t2-t1+1)/12 PL]),1));
[X,Y] = ndgrid(x,y);
axes2 = subplot(3,1,2);
contourf(X,Y,v,[-1e50 Tcxs(1):Tcxs(2):Tcxs(3) 1e50],'linestyle','none');
caxis([Tcxs(1) Tcxs(3)]);
ylim(tlim);
xlim(xlims);
colormap(redblue);
% $$$     xlabel('Year');
set(gca,'xticklabel',[]);
ylabel(tlab);
cb = colorbar;
ylabel(cb,'Temperature Anomalies ($^\circ$C)','Interpreter','latex');
text(xlims(1)+(xlims(2)-xlims(1))*0.01,0.83*tlim(2),['(b) ' Tlab ' anomalies'],'Backgroundcolor','w');
set(gca,'Position',[0.13 0.38 0.75 0.28]);
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
    Zticks = [18 12 8:-2:4 3:-1:0];
% $$$         Zticks = [30 22 18 16:-2:8];
% $$$         Zticks = [2:-0.25:-0.5 -2];
    Pticks = zeros(size(Zticks));
    for ii=1:length(Zticks)
        Pticks(ii) = interp1(TvP.Tp_mean,P,Zticks(ii),'linear');
    end
    Pticks(end) = 100;
    set(ax2,'ytick',Pticks);
    set(ax2,'yticklabel',Zticks);
    ylabel(ax2,'Temperature ($^\circ$C)');
    ylim(ax2,zlim);
    set(ax2,'ydir','reverse');
end

x = mean(reshape(time(t1:t2),[12 (t2-t1+1)/12]),1)';
y = y_YvPar;
v = squeeze(mean(reshape(yvar(:,t1:t2)',[12 (t2-t1+1)/12 PL]),1));
[X,Y] = ndgrid(x,y);
axes3 = subplot(3,1,3);
contourf(X,Y,v,[-1e50 ycxs(1):ycxs(2):ycxs(3) 1e50],'linestyle','none');
caxis([ycxs(1) ycxs(3)]);
ylim(yylim);
xlim(xlims);
xlabel('Year');
ylabel(ylab);
% $$$     cb = colorbar;
% $$$     ylabel(cb,'Temperature Anomalies ($^\circ$C)','Interpreter','latex');
text(xlims(1)+(xlims(2)-xlims(1))*0.01,7,['(c) ' Ylab ' anomalies'],'Backgroundcolor','w');
set(gca,'Position',[0.13 0.07 0.75 0.28]);
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

% Inset zoom plots:
ZvPar = ZvP.Tp;TvPar = TvP.Tp;yvar = YvP.Tp;
Zcxs = [-0.2 0.01 0.2];Tcxs = Zcxs;ycxs = Zcxs;
zlim = [0 10];
tlim = [0 10];
yylim = [40 80];
xlims = [305 325];

axes(axes1);
hold on;
plot([xlims(1) xlims(1) xlims(2) xlims(2) xlims(1)],...
     [zlim(1) zlim(2) zlim(2) zlim(1) zlim(1)],'-k');
axes(axes2);
hold on;
plot([xlims(1) xlims(1) xlims(2) xlims(2) xlims(1)],...
     [tlim(1) tlim(2) tlim(2) tlim(1) tlim(1)],'-k');
axes(axes3);
hold on;
plot([xlims(1) xlims(1) xlims(2) xlims(2) xlims(1)],...
     [yylim(1) yylim(2) yylim(2) yylim(1) yylim(1)],'-k');

[tmp t1] = min(abs(time-xlims(1)));
[tmp t2] = min(abs(time-xlims(2)));
x = time(t1:t2);
y = Z_YvPar;
v = ZvPar(:,t1:t2)';
[X,Y] = ndgrid(x,y);
axes1i = axes('Position',[0.655    0.7418    0.2143    0.12]);
axes(axes1i);
contourf(X,Y,v,[-1e50 Zcxs(1):Zcxs(2):Zcxs(3) 1e50],'linestyle','none');
caxis([Zcxs(1) Zcxs(3)]);
ylim(zlim);
xlim(xlims);
set(gca,'xtick',[xlims(1) mean(xlims) xlims(2)]);
cb = colorbar;
set(gca,'ydir','reverse');
set(gca,'Position',[0.655    0.7418    0.2143    0.12]);

hold on;
plot(time(t1:t2),-(CIN.N34(t1:t2)-mean(CIN.N34))+7.5,'-k');
plot([time(t1) time(t2)],[7.5 7.5],'--k');

if (~remap_to_T)
    ax1=gca;
    ax1_pos = ax1.Position; % position of first axes
    ax1_pos(1) = ax1_pos(1)-0.02;
    ax1_pos(3) = 0.00001;
    ax2 = axes('Position',ax1_pos,...
               'Color','none');
    set(ax2,'FontSize',15);
    Zticks = [0:100:300];
% $$$         Zticks = [0:50:800];
% $$$         Zticks = [3000:250:5500];
    Pticks = zeros(size(Zticks));
    for ii=1:length(Zticks)
        Pticks(ii) = interp1(zofP_mean,P,-Zticks(ii),'linear');
    end
    Pticks(1) = 0;
    set(ax2,'ytick',Pticks);
    set(ax2,'yticklabel',Zticks);
% $$$     ylabel(ax2,'Depth (m)');
    ylim(ax2,zlim);
    set(ax2,'ydir','reverse');
end        

% $$$ axes1iN = axes('Position',[0.655    0.87    0.2143    0.07]);
% $$$ axes(axes1iN);

y = T_YvPar;
v = TvPar(:,t1:t2)';
[X,Y] = ndgrid(x,y);
axes2i = axes('Position',[0.655    0.4448    0.2143    0.12]);
axes(axes2i);
contourf(X,Y,v,[-1e50 Tcxs(1):Tcxs(2):Tcxs(3) 1e50],'linestyle','none');
caxis([Tcxs(1) Tcxs(3)]);
ylim(tlim);
xlim(xlims);
set(gca,'xtick',[xlims(1) mean(xlims) xlims(2)]);
% $$$ cb = colorbar;
set(gca,'ydir','reverse');
set(gca,'Position',[0.655    0.4448    0.2143    0.12]);
if (~remap_to_T)
    set(gca,'ydir','reverse');
    ax1=gca;
    ax1_pos = ax1.Position; % position of first axes
    ax1_pos(1) = ax1_pos(1)-0.02;
    ax1_pos(3) = 0.00001;
    ax2 = axes('Position',ax1_pos,...
               'Color','none');
    set(ax2,'FontSize',15);
    Zticks = [30 20 14 12];
    Pticks = zeros(size(Zticks));
    for ii=1:length(Zticks)
        Pticks(ii) = interp1(TvP.Tp_mean,P,Zticks(ii),'linear');
    end
    set(ax2,'ytick',Pticks);
    set(ax2,'yticklabel',Zticks);
% $$$     ylabel(ax2,'Temperature ($^\circ$C)');
    ylim(ax2,zlim);
    set(ax2,'ydir','reverse');
end

y = y_YvPar;
v = yvar(:,t1:t2)';
[X,Y] = ndgrid(x,y);
axes3i = axes('Position',[0.655    0.1271    0.2143    0.12]);
axes(axes3i);
contourf(X,Y,v,[-1e50 ycxs(1):ycxs(2):ycxs(3) 1e50],'linestyle','none');
caxis([ycxs(1) ycxs(3)]);
ylim(yylim);
xlim(xlims);
set(gca,'xtick',[xlims(1) mean(xlims) xlims(2)]);
cb = colorbar;
set(gca,'Position',[0.655    0.1271    0.2143    0.12]);
if (~remap_to_T)
    ax1=gca;
    ax1_pos = ax1.Position; % position of first axes
    ax1_pos(1) = ax1_pos(1)-0.02;
    ax1_pos(3) = 0.00001;
    ax2 = axes('Position',ax1_pos,...
               'Color','none');
    set(ax2,'FontSize',15);
    Zticks = [-15 0 15];%:15:60];
    Pticks = zeros(size(Zticks));
    for ii=1:length(Zticks)
        Pticks(ii) = interp1(yofP_mean,P,Zticks(ii),'linear');
    end
    set(ax2,'ytick',Pticks);
    set(ax2,'yticklabel',Zticks);
% $$$     ylabel(ax2,'Latitude ($^\circ$N)');
    ylim(ax2,yylim);
end

colormap(redblue);


% $$$ % ENSO Time series to go with that:
% $$$ figure;
% $$$ set(gcf,'Position',[1 41 2560 1330]);
% $$$ set(gcf,'defaulttextfontsize',15);
% $$$ set(gcf,'defaultaxesfontsize',15);
% $$$ subplot(2,1,1);
% $$$ hold on;
% $$$ plot(time(t1:t2),GMSST(t1:t2)'*10,'-r');
% $$$ plot(time(t1:t2),OHC(t1:t2)/1e22,'-','color',[0.5 0.5 0.5]);
% $$$ legend('Nino 3.4','Global SST * 10','Global OHC anomaly ($/10^{22} J$)');

%%%%%% Plot P-t of anomalies + budget terms %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xlims = [yr1 yr1+nyrs];
xlims = [320 350];

[tmp t1] = min(abs(time-xlims(1)));
[tmp t2] = min(abs(time-xlims(2)));
bvars = {'TENp','ADVp','FORp','RMIXp','VMIXp'};
names = {'Tendency','Advection','Forcing','Redi Mixing', ...
         'Vertical Mixing'};

poss = [0.0539 0.7073 0.4 0.26; ...
        0.5168 0.7073 0.4 0.26; ...
        0.0539 0.3979 0.4 0.26; ...
        0.5168 0.3979 0.4 0.26; ...
        0.0539 0.0792 0.4 0.26; ...
        0.5168 0.0792 0.4 0.26;];

% Take time derivative?
tder = 1; 

% Depth percentile:
Zlab = '$\Theta(p_z)$';
    zlim = [0 15];zlab = 'Depth Percentile $p_z$';    
if (tder)
    Zcxs = [-0.25e-7 0.1e-8 0.25e-7];
else
    Zcxs = [-0.3 0.01 0.3];
end

figure;
set(gcf,'Position',[1          36        1920         970]);
set(gcf,'defaulttextfontsize',15);
set(gcf,'defaultaxesfontsize',15);
if (~tder)
[X,Y] = ndgrid(time(t1:t2),P);
subplot(3,2,1);
contourf(X,Y,ZvP.Tp(:,t1:t2)',[-1e50 Zcxs(1):Zcxs(2):Zcxs(3) 1e50],'linestyle','none');
caxis([Zcxs(1) Zcxs(3)]);
ylim(zlim);
ylabel(zlab);
set(gca,'ydir','reverse');
set(gca,'xticklabel',[]);
title([Zlab ' Anomalies']);
set(gca,'Position',poss(1,:));
end
for vi = 1:length(bvars)
    [X,Y] = ndgrid(time(t1:t2),Pe);
    eval(['var = ZvP.' bvars{vi} '(:,t1:t2);']);

    if (tder)
        var = diff(var,[],2)./repmat(DT_A(t1+1:t2)',[PL 1]);
        var = cat(2,zeros(PL,1),var);
    end

    subplot(3,2,1+vi);
    contourf(avg(X,2),avg(Y,2),var',[-1e50 Zcxs(1):Zcxs(2):Zcxs(3) 1e50],'linestyle','none');
    caxis([Zcxs(1) Zcxs(3)]);
    ylim(zlim);
    if (vi<4)
        set(gca,'xticklabel',[]);
    else
        xlabel('Year');
    end
    if (vi == 2 | vi == 4)
        ylabel(zlab);
    else
        cb = colorbar;
        ylabel(cb,'$^\circ$C');
    end
    set(gca,'ydir','reverse');
    title(names{vi});
    set(gca,'Position',poss(1+vi,:));
end
colormap(redblue);

% Temperature Percentile:
if (tder)
    Tcxs = [-0.25e-7 0.1e-8 0.25e-7];
else
    Tcxs = [-0.2 0.01 0.2];
end
tlim = [0 5];
Tlab = '$\Theta(p_\Theta)$';
tlab = 'Temperature Percentile $p_\Theta$';

figure;
set(gcf,'Position',[1          36        1920         970]);
set(gcf,'defaulttextfontsize',15);
set(gcf,'defaultaxesfontsize',15);

if (~tder)
[X,Y] = ndgrid(time(t1:t2),P);
subplot(3,2,1);
contourf(X,Y,TvP.Tp(:,t1:t2)',[-1e50 Tcxs(1):Tcxs(2):Tcxs(3) 1e50],'linestyle','none');
caxis([Tcxs(1) Tcxs(3)]);
ylim(tlim);
ylabel(tlab);
set(gca,'ydir','reverse');
set(gca,'xticklabel',[]);
title([Tlab ' Anomalies']);
set(gca,'Position',poss(1,:));
end
for vi = 1:length(bvars)
    [X,Y] = ndgrid(time(t1:t2),Pe);
    eval(['var = TvP.' bvars{vi} '(:,t1:t2);']);
    if (tder)
        var = diff(var,[],2)./repmat(DT_A(t1+1:t2)',[PL 1]);
        var = cat(2,zeros(PL,1),var);
    end
    subplot(3,2,1+vi);
    contourf(avg(X,2),avg(Y,2),var',[-1e50 Tcxs(1):Tcxs(2):Tcxs(3) 1e50],'linestyle','none');
    caxis([Tcxs(1) Tcxs(3)]);
    ylim(tlim);
    if (vi<4)
        set(gca,'xticklabel',[]);
    else
        xlabel('Year');
    end
    if (vi == 2 | vi == 4)
        ylabel(tlab);
    else
        cb = colorbar;
        ylabel(cb,'$^\circ$C');
    end
    set(gca,'ydir','reverse');
    title(names{vi});
    set(gca,'Position',poss(1+vi,:));
end
colormap(redblue);

% Latitude percentile:
Ycxs = [-0.1 0.005 0.1];Ylab = '$\Theta(p_\phi)$';
yylim = [0 100];
ylab = 'Latitude Percentile $p_\phi$';

figure;
set(gcf,'Position',[1          36        1920         970]);
set(gcf,'defaulttextfontsize',15);
set(gcf,'defaultaxesfontsize',15);

[X,Y] = ndgrid(time(t1:t2),P);
subplot(3,2,1);
contourf(X,Y,YvP.Tp(:,t1:t2)',[-1e50 Ycxs(1):Ycxs(2):Ycxs(3) 1e50],'linestyle','none');
caxis([Ycxs(1) Ycxs(3)]);
ylim(yylim);
ylabel(ylab);
set(gca,'ydir','reverse');
set(gca,'xticklabel',[]);
title([Ylab ' Anomalies']);
set(gca,'Position',poss(1,:));

for vi = 1:length(bvars)
    [X,Y] = ndgrid(time(t1:t2),Pe);
    eval(['var = YvP.' bvars{vi} '(:,t1:t2);']);
    if (tder)
        var = diff(var,[],2)./repmat(DT_A(t1+1:t2)',[PL 1]);
        var = cat(2,zeros(PL,1),var);
    end
    subplot(3,2,1+vi);
    contourf(avg(X,2),avg(Y,2),var',[-1e50 Ycxs(1):Ycxs(2):Ycxs(3) 1e50],'linestyle','none');
    caxis([Ycxs(1) Ycxs(3)]);
    ylim(yylim);
    if (vi<4)
        set(gca,'xticklabel',[]);
    else
        xlabel('Year');
    end
    if (vi == 2 | vi == 4)
        ylabel(ylab);
    else
        cb = colorbar;
        ylabel(cb,'$^\circ$C');
    end
    set(gca,'ydir','reverse');
    title(names{vi});
    set(gca,'Position',poss(1+vi,:));
end
colormap(redblue);


%%%%% Spectral analysis: %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs = 12;
nw = 10;
[pxx,f] = pmtm(ZvP.Tp(15,:)',nw,[],fs);
fL = length(f);
df = f(2)-f(1);

TorH = 1;

if (TorH)
vars = {'Tp','Tp','TENp','ADVp','FORp','RMIXp','VMIXp'};
else
vars = {'Hp','Hp','TEN_c','ADV_c','FOR_c','RMIX_c','VMIX_c'};
end
tder = [0 1 1 1 1 1 1]; % Take time-derivative of budget terms first

tnames = {'Temperature','Depth','Latitude'};
tcols = {'r','k','b'};

for ti=1:length(typs)
    for vi=1:length(vars)
        vars{vi}
        eval(['var = ' typs{ti} 'vP.' vars{vi} ';']);
        if (tder(vi))
                var = diff(var,[],2)./repmat(DT_A(2:end)',[length(var(:,1)) 1]);
                var = cat(2,zeros(length(var(:,1)),1),var);
        end
        lg = length(var(:,1));
        spec = zeros(fL,lg);
        for pi=1:lg
            [spec(:,pi),~] = pmtm(var(pi,:)',nw,[],fs);
        end
        if (tder(vi))
            eval([typs{ti} 'vP.' vars{vi} 'd_spec = spec;']);
        else
            eval([typs{ti} 'vP.' vars{vi} '_spec = spec;']);
        end            
    end
end

poss = [0.07    0.55    0.42    0.4; ...
        0.55     0.55    0.42    0.4; ...
        0.07    0.08    0.42    0.4; ...
        0.55     0.08    0.42    0.4];

pranges = [10 40;
           10 40;
           0 40;]
pranges = [0 100;
           0 100;
           0 100;]
prangesi = pranges;
for pi=1:length(pranges(:))
    [tmp prangesi(pi)] = min(abs(P-pranges(pi)));
end

% Average over frequency ranges and print
sfr = [5 1; ...
       100 5];
for si = 1:length(sfr(:,1))
    [tmp ind1] = min(abs(f-1/sfr(si,1)));
    [tmp ind2] = min(abs(f-1/sfr(si,2)));
    [sprintf('Z var (%01d-%01d years) = %3.2f',sfr(si,2),sfr(si,1), ...
             10^(4)*mean(sum(ZvP.Tp_spec(ind1:ind2,:)*df,1),2)) '$10^{-4}^\circ$C$^2$']
    [sprintf('T var (%01d-%01d years) = %3.2f',sfr(si,2),sfr(si,1), ...
             10^(4)*mean(sum(TvP.Tp_spec(ind1:ind2,:)*df,1),2)) '$10^{-4}^\circ$C$^2$']
    [sprintf('Y var (%01d-%01d years) = %3.2f',sfr(si,2),sfr(si,1), ...
             10^(4)*mean(sum(YvP.Tp_spec(ind1:ind2,:)*df,1),2)) '$10^{-4}^\circ$C$^2$']
end

% Spectral average all terms:
figure;
set(gcf,'Position',[1 1 1213.3*2 614.7*2]);
subplot(2,2,1);
if (TorH)
plot(f,log10(mean(ZvP.Tp_spec(:,prangesi(1,1):prangesi(1,2)),2)),'-k','linewidth',2);
hold on;
plot(f,log10(mean(TvP.Tp_spec(:,prangesi(2,1):prangesi(2,2)),2)),'-r','linewidth',2);
plot(f,log10(mean(YvP.Tp_spec(:,prangesi(3,1):prangesi(3,2)),2)),'-b','linewidth',2);
ylim([-5 -2]);
else
plot(f,log10(mean(ZvP.Hp_spec(:,prangesi(1,1):prangesi(1,2)),2)),'-k','linewidth',2);
hold on;
plot(f,log10(mean(TvP.Hp_spec(:,prangesi(2,1):prangesi(2,2)),2)),'-r','linewidth',2);
plot(f,log10(mean(YvP.Hp_spec(:,prangesi(3,1):prangesi(3,2)),2)),'-b','linewidth',2);
end
xlim([1/250 1]);
xlab = [300 100 50 10 7 5 3 2 1];
set(gca,'xscale','log');
set(gca,'xtick',1./xlab);
set(gca,'xticklabel',[]);
% $$$ xlabel('Period (years)');
legend('$\Theta_z$','$\Theta_\Theta$','$\Theta_\phi$','FontSize',10);
ylabel('Power $\log_{10}$($^\circ$C$^2$ year)');
grid on;
% $$$ title('ACCESS-CM2 PI-control Percentile-Averaged Spectra');
set(gca,'Position',poss(1,:));

% Spectral average (budget terms):
colors = {'m','b','k','r',[0 0.5 0]};
names = {{'$\partial\Theta_\Theta/\partial t$','Tendency','Numerical Mixing','Forcing','Neutral Mixing','Vertical Mixing'}, ...
         {'$\partial\Theta_z/\partial t$','Tendency','Advection','Forcing','Neutral Mixing','Vertical Mixing'}, ...
         {'$\partial\Theta_\phi/\partial t$','Tendency','Advection','Forcing','Neutral Mixing','Vertical Mixing'}};
for ti=1:length(typs)
% $$$     figure;
% $$$     set(gcf,'Position',[421.7 537.7 1213.3 614.7]);
    subplot(2,2,ti+1);
    if (TorH)
        eval(['var = ' typs{ti} 'vP.Tpd_spec;']);
    else
        eval(['var = ' typs{ti} 'vP.Hpd_spec;']);
    end
    plot(f,log10(mean(var(:,prangesi(ti,1):prangesi(ti,2)),2)),'--','linewidth',2,'color','m');
    hold on;
    for vi=3:length(vars)
        eval(['var = ' typs{ti} 'vP.' vars{vi} 'd_spec;']);
        plot(f,log10(mean(var(:,prangesi(ti,1):prangesi(ti,2)),2)),'-','linewidth',2,'color',colors{vi-2});
        hold on;
    end
xlim([1/250 1]);    
xlab = [300 100 50 10 7 5 3 2 1];
set(gca,'xscale','log');
grid on;
set(gca,'xtick',1./xlab);
if (ti>=2)
set(gca,'xticklabel',xlab);
xlabel('Period (years)');
else
    set(gca,'xticklabel',[]);
end
legend(names{ti},'FontSize',12,'Location','SouthWest');
if (TorH)
    ylim([-21 -17.25]);
else
    ylim([24.5 29.5]);
end
if (1)
    ylabel('Power $\log_{10}$($^\circ$C$^2$s$^{-2}$ year)');
% $$$     legend({'d$\Theta(p)$/dt',names{:}});
else
    ylabel('Power $\log_{10}$($^\circ$C$^2$ year)');
% $$$     legend({'$\Theta(p)$',names{:}});
end
set(gca,'Position',poss(ti+1,:));
title([tnames{ti} '-percentile budget'],'Color',tcols{ti});
end


%%%% Time-of-emergence
his = load([baseMAT 'hisPP.mat']);
hisNATe1 = load([baseMAT 'hisNATe1PP.mat']);

typs = {'T','Z','Y'};
for ti = 1:length(typs)
    tn = typs{ti};
    for pi=1:PL
        eval(['stdv = std(hisNATe1.' tn 'vP.Tp(pi,:));']);
        eval(['varv = his.' tn 'vP.Tp(pi,:);']);
        ind = find(varv<=2*stdv,1,'last');
        if (ind==length(varv))
            toe = NaN;
        else
            toe = time(ind);
        end
        eval([tn 'vP.ToE(pi) = toe;']);
    end
end

figure;
set(gcf,'Position',[23         208        1807         716]);
set(gcf,'defaulttextfontsize',15);
set(gcf,'defaultaxesfontsize',15);
subplot(1,3,1);
plot(std(hisNATe1.ZvP.Tp,[],2),P,'-k','linewidth',2);
hold on;
plot(std(hisNATe1.TvP.Tp,[],2),P,'-r','linewidth',2);
plot(std(hisNATe1.YvP.Tp,[],2),P,'-b','linewidth',2);
xlabel('$\sigma$ hist-NAT ($^\circ$C)');
xlim([0 0.1]);
ylabel('Percentile');
ylim([0 100]);
legend('$\Theta_z(p_z)$','$\Theta_\Theta(p_\Theta)$','$\Theta_\phi(p_\phi)$','Location','SouthEast');
set(gca,'ydir','reverse')
pos1 = [0.1759    0.1030    0.2134    0.8150];
set(gca,'Position',pos1);
text(28,5,'(a)');

ax2_pos = pos1;
ax2_pos(1) = pos1(1)-pos1(1)/4;
ax2_pos(3) = 0.00001;
ax3_pos = ax2_pos;
ax3_pos(1) = pos1(1)-1.9*pos1(1)/4;
ax4_pos = ax2_pos;
ax4_pos(1) = pos1(1)-3.1*pos1(1)/4;

ax2 = axes('Position',ax2_pos,'XColor','r','YColor','r','FontSize',15);
set(ax2,'FontSize',15);
Zticks = [30 18 12:-2:4 3:-1:0];
Pticks = zeros(size(Zticks));
for ii=1:length(Zticks)
    Pticks(ii) = interp1(TvP.Tp_mean,P,Zticks(ii),'linear');
end
set(ax2,'ytick',Pticks);
set(ax2,'yticklabel',Zticks);
ylabel(ax2,'Temperature ($^\circ$C)');
ylim(ax2,[0 100]);
set(ax2,'ydir','reverse');

ax3 = axes('Position',ax3_pos,'XColor','k','YColor','k','FontSize',15);
Zticks = [0:500:4000 5000];
Pticks = zeros(size(Zticks));
for ii=1:length(Zticks)
    Pticks(ii) = interp1(zofP_mean,P,-Zticks(ii),'linear');
end
Pticks(1) = 0;
set(ax3,'ytick',Pticks);
set(ax3,'yticklabel',Zticks);
ylabel(ax3,'Depth (m)');
ylim(ax3,[0 100]);
set(ax3,'ydir','reverse');

ax4 = axes('Position',ax4_pos,'XColor','b','YColor','b','FontSize',15);
Zticks = [-75:15:60];
Pticks = zeros(size(Zticks));
for ii=1:length(Zticks)
    Pticks(ii) = interp1(yofP_mean,P,Zticks(ii),'linear');
end
Pticks(1) = 0;
set(ax4,'ytick',Pticks);
set(ax4,'yticklabel',Zticks);
ylabel(ax4,'Latitude ($^\circ$N)');
ylim(ax4,[0 100]);
set(gca,'ydir','reverse');

subplot(1,3,2);
plot(ZvP.ToE,P,'-k','linewidth',2);
hold on;
plot(TvP.ToE,P,'-r','linewidth',2);
plot(YvP.ToE,P,'-b','linewidth',2);
ylim([0 100]);
xlabel('Year-of-emergence');
xlim([1920 2020]);
set(gca,'ydir','reverse');
set(gca,'yticklabel',[]);
set(gca,'Position',[0.4108    0.1030    0.2134    0.8150]);
text(1860,5,'(b)');




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
    
    
    
%%% Cross correlation maps:    

Tzcc = zeros(PL,PL);
Tycc = zeros(PL,PL);
zycc = zeros(PL,PL);

for Pi = 1:PL
    TS = TvP.Tp(Pi,:)';
    TS = TS/sqrt(sum(TS.^2));
    Tzcc(Pi,:) = ZvP.Tp*TS./sqrt(sum(ZvP.Tp.^2,2));
    Tycc(Pi,:) = YvP.Tp*TS./sqrt(sum(YvP.Tp.^2,2));
    
    TS = ZvP.Tp(Pi,:)';
    TS = TS/sqrt(sum(TS.^2));
    zycc(Pi,:) = YvP.Tp*TS./sqrt(sum(YvP.Tp.^2,2));
end

figure;
[X,Y] = ndgrid(P,P);
subplot(2,2,1);
pcolPlot(X,Y,Tzcc);
hold on;
contour(X,Y,Tzcc,[-0.4 0.4],'-k');
xlabel('Temperature Percentile');
ylabel('Depth Percentile');
caxis([-1 1]);
colorbar;
subplot(2,2,2);
pcolPlot(X,Y,Tycc);
hold on;
contour(X,Y,Tycc,[-0.4 0.4],'-k');
xlabel('Temperature Percentile');
ylabel('Latitude Percentile');
colorbar;
caxis([-1 1]);
subplot(2,2,3);
pcolPlot(X,Y,zycc);
hold on;
contour(X,Y,zycc,[-0.4 0.4],'-k');
xlabel('Depth Percentile');
ylabel('Latitude Percentile');
colorbar;
caxis([-1 1]);
colormap('redblue');

%%% Tyz plotting:

% Regressions:
ts = CIN.N34;
label = 'ENSO';

[tmp ind] = min(abs(P-20));
ts = ZvP.Tp(ind,:)';

ts = mean(ZvP.Tp(10:30,:),1)';
ts = mean(TvP.Tp(90:100,:),1)';
ts = filter_field(CIN.AMOC,121,'-t');
ts(isnan(ts)) = 0;

ts = mean(TvP.Tp(18:22,:),1)';
ts_dec = filter_field(ts,30*12+1,'-t');
ts_dec(isnan(ts_dec)) = 0;
ts_lp = ts-ts_dec;
ts_lp(isnan(ts_lp)) = 0;
ts = ts_lp;

Tyz_reg = reshape(reshape(Tyz,[yL*zL tL])*ts/(sum(ts.^2)),[yL zL]);

figure;
[X,Y] = ndgrid(latv,Z);
subplot(2,2,1);
pcolPlot(X,Y,Tyz_reg);
set(gca,'ydir','reverse');
ylim([0 4000]);
% $$$ caxis([-1 1]);%-0.5 0.5]);
colormap('redblue');

% Time-average:

[tmp t1] = min(abs(time-225));
[tmp t2] = min(abs(time-275));

figure;
[X,Y] = ndgrid(latv,Z);
pcolPlot(X,Y,mean(Tyz(:,:,t1:t2),3));
set(gca,'ydir','reverse');
ylim([0 4000]);
% $$$ caxis([-0.3 0.3]);%-0.5 0.5]);
colormap('redblue');

% Nice decadal plot:
figure;
set(gcf,'Position',[1          36        1920         970]);
set(gcf,'defaulttextfontsize',15);
set(gcf,'defaultaxesfontsize',15);
ts1 = mean(ZvP.Tp(10:30,:),1)';
ts2 = mean(TvP.Tp(10:30,:),1)';
ts1 = filter_field(ts1,12*25+1,'-t');
ts1(isnan(ts1)) = 0;
ts2 = filter_field(ts2,12*25+1,'-t');
ts2(isnan(ts2)) = 0;
subplot(5,1,1);
plot(time,ts1,'-k');
hold on;
plot(time,ts2,'-r');
legend('$\overline{\Theta_z(10\%<p_z<30\%,t)}$','$\overline{\Theta_\Theta(10\%<p_\Theta<30\%,t)}$');
hold on;
plot([225 225],[-0.02 0.02],'--k');
plot([275 275],[-0.02 0.02],'--k');
xlim([50 600]);
ylim([-0.03 0.03]);
set(gca,'Position',[0.1300    0.8172    0.7058    0.1672]);
xlabel('Year');
ylabel('$\Theta$ Anomaly $(^\circ$C)');

Tyz_reg = reshape(reshape(Tyz,[yL*zL tL])*ts1/(sum(ts1.^2)),[yL zL]);
[X,Y] = ndgrid(latv,Z);

subplot(5,1,[2 3]);
contourf(X,Y,Tyz_reg,[-1e50 -0.8e22:0.025e22:0.8e22 1e50],'linestyle','none');
hold on;
[c,h] = contour(X,Y,Tyz_mean,[-2:2:40],'-k');
clabel(c,h);
caxis([-0.8 0.8]*1e22);
cb = colorbar;
ylabel(cb,'J/m/$^\circ$-latitude/$^\circ$C');
set(gca,'ydir','reverse');
% $$$ xlabel('Latitude ($^\circ$N)');
ylabel('Depth (m)');
ylim([0 4000]);
xlim([-80 70]);

[tmp t1] = min(abs(time-225));
[tmp t2] = min(abs(time-275));

subplot(5,1,[4 5]);
contourf(X,Y,mean(Tyz(:,:,t1:t2),3),[-1e50 -1e20:0.025e20:1e20 1e50],'linestyle','none');
hold on;
[c,h] = contour(X,Y,Tyz_mean,[-2:2:40],'-k');
clabel(c,h);
caxis([-0.1 0.1]*1e21);
cb = colorbar;
ylabel(cb,'J/m/$^\circ$-latitude');
set(gca,'ydir','reverse');
xlabel('Latitude ($^\circ$N)');
ylabel('Depth (m)');
ylim([0 4000]);
xlim([-80 70]);
colormap('redblue');






%%%% Plot Lag Regressions:
    
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


    % ENSO:
    ts = CIN.N34-mean(CIN.N34);
    ts = ts/std(ts);
    label = 'Ni\~{n}o 3.4';
    lags = [-12*2:1:12*2]; %months lag
    zlim = [0 7];
    tlim = [0 7];
    yylim = [40 85];
    Zcxs = [-0.1 0.001 0.1]; Tcxs = Zcxs;Ycxs=[-0.02 0.0005 0.02];
    
    ll = length(lags);

    ZvPar_lr = zeros(PL,ll);
    TvPar_lr = zeros(PL,ll);
    yvar_lr = zeros(PL,ll);
    for ii=1:ll
        lag = lags(ii);
        TS = zeros(size(ts));
        if (lag<0)
            TS(1:(end+lag)) = ts((-lag+1):end);
            TS((end+lag+1):end) = 0;
        elseif (lag == 0)
            TS = ts;
        else
            TS(1:(lag)) = 0;
            TS((lag+1):end) = ts(1:(end-lag));
        end
        
        ZvPar_lr(:,ii) = ZvPar*TS/(sum(TS.^2));
        TvPar_lr(:,ii) = TvPar*TS/(sum(TS.^2));
        yvar_lr(:,ii) = yvar*TS/(sum(TS.^2));
    end

    figure;
    set(gcf,'Position',[3          40        1130         963]);
    set(gcf,'defaulttextfontsize',15);
    set(gcf,'defaultaxesfontsize',15);

    subplot(3,1,1);
    [X,Y] = ndgrid(lags/12,Z_YvPar);
    contourf(X,Y,ZvPar_lr',[-1e50 Zcxs(1):Zcxs(2):Zcxs(3) 1e50],'linestyle','none');
    caxis([Zcxs(1) Zcxs(3)]);
% $$$     xlabel('Lag (years)');
    set(gca,'xticklabel',[]);
    ylabel(zlab);
    ylim(zlim);
    text(-1.95,6.5,['(a) $\Theta_z(p_z)$ ' label ' regression'],'BackgroundColor','w');
    cb = colorbar;
    ylabel(cb,'Temperature Anomaly ($^\circ$C / $\sigma$)','Interpreter','latex');
    set(gca,'Position',[0.13 0.7 0.75 0.28]);
    if (~remap_to_T)
        set(gca,'ydir','reverse');
        ax1=gca;
        ax1_pos = ax1.Position; % position of first axes
        ax1_pos(1) = ax1_pos(1)*0.5;
        ax1_pos(3) = 0.00001;
        ax2 = axes('Position',ax1_pos,...
                   'Color','none');
        set(ax2,'FontSize',15);
        Zticks = [0:500:4000 5000];
% $$$         Zticks = [0:200:2000];
% $$$         Zticks = [0:100:800];
        Zticks = [0:50:250];
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
    [X,Y] = ndgrid(lags/12,T_YvPar);
    contourf(X,Y,TvPar_lr',[-1e50 Tcxs(1):Tcxs(2):Tcxs(3) 1e50],'linestyle','none');
    caxis([Tcxs(1) Tcxs(3)]);
% $$$     xlabel('Lag (years)');
    set(gca,'xticklabel',[]);
    ylabel(tlab);
    text(-1.95,6.5,['(b) $\Theta_\Theta(p_\Theta)$ ' label ' regression']);
    cb = colorbar;
    ylabel(cb,'Temperature Anomaly ($^\circ$C / $\sigma$)','Interpreter','latex');
    colormap(redblue);
    ylim(tlim);
    set(gca,'Position',[0.13 0.38 0.75 0.28]);
    if (~remap_to_T)
        set(gca,'ydir','reverse');
        ax1=gca;
        ax1_pos = ax1.Position; % position of first axes
        ax1_pos(1) = ax1_pos(1)*0.5;
        ax1_pos(3) = 0.00001;
        ax2 = axes('Position',ax1_pos,...
                   'Color','none');
        set(ax2,'FontSize',15);
        Zticks = [30 25 20 18 16 14];
% $$$         Zticks = [30 26 22 18:-2:10];
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

    subplot(3,1,3);
    [X,Y] = ndgrid(lags/12,y_YvPar);
    contourf(X,Y,yvar_lr',[-1e50 Ycxs(1):Ycxs(2):Ycxs(3) 1e50],'linestyle','none');
    caxis([Ycxs(1) Ycxs(3)]);
    xlabel('Lag (years)');
    ylabel(ylab);
    ylim(yylim);
    text(-1.95,45,['(b) $\Theta_\phi(p_\phi)$ ' label ' regression']);
    cb = colorbar;
    ylabel(cb,'Temperature Anomaly ($^\circ$C / $\sigma$)','Interpreter','latex');
    set(gca,'Position',[0.13 0.07 0.75 0.28]);
    if (~remap_to_T)
        ax1=gca;
        ax1_pos = ax1.Position; % position of first axes
        ax1_pos(1) = ax1_pos(1)*0.5;
        ax1_pos(3) = 0.00001;
        ax2 = axes('Position',ax1_pos,...
                   'Color','none');
        set(ax2,'FontSize',15);
        Zticks = [-75:15:75];
% $$$         Zticks = [-70:10:80];
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
