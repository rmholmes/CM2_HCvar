% Plotting scripts for Holmes, Sohail and Zika J Climate J-CLI-D-21-0695

% Load .mat data created by PostProcess_ACCESSCM2.m
clear all;
baseMAT = 'C:/Users/rhol9417/data/access-cm2/';
load([baseMAT 'PIcontrolTb05PP_Tint.mat']);

% For proper text labels:
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

%%%% Plot mean profiles and conversions (Figure 2):

figure;
set(gcf,'Position',[118         316        1673         635]);
subplot(1,4,1);
plot(Te,mean(Tv.P,2),'linewidth',2);
set(gca,'ydir','reverse');
ylabel('Time-mean Temperature Percentile ($\overline{p_\Theta}(\Theta)$)');
xlabel('Temperature ($^\circ$C)');
xlim([-3 34]);
ylim([0 100]);
grid on;
text(-2,3,'(a)','FontSize',15);
set(gca,'Position',[0.0568    0.1006    0.1809    0.8150]);

subplot(1,4,2);
plot(Ze,mean(Zv.P,2),'linewidth',2);
set(gca,'ydir','reverse');
ylabel('Time-mean Depth Percentile ($\overline{p_z}(z)$)');
xlabel('Depth (m)');
xlim([0 5000]);
ylim([0 100]);
grid on;
text(4300,3,'(b)','FontSize',15);

subplot(1,4,3);
plot(latv_edges,mean(Yv.P,2),'linewidth',2);
set(gca,'ydir','reverse');
ylabel('Time-mean Latitude Percentile ($\overline{p_\phi}(\phi)$)');
xlabel('Latitude ($^\circ$N)');
xlim([-90 90]);
ylim([0 100]);
grid on;
text(65,3,'(c)','FontSize',15);

subplot(1,4,4);
plot(TvP.Tp_mean,P,'-r','linewidth',2);
hold on;
plot(ZvP.Tp_mean,P,'-k','linewidth',2);
plot(YvP.Tp_mean,P,'-b','linewidth',2);
set(gca,'ydir','reverse');
legend('$\overline{\Theta_\Theta}(p)$','$\overline{\Theta_z}(p)$', ...
       '$\overline{\Theta_\phi}(p)$','Location','southeast');
ylabel('Percentile');
xlabel('Temperature ($^\circ$C)');
grid on;
text(2,3,'(d)','FontSize',15);

%%%%%% All anomalies P-t (Figure 4):

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
    zlab = 'Depth Percentile $p$';

    T_YvPar = P;
    tlim = zlim;
    tlab = 'Temperature Percentile $p_\Theta$';
    
    y_YvPar = P;
    yylim = zlim;
    ylab = 'Latitude Percentile $p_\phi$';
end        

ZvPar = ZvP.Tp;TvPar = TvP.Tp;yvar = YvP.Tp;
Zcxs = [-0.06 0.002 0.06]; % PI control anomalies
Zlab = '$\Theta_z(p)$';
Tcxs = Zcxs;
Tlab = '$\Theta_\Theta(p_\Theta)$';
ycxs = Zcxs;
Ylab = '$\Theta_\phi(p_\phi)$';

zlim = [0 100];
tlim = [0 100];
yylim = [0 100];
xlims = [yr1 yr1+nyrs];

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
cb = colorbar;
ylabel(cb,'Temperature Anomalies ($^\circ$C)','Interpreter','latex');
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
colormap(hot);
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
cb = colorbar;
ylabel(cb,'Temperature Anomalies ($^\circ$C)','Interpreter','latex');
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

EqP = interp1(yofP_mean,P,0,'linear');
hold on;
plot([50 600],[EqP EqP],'--k');

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
axes1i = axes('Position',[0.6  0.738    0.2    0.12]);
axes(axes1i);
contourf(X,Y,v,[-1e50 Zcxs(1):Zcxs(2):Zcxs(3) 1e50],'linestyle','none');
caxis([Zcxs(1) Zcxs(3)]);
ylim(zlim);
xlim(xlims);
set(gca,'xtick',[xlims(1) mean(xlims) xlims(2)]);
cb = colorbar;
set(gca,'ydir','reverse');
set(gca,'Position',[0.63    0.738    0.2    0.12]);

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

% $$$ axes1iN = axes('Position',[0.655    0.87    0.2    0.07]);
% $$$ axes(axes1iN);

y = T_YvPar;
v = TvPar(:,t1:t2)';
[X,Y] = ndgrid(x,y);
axes2i = axes('Position',[0.63    0.4448    0.2    0.12]);
axes(axes2i);
contourf(X,Y,v,[-1e50 Tcxs(1):Tcxs(2):Tcxs(3) 1e50],'linestyle','none');
caxis([Tcxs(1) Tcxs(3)]);
ylim(tlim);
xlim(xlims);
set(gca,'xtick',[xlims(1) mean(xlims) xlims(2)]);
cb = colorbar;
set(gca,'ydir','reverse');
set(gca,'Position',[0.63    0.4448    0.2    0.12]);
if (~remap_to_T)
    set(gca,'ydir','reverse');
    ax1=gca;
    ax1_pos = ax1.Position; % position of first axes
    ax1_pos(1) = ax1_pos(1)-0.02;
    ax1_pos(3) = 0.00001;
    ax2 = axes('Position',ax1_pos,...
               'Color','none');
    set(ax2,'FontSize',15);
    Zticks = [28 20 14 12];
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
axes3i = axes('Position',[0.63    0.135    0.2    0.12]);
axes(axes3i);
contourf(X,Y,v,[-1e50 ycxs(1):ycxs(2):ycxs(3) 1e50],'linestyle','none');

caxis([ycxs(1) ycxs(3)]);
ylim(yylim);
xlim(xlims);
set(gca,'xtick',[xlims(1) mean(xlims) xlims(2)]);
cb = colorbar;

EqP = interp1(yofP_mean,P,0,'linear');
hold on;
plot(xlims,[EqP EqP],'--k');

set(gca,'Position',[0.63    0.135    0.2    0.12]);
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

colormap(hot);

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

%%% Temperature anomalies standard deviation (Figure 5):
figure;
set(gcf,'Position',[23         208        1896         716]);
set(gcf,'defaulttextfontsize',15);
set(gcf,'defaultaxesfontsize',15);

pos1 = [0.1759    0.1030    0.18    0.8150];
axs = axes('Position',pos1);
plot(std(ZvP.Tp,[],2),P,'-k','linewidth',2);
hold on;
plot(std(TvP.Tp,[],2),P,'-r','linewidth',2);
plot(std(YvP.Tp,[],2),P,'-b','linewidth',2);
xlabel('Standard Deviation of $\Theta$ ($^\circ$C)');
ylim([0 100]);
xlim([0 0.1]);
set(gca,'ydir','reverse');

% Averages:
text(0.04,30,['$\overline{\sigma_z}$ = ' sprintf('%3.3f',mean(std(ZvP.Tp,[],2))) '$^\circ$C'],'color','k');
text(0.04,35,['$\overline{\sigma_\Theta}$ = ' sprintf('%3.3f',mean(std(TvP.Tp,[],2))) '$^\circ$C'],'color','r');
text(0.04,40,['$\overline{\sigma_\phi}$ = ' sprintf('%3.3f',mean(std(YvP.Tp,[],2))) '$^\circ$C'],'color','b');

ax2_pos = pos1;
ax2_pos(1) = pos1(1)-pos1(1)/4;
ax2_pos(3) = 0.00001;
ax3_pos = ax2_pos;
ax3_pos(1) = pos1(1)-1.9*pos1(1)/4;
ax4_pos = ax2_pos;
ax4_pos(1) = pos1(1)-3.1*pos1(1)/4;

ax2 = axes('Position',ax2_pos,'XColor','r','YColor','r','FontSize',15);
set(ax2,'FontSize',15);
Zticks = [28 18 14 12:-2:4 3:-1:-1];
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

%%%% ENSO Lag Regressions (Figures 6-7):
   
    % ENSO:
    ts = CIN.N34-mean(CIN.N34);
    ts = ts/std(ts);
    label = 'Ni\~{n}o 3.4';
    lags = [-12*2:1:12*2]; %months lag
    zlim = [0 7]; tlim = [0 7]; yylim = [40 85];

% Variables:
    Z_YvPar = P;    zlab = 'Depth Percentile $p$';
    T_YvPar = P;    tlab = 'Temperature Percentile $p$';
    y_YvPar = P;    ylab = 'Latitude Percentile $p$';

    vars = {'Tp','Tp','TENp','ADVp','FORp','VMIXp','RMIXp'};
    tder = [0 1 1 1 1 1 1 1 1]; % Take time-derivative of budget terms first

    remap_to_T = 0;
    ll = length(lags);
    reg_or_corr = 1; % Regression (1) or correlation (0)

    % Calculate regressions:
    for ti=1:length(typs)
        for vi=1:length(vars)
            vars{vi}
            eval(['var = ' typs{ti} 'vP.' vars{vi} ';']);
            if (tder(vi))
                var = diff(var,[],2)./repmat(DT_A(2:end)',[length(var(:,1)) 1]);
                var = cat(2,zeros(length(var(:,1)),1),var);
            end
            lg = length(var(:,1));
            lr = zeros(lg,ll);
            
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
                if (reg_or_corr)
                    lr(:,ii) = var*TS/(sum(TS.^2));
                else
                    varnor = (var-repmat(mean(var,2),[1 tL]))./ ...
                            repmat(std(var,[],2),[1 tL]);
                    lr(:,ii) = var*TS/sqrt(sum(TS.^2))./ ...
                        sqrt(sum(var.^2,2));
                end
            end
            if (tder(vi))
                eval([typs{ti} 'vP.' vars{vi} 'd_lr = lr;']);
            else
                eval([typs{ti} 'vP.' vars{vi} '_lr = lr;']);
            end            
        end
    end
    
    % Plot temperatures:
    Zcxs = [-0.08 0.0025 0.08]; Tcxs = Zcxs;Ycxs=[-0.02 0.001 0.02];

    figure;
    set(gcf,'Position',[3          40        1130         963]);
    set(gcf,'defaulttextfontsize',15);
    set(gcf,'defaultaxesfontsize',15);

    subplot(3,1,1);
    [X,Y] = ndgrid(lags/12,Z_YvPar);
    contourf(X,Y,ZvP.Tp_lr',[-1e50 Zcxs(1):Zcxs(2):Zcxs(3) 1e50],'linestyle','none');
    hold on;
    contour(X,Y,ZvP.Tp_lr',[0.02:0.02:0.08],'-k');
    contour(X,Y,ZvP.Tp_lr',[-0.08:0.02:-0.02],'--k');
    caxis([Zcxs(1) Zcxs(3)]);
% $$$     xlabel('Lag (years)');
    set(gca,'xticklabel',[]);
    ylabel(zlab);
    ylim(zlim);
    text(-1.95,6.5,['(a) $\Theta_z(p)$ ' label ' regression'],'BackgroundColor','w');
    cb = colorbar;
    set(cb,'ytick',[-0.08:0.04:0.08]);
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
    contourf(X,Y,TvP.Tp_lr',[-1e50 Tcxs(1):Tcxs(2):Tcxs(3) 1e50],'linestyle','none');
    hold on;
    contour(X,Y,TvP.Tp_lr',[0.02:0.02:0.08],'-k');
    contour(X,Y,TvP.Tp_lr',[-0.08:0.02:-0.02],'--k');
    caxis([Tcxs(1) Tcxs(3)]);
% $$$     xlabel('Lag (years)');
    set(gca,'xticklabel',[]);
    ylabel(tlab);
    text(-1.95,6.5,['(b) $\Theta_\Theta(p)$ ' label ' regression']);
    cb = colorbar;
    set(cb,'ytick',[-0.08:0.04:0.08]);
    ylabel(cb,'Temperature Anomaly ($^\circ$C / $\sigma$)','Interpreter','latex');
    colormap(hot);
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
        Zticks = [25 20 18 16 14];
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
    contourf(X,Y,YvP.Tp_lr',[-1e50 Ycxs(1):Ycxs(2):Ycxs(3) 1e50],'linestyle','none');
    hold on;
    contour(X,Y,YvP.Tp_lr',[0.02:0.02:0.08],'-k');
    contour(X,Y,YvP.Tp_lr',[-0.08:0.02:-0.02],'--k');
    caxis([Ycxs(1) Ycxs(3)]);
    xlabel('Lag (years)');
    ylabel(ylab);
    ylim(yylim);
    text(-1.95,45,['(b) $\Theta_\phi(p)$ ' label ' regression']);
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
    
    % Plot budget terms:
    vars = {'Tp','TENp','ADVp','FORp','VMIXp','RMIXp'};
    names = {'d$\Theta$/dt','Tendency','Numerical Mixing','Forcing','Vertical Mixing','Neutral Mixing'};
    names = {'d$\Theta$/dt','Tendency','Advection','Forcing','Vertical Mixing','Neutral Mixing'};

% $$$     vars = {'FORp','VMIXp','ADVp'};
% $$$     names = {'Surface Forcing','Vertical Mixing','Advection'};
% $$$     names = {'Surface Forcing','Vertical Mixing','Numerical Mixing'};
    txtlabs = {'(a)','(b)','(c)','(d)','(e)','(f)'};
    
    figure;
    set(gcf,'Position',[1          36        1920         970]);
    set(gcf,'defaulttextfontsize',15);
    set(gcf,'defaultaxesfontsize',15);

    Zcxs = [-5e-9 1e-10 5e-9];
% $$$     Zcxs = [-1e-10 0.5e-11 1e-10];
% $$$     Zcxs = [-0.25e-10 0.1e-11 0.25e-10];
    [X,Y] = ndgrid(lags/12,Z_YvPar);
    for vi=1:length(vars)
        eval(['var = ZvP.' vars{vi} 'd_lr;']);
        subplot(3,2,vi);
% $$$         subplot(3,2,2*vi-1);
        contourf(X,Y,var',[-1e50 Zcxs(1):Zcxs(2):Zcxs(3) 1e50],'linestyle','none');
        caxis([Zcxs(1) Zcxs(3)]);
        xlabel('Lag (years)');
% $$$         set(gca,'xticklabel',[]);
        ylabel(zlab);
        ylim(zlim);
% $$$         text(-1.95,6.5,[names{vi} ' ' label ' regression'],'BackgroundColor','w');
% $$$         text(-1.95,6.5,[txtlabs{vi} ' ' names{vi}],'BackgroundColor','w');
        text(-49,6.5,[txtlabs{vi} ' ' names{vi}]);
        cb = colorbar;
        ylabel(cb,'$\Theta_z$ Tendency ($^\circ$Cs$^{-1}$ / $\sigma$)','Interpreter','latex');
        set(gca,'ydir','reverse');
    end
        colormap(hot);

    figure;
    set(gcf,'Position',[1          36        1920         970]);
    set(gcf,'defaulttextfontsize',15);
    set(gcf,'defaultaxesfontsize',15);

    Tcxs = [-3e-9 0.2e-10 3e-9];
% $$$     Tcxs = [-0.25e-10 0.1e-11 0.25e-10];
    [X,Y] = ndgrid(lags/12,T_YvPar);
    for vi=1:length(vars)
        eval(['var = TvP.' vars{vi} 'd_lr;']);
        subplot(3,2,vi);
% $$$         subplot(3,2,2*vi);
        contourf(X,Y,var',[-1e50 Tcxs(1):Tcxs(2):Tcxs(3) 1e50],'linestyle','none');
        caxis([Tcxs(1) Tcxs(3)]);
        xlabel('Lag (years)');
        ylabel(tlab);
% $$$         text(-1.95,6.5,[names{vi} ' ' label ' regression']);
% $$$         text(-1.95,6.5,[txtlabs{vi+3} ' ' names{vi}]);
        text(-49,6.5,[txtlabs{vi} ' ' names{vi}]);
        cb = colorbar;
        ylabel(cb,'$\Theta_\Theta$ Tendency ($^\circ$Cs$^{-1}$ / $\sigma$)','Interpreter','latex');
        colormap(hot);
        ylim(tlim);
        set(gca,'ydir','reverse');
    end
        colormap(hot);

    figure;
    set(gcf,'Position',[1          36        1920         970]);
    set(gcf,'defaulttextfontsize',15);
    set(gcf,'defaultaxesfontsize',15);

    Ycxs = [-1e-9 1e-11 1e-9];
    [X,Y] = ndgrid(lags/12,y_YvPar);
    for vi=1:length(vars)
        eval(['var = YvP.' vars{vi} 'd_lr;']);
        subplot(3,2,vi);
        contourf(X,Y,var',[-1e50 Ycxs(1):Ycxs(2):Ycxs(3) 1e50],'linestyle','none');
        caxis([Ycxs(1) Ycxs(3)]);
        xlabel('Lag (years)');
        ylabel(ylab);
        ylim(yylim);
% $$$         text(-1.95,45,[names{vi} ' ' label ' regression']);
        text(-1.95,45,names{vi});
        cb = colorbar;
        ylabel(cb,'$\Theta_\phi$ Tendency ($^\circ$Cs$^{-1}$ / $\sigma$)','Interpreter','latex');
    end
        colormap(hot);


%%% Heat content anomalies and correlations (Figure 12):
axes(axs); cla;
nfilt = 12*5-1;
Zvar = ZvP.Hp;
Tvar = TvP.Hp;
Yvar = YvP.Hp;

ZvarLP = filter_field(Zvar,nfilt,'-t');
TvarLP = filter_field(Tvar,nfilt,'-t');
YvarLP = filter_field(Yvar,nfilt,'-t');
ZvarLP(isnan(ZvarLP)) = 0;
TvarLP(isnan(TvarLP)) = 0;
YvarLP(isnan(YvarLP)) = 0;

ZvarHP = Zvar - filter_field(Zvar,nfilt,'-t');
TvarHP = Tvar - filter_field(Tvar,nfilt,'-t');
YvarHP = Yvar - filter_field(Yvar,nfilt,'-t');
ZvarHP(isnan(ZvarHP)) = 0;
TvarHP(isnan(TvarHP)) = 0;
YvarHP(isnan(YvarHP)) = 0;

plot(std(Zvar/1e22,[],2),P,'-k','linewidth',2);
hold on;
plot(std(Tvar/1e22,[],2),P,'-r','linewidth',2);
plot(std(Yvar/1e22,[],2),P,'-b','linewidth',2);

plot(std(ZvarLP/1e22,[],2),P,'--k','linewidth',1);
plot(std(TvarLP/1e22,[],2),P,'--r','linewidth',1);
plot(std(YvarLP/1e22,[],2),P,'--b','linewidth',1);

plot(std(ZvarHP/1e22,[],2),P,':k','linewidth',1);
plot(std(TvarHP/1e22,[],2),P,':r','linewidth',1);
plot(std(YvarHP/1e22,[],2),P,':b','linewidth',1);

xlabel('Standard Deviation of $\mathcal{H}$ ($10^{22}$ J)');
ylim([0 100]);
xlim([0 1.6]);
set(gca,'ydir','reverse');
% $$$ set(gca,'yticklabel',[]);
ylabel('Percentile');
text(0.001,2,'(a)');

text(0.1,30,['$\overline{\sigma_z}$ = ' sprintf('%3.2f',mean(std(ZvP.Hp,[],2))/1e22) ' $\times10^{22}J$'],'color','k');
text(0.1,35,['$\overline{\sigma_\Theta}$ = ' sprintf('%3.2f',mean(std(TvP.Hp,[],2))/1e22) ' $\times10^{22}J$'],'color','r');
text(0.1,40,['$\overline{\sigma_\phi}$ = ' sprintf('%3.2f',mean(std(YvP.Hp,[],2))/1e22) ' $\times10^{22}J$'],'color','b');

% Full -> full

TS = (OHC-mean(OHC))/std(OHC);

ZvPar = ZvP.Hp;TvPar = TvP.Hp;yvar = YvP.Hp;
Z_YvPar = P;    zlab = 'Depth Percentile $p$';
T_YvPar = P;    tlab = 'Temperature Percentile $p_\Theta$';
y_YvPar = P;    ylab = 'Latitude Percentile $p_\phi$';
    
Znor = (ZvPar-repmat(mean(ZvPar,2),[1 tL]))./ ...
       repmat(std(ZvPar,[],2),[1 tL]);
Tnor = (TvPar-repmat(mean(TvPar,2),[1 tL]))./ ...
       repmat(std(TvPar,[],2),[1 tL]);
Ynor = (yvar-repmat(mean(yvar,2),[1 tL]))./ ...
       repmat(std(yvar,[],2),[1 tL]);

% $$$     % Apply low-pass filter:
% $$$     nfilt = 12*9-1;
% $$$     Znor = filter_field(Znor,nfilt,'-t');
% $$$     Tnor = filter_field(Tnor,nfilt,'-t');
% $$$     Ynor = filter_field(Ynor,nfilt,'-t');
% $$$     TS = filter_field(TS,nfilt,'-t');
% $$$     Znor(isnan(Znor)) = 0;
% $$$     Tnor(isnan(Tnor)) = 0;
% $$$     Ynor(isnan(Ynor)) = 0;
% $$$     TS(isnan(TS)) = 0;
% $$$ 
% $$$     % Apply high-pass filter:
% $$$     nfilt = 12*9-1;
% $$$     Znor = Znor - filter_field(Znor,nfilt,'-t');
% $$$     Tnor = Tnor - filter_field(Tnor,nfilt,'-t');
% $$$     Ynor = Ynor - filter_field(Ynor,nfilt,'-t');
% $$$     TS = TS - filter_field(TS,nfilt,'-t');
% $$$     Znor(isnan(Znor)) = 0;
% $$$     Tnor(isnan(Tnor)) = 0;
% $$$     Ynor(isnan(Ynor)) = 0;
% $$$     TS(isnan(TS)) = 0;
    

Znor = (Znor-repmat(mean(Znor,2),[1 tL]))./ ...
       repmat(std(Znor,[],2),[1 tL]);
Tnor = (Tnor-repmat(mean(Tnor,2),[1 tL]))./ ...
       repmat(std(Tnor,[],2),[1 tL]);
Ynor = (Ynor-repmat(mean(Ynor,2),[1 tL]))./ ...
       repmat(std(Ynor,[],2),[1 tL]);
TS = (TS-mean(TS))/std(TS);
    
ZvPar_lr = Znor*TS/sqrt(sum(TS.^2))./ ...
    sqrt(sum(Znor.^2,2));
TvPar_lr = Tnor*TS/sqrt(sum(TS.^2))./ ...
    sqrt(sum(Tnor.^2,2));
yvar_lr = Ynor*TS/sqrt(sum(TS.^2))./ ...
          sqrt(sum(Ynor.^2,2));

% +-60 -> full
baseMAT = 'D:/DATA/access-cm2/';
PM60data = load([baseMAT 'PIcontrolTb05PP_ypm60_Tint.mat'],'ZvP','TvP','YvP');
    
Znor60 = (PM60data.ZvP.Hp-repmat(mean(PM60data.ZvP.Hp,2),[1 tL]))./ ...
       repmat(std(PM60data.ZvP.Hp,[],2),[1 tL]);
Tnor60 = (PM60data.TvP.Hp-repmat(mean(PM60data.TvP.Hp,2),[1 tL]))./ ...
       repmat(std(PM60data.TvP.Hp,[],2),[1 tL]);
Ynor60 = (PM60data.YvP.Hp-repmat(mean(PM60data.YvP.Hp,2),[1 tL]))./ ...
       repmat(std(PM60data.YvP.Hp,[],2),[1 tL]);

% $$$     % Apply low-pass filter:
% $$$     nfilt = 12*9-1;
% $$$     Znor = filter_field(Znor,nfilt,'-t');
% $$$     Tnor = filter_field(Tnor,nfilt,'-t');
% $$$     Ynor = filter_field(Ynor,nfilt,'-t');
% $$$     TS = filter_field(TS,nfilt,'-t');
% $$$     Znor(isnan(Znor)) = 0;
% $$$     Tnor(isnan(Tnor)) = 0;
% $$$     Ynor(isnan(Ynor)) = 0;
% $$$     TS(isnan(TS)) = 0;
% $$$ 
% $$$     % Apply high-pass filter:
% $$$     nfilt = 12*9-1;
% $$$     Znor = Znor - filter_field(Znor,nfilt,'-t');
% $$$     Tnor = Tnor - filter_field(Tnor,nfilt,'-t');
% $$$     Ynor = Ynor - filter_field(Ynor,nfilt,'-t');
% $$$     TS = TS - filter_field(TS,nfilt,'-t');
% $$$     Znor(isnan(Znor)) = 0;
% $$$     Tnor(isnan(Tnor)) = 0;
% $$$     Ynor(isnan(Ynor)) = 0;
% $$$     TS(isnan(TS)) = 0;
    

Znor60 = (Znor60-repmat(mean(Znor60,2),[1 tL]))./ ...
       repmat(std(Znor60,[],2),[1 tL]);
Tnor60 = (Tnor60-repmat(mean(Tnor60,2),[1 tL]))./ ...
       repmat(std(Tnor60,[],2),[1 tL]);
Ynor60 = (Ynor60-repmat(mean(Ynor60,2),[1 tL]))./ ...
       repmat(std(Ynor60,[],2),[1 tL]);
TS = (TS-mean(TS))/std(TS);
    
ZvPar_lr60 = Znor60*TS/sqrt(sum(TS.^2))./ ...
    sqrt(sum(Znor60.^2,2));
TvPar_lr60 = Tnor60*TS/sqrt(sum(TS.^2))./ ...
    sqrt(sum(Tnor60.^2,2));
yvar_lr60 = Ynor60*TS/sqrt(sum(TS.^2))./ ...
          sqrt(sum(Ynor60.^2,2));

pos1 = [0.3808    0.1030    0.18    0.8150];
axs = axes('Position',pos1);%[0.8    0.1030    0.18    0.8150]);
plot(ZvPar_lr,P,'-k','linewidth',2);
hold on;
plot(TvPar_lr,P,'-r','linewidth',2);
plot(yvar_lr,P,'-b','linewidth',2);
plot(ZvPar_lr60,P,'--k','linewidth',1);
plot(TvPar_lr60,P,'--r','linewidth',1);
plot(yvar_lr60,P,'--b','linewidth',1);
xlabel('Correlation with total OHC');
xlim([-0.45 1]);
set(gca,'yticklabel',[]);
ylim([0 100]);
legend('$H_z(p)$','$H_\Theta(p_\Theta)$','$H_\phi(p_\phi)$','Location','SouthWest');
set(gca,'ydir','reverse')
text(0.001,2,'(b)');
% $$$ text(0.88,2,'(c)');

% Linear trends:
%----------------

% $$$ % Testing on total OHC:
% $$$ t1 = 1;
% $$$ t2 = tL;
% $$$ timean = mean(reshape(time(t1:t2),[12 (t2-t1+1)/12]),1)';
% $$$ OHCan = squeeze(mean(reshape(OHC(t1:t2),[12 (t2-t1+1)/12]),1));
% $$$ 
% $$$ windows = [3:4:100];
% $$$ trends_std = zeros(size(windows));
% $$$ for wi = 1:length(windows)
% $$$     trends = lintrends(timean,OHCan',windows(wi));
% $$$     trends_std(wi) = nanstd(trends);
% $$$ end
% $$$ plot(windows,trends_std);
% $$$ % Checks out. 

% Calculate annual average:
t1 = 1;
t2 = tL;
timean = mean(reshape(time(t1:t2),[12 (t2-t1+1)/12]),1)';
Zvaran = squeeze(mean(reshape(ZvP.Hp(:,t1:t2)',[12 (t2-t1+1)/12 PL]),1));
Tvaran = squeeze(mean(reshape(TvP.Hp(:,t1:t2)',[12 (t2-t1+1)/12 PL]),1));
Yvaran = squeeze(mean(reshape(YvP.Hp(:,t1:t2)',[12 (t2-t1+1)/12 PL]),1));

% Calculate trends over specified windows:
windows = [5 15 51];

Ztrends = zeros(length(windows),PL);
Ttrends = zeros(length(windows),PL);
Ytrends = zeros(length(windows),PL);

for wi = 1:length(windows)
    window = windows(wi)
    for pi=1:PL
        trends = lintrends(timean,Zvaran(:,pi),window);
        Ztrends(wi,pi) = nanstd(trends);
        trends = lintrends(timean,Tvaran(:,pi),window);
        Ttrends(wi,pi) = nanstd(trends);
        trends = lintrends(timean,Yvaran(:,pi),window);
        Ytrends(wi,pi) = nanstd(trends);
    end
end

%%% Plot them:

figure;
set(gcf,'Position',[23         208        1896         716]); % total OHC correlation.
set(gcf,'defaulttextfontsize',15);
set(gcf,'defaultaxesfontsize',15);

poss = [0.1759    0.1030    0.18    0.8150; ...
        0.3808    0.1030    0.18    0.8150; ...
        0.586    0.1030    0.18    0.8150];


pos1 = poss(1,:);
axs = axes('Position',pos1);
plot(Ztrends(1,:)/1e20,P,'-k','linewidth',2);
hold on;
plot(Ttrends(1,:)/1e20,P,'-r','linewidth',2);
plot(Ytrends(1,:)/1e20,P,'-b','linewidth',2);
xlabel(['' num2str(windows(1)) '-year trends (J/year)']);
ylim([0 100]);
% $$$ xlim([0 0.1]);
set(gca,'ydir','reverse');

% Add extra axes:
ax2_pos = pos1;
ax2_pos(1) = pos1(1)-pos1(1)/4;
ax2_pos(3) = 0.00001;
ax3_pos = ax2_pos;
ax3_pos(1) = pos1(1)-1.9*pos1(1)/4;
ax4_pos = ax2_pos;
ax4_pos(1) = pos1(1)-3.1*pos1(1)/4;

ax2 = axes('Position',ax2_pos,'XColor','r','YColor','r','FontSize',15);
set(ax2,'FontSize',15);
Zticks = [28 18 14 12:-2:4 3:-1:-1];
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

pos1 = poss(2,:);
axs = axes('Position',pos1);
plot(Ztrends(2,:)/1e20,P,'-k','linewidth',2);
hold on;
plot(Ttrends(2,:)/1e20,P,'-r','linewidth',2);
plot(Ytrends(2,:)/1e20,P,'-b','linewidth',2);
xlabel(['' num2str(windows(2)) '-year trends (J/year)']);
ylim([0 100]);
% $$$ xlim([0 0.1]);
set(gca,'ydir','reverse');

pos1 = poss(3,:);
axs = axes('Position',pos1);
plot(Ztrends(3,:)/1e20,P,'-k','linewidth',2);
hold on;
plot(Ttrends(3,:)/1e20,P,'-r','linewidth',2);
plot(Ytrends(3,:)/1e20,P,'-b','linewidth',2);
xlabel(['' num2str(windows(3)) '-year trends (J/year)']);
ylim([0 100]);
% $$$ xlim([0 0.1]);
set(gca,'ydir','reverse');

%%% Tyz and MOCyz plotting (Figures 8-10):

% Load Tyz and MOCyz data:
baseMAT = 'D:/DATA/access-cm2/';
load([baseMAT 'PIcontrolTb05PP_Tint_MOC.mat']);
load([baseMAT 'PIcontrolTb05PP_Tint_Tyz.mat']);

[tmp ind1] = min(abs(P-10));[tmp ind2] = min(abs(P-30));
ts = mean(ZvP.Tp(ind1:ind2,:),1)';
ts = filter_field(ts,12*10+1,'-t');
ts(isnan(ts)) = 0;
label = '$\Theta_z(10<p<30,t';
name = 'Tzp10to30_MOCfull_TyzReg_10yrSmoothing';

[tmp ind1] = min(abs(P-90));[tmp ind2] = min(abs(P-100));
ts = mean(TvP.Tp(ind1:ind2,:),1)';
ts = filter_field(ts,12*1+1,'-t'); % Filtering
ts(isnan(ts)) = 0;
label = '$\Theta_\Theta(90<p_\Theta<100,t';
name = 'TTp90to100_MOCfull_TyzReg_10yrSmoothing';

[tmp ind1] = min(abs(P-95));[tmp ind2] = min(abs(P-100));
ts = mean(YvP.Tp(ind1:ind2,:),1)';
ts = filter_field(ts,1*12+1,'-t');
ts(isnan(ts)) = 0;
label = '$\Theta_\phi(95<p<100,t';
name = 'Tzp90to100_MOCfull_TyzReg';

% $$$ lags = [-100:10:100];
% $$$ lags = [0]
lags = [-8:1:8];

for li = 1:length(lags)
    lag = lags(li)*12;

% One time series regressions:
figure;
set(gcf,'defaulttextfontsize',15);
set(gcf,'defaultaxesfontsize',15);
set(gcf,'Position',[963    43   956   953]);

ts_lagged = zeros(size(ts));
if (lag<0)
    ts_lagged(1:(end+lag)) = ts((-lag+1):end);
    ts_lagged((end+lag+1):end) = 0;
elseif (lag == 0)
    ts_lagged = ts;
else
    ts_lagged(1:(lag)) = 0;
    ts_lagged((lag+1):end) = ts(1:(end-lag));
end

ts_norm = (ts_lagged-nanmean(ts_lagged))/std(ts_lagged);

MOCfull = MOC+MOCgm;
MOCyz_reg = reshape(reshape(MOCfull,[yL*zL tL])*ts_norm/(sum(ts_norm.^2)),[yL zL]);
Tyz_reg = reshape(reshape(Tyz,[yL*zL tL])*ts_norm/(sum(ts_norm.^2)),[yL zL]);

subplot(3,2,[1 2]);
plot(time,ts_lagged,'-k');
hold on;
plot(time,ts,'--k');
legend([label '-' num2str(lag/12) ')$'],[label ')$']);
% $$$ legend([label ')$']);
xlim([50 600]);
ylim([-0.03 0.03]);
xlabel('Year');
ylabel('Temperature $(^\circ$C)');
grid on;
text(53,0.026,'(a)');
set(gca,'Position',[0.0973 0.7135 0.7464 0.2157]);

[X,Y] = ndgrid(latv,Z);

subplot(3,2,[5 6]);
contourf(X,Y,Tyz_reg,[-1e50 -0.5e20:0.25e19:0.5e20 1e50],'linestyle','none');
hold on;
[c,h] = contour(X,Y,Tyz_mean,[-2:2:40],'-k');
clabel(c,h);
caxis([-0.5 0.5]*1e20);
cb = colorbar;
ylabel(cb,'J/m/$^\circ$-latitude/$\sigma$','Interpreter','latex');
set(gca,'ydir','reverse');
% $$$ xlabel('Latitude ($^\circ$N)');
set(gca,'xticklabel',[]);
ylabel('Depth (m)');
ylim([0 4000]);
xlim([-80 70]);
colormap('hot');
text(-79,200,'(b)');
set(gca,'Position',[0.0973 0.3788 0.7464 0.2739]);

subplot(3,2,[5 6]);
% $$$ contourf(X,Y,MOCyz_reg,[-1e50 -0.5:0.025:0.5 1e50],'linestyle','none');
contourf(X,Y,MOCyz_reg,[-1e50 -1:0.05:1 1e50],'linestyle','none');
hold on;
[c,h] = contour(X,Y,MOC_mean+MOCgm_mean,[3:3:100],'-k');
clabel(c,h);
[c,h] = contour(X,Y,MOC_mean+MOCgm_mean,[-99:3:-3],'--k');
clabel(c,h);
caxis([-1 1]);
cb = colorbar;
ylabel(cb,'Sv/$\sigma$','Interpreter','latex');
set(gca,'ydir','reverse');
xlabel('Latitude ($^\circ$N)');
% $$$ set(gca,'xticklabel',[]);
ylabel('Depth (m)');
ylim([0 4000]);
xlim([-80 70]);
text(-79,200,'(c)');
set(gca,'Position',[0.0973 0.0703 0.7464 0.2739]);

% $$$ if (lag>0)
% $$$     saveas(gcf,sprintf(['Decadal/' name '_%02d.png'],li));
% $$$ else
% $$$     saveas(gcf,['Decadal/ZTp10to30_MOC_TyzReg_10yrSmoothing_' ...
% $$$                 num2str(-lag/12) 'lag.png']);
% $$$ end    
end


%%%%% Spectral analysis (Figure 11): %%%%%%%%%%%%%%%%%%%%%%%%%

fs = 12;
nw = 10;
[pxx,f] = pmtm(ZvP.Tp(15,:)',nw,[],fs);
fL = length(f);
df = f(2)-f(1);

TorH = 0;

% Combination variables:
ZvP.FORpVMIXp = ZvP.FORp+ZvP.VMIXp;
TvP.FORpVMIXp = TvP.FORp+TvP.VMIXp;
YvP.FORpVMIXp = YvP.FORp+YvP.VMIXp;

ZvP.ADVpVMIXp = ZvP.ADVp+ZvP.VMIXp;
TvP.ADVpVMIXp = TvP.ADVp+TvP.VMIXp;
YvP.ADVpVMIXp = YvP.ADVp+YvP.VMIXp;

ZvP.MIXp = ZvP.RMIXp+ZvP.VMIXp;
TvP.MIXp = TvP.RMIXp+TvP.VMIXp+TvP.ADVp;
YvP.MIXp = YvP.RMIXp+YvP.VMIXp;

ZvP.FORpMIXp = ZvP.FORp+ZvP.MIXp;
TvP.FORpMIXp = TvP.FORp+TvP.MIXp;
YvP.FORpMIXp = YvP.FORp+YvP.MIXp;

ZvP.ADVpFORp = ZvP.FORp+ZvP.ADVp;
TvP.ADVpFORp = TvP.FORp+TvP.ADVp;
YvP.ADVpFORp = YvP.FORp+YvP.ADVp;

% $$$ ZvP.SUMp = ZvP.ADVpFORp+ZvP.MIXp;
% $$$ TvP.SUMp = TvP.FORp+TvP.MIXp;
% $$$ YvP.SUMp = YvP.ADVpFORp+YvP.MIXp;
% $$$ 
if (TorH)
% $$$ vars = {'Tp','Tp','TENp','ADVp','FORp','RMIXp','VMIXp','FORpVMIXp','ADVpVMIXp'};
% $$$ vars = {'Tp','Tp','TENp','ADVp','FORp','MIXp'};%,'VMIXp','FORpVMIXp','ADVpVMIXp'};
vars = {'Tp','Tp','TENp','ADVp','FORp','RMIXp','VMIXp'};%,'VMIXp','FORpVMIXp','ADVpVMIXp'};
% $$$ vars = {'Tp','Tp','TENp','ADVp','FORp','RMIXp','VMIXp','FORpVMIXp','ADVpVMIXp','MIXp','FORpMIXp','ADVpFORp'};%,'SUMp'};
tder = [0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; % Take time-derivative of budget terms first
else
vars = {'Hp','TEN_c','ADV_c','FOR_c','RMIX_c','VMIX_c'};%,'Hp','TEN_c','ADV_c','FOR_c','RMIX_c','VMIX_c'};
tder = [0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; % Take time-derivative of budget terms first
end

% $$$ vars = {'FORpVMIXp','ADVpVMIXp','MIXp','FORpMIXp','ADVpFORp'};
% $$$ tder = [1 1 1 1 1 1 1 1 1]; % Take time-derivative of budget terms first

tnames = {'$\Theta_\Theta(p_\Theta,t)$','$\Theta_z(p,t)$','$\Theta_\phi(p_\phi,t)$'};
tcols = {'r','k','b'};

names = {{'$\partial\Theta_\Theta/\partial t$','Tendency','Numerical Mixing','Forcing','Neutral Mixing','Vertical Mixing', ...
         'FOR+VMIX','ADV+VMIX','MIX','FOR+MIX','ADV+FOR'}, ...
         {'$\partial\Theta_z/\partial t$','Tendency','Advection','Forcing','Neutral Mixing','Vertical Mixing', ...
         'FOR+VMIX','ADV+VMIX','MIX','FOR+MIX','ADV+FOR'}, ...
         {'$\partial\Theta_\phi/\partial t$','Tendency','Advection','Forcing','Neutral Mixing','Vertical Mixing', ...
         'FOR+VMIX','ADV+VMIX','MIX','FOR+MIX','ADV+FOR'}};

vars = {'ADVp'};

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

% $$$ pranges = [10 90;
% $$$            10 90;
% $$$            0 100;]
% $$$ pranges = [0 10;
% $$$            0 10;
% $$$            0 100;]
% $$$ 
pranges = [0 100;
           0 100;
           0 100;]
% $$$ pranges = [90 100;
% $$$            90 100;
% $$$            90 100;]
pranges = [0 95;
           0 95;
           0 95;]
prangesi = pranges;
for pi=1:length(pranges(:))
    [tmp prangesi(pi)] = min(abs(P-pranges(pi)));
end

letlab = {'(b)','(c)','(d)'};

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
ylim([-5 -1]);
% $$$ ylim([-5 -2]);
% $$$ ylim([-6 -2]);
else
    pranges = [5 10; ...
               10 20; ...
               20 50; ...
              ];
    lnst = {'-','--',':','-.'};
    for i=1:length(pranges)
        [tmp Pilw] = min(abs(P-pranges(i,1)));
        [tmp Pihg] = min(abs(P-pranges(i,2)));
        plot(f,log10(mean(ZvP.Hp_spec(:,Pilw:Pihg),2)),'k','linewidth',2,'linestyle',lnst{i});
        hold on;
        plot(f,log10(mean(TvP.Hp_spec(:,Pilw:Pihg),2)),'r','linewidth',2,'linestyle',lnst{i});
% $$$         plot(f,log10(mean(TvP.Hp_spec(:,prangesi(2,1):prangesi(2,2)),2)),'-r','linewidth',2);
% $$$         plot(f,log10(mean(YvP.Hp_spec(:,prangesi(3,1):prangesi(3,2)),2)),'-b','linewidth',2);
    end
end
xlim([1/250 1]);
xlab = [300 100 50 10 7 5 3 2 1];
set(gca,'xscale','log');
set(gca,'xtick',1./xlab);
set(gca,'xticklabel',[]);
% $$$ xlabel('Period (years)');
legend('$\Theta_z$','$\Theta_\Theta$','$\Theta_\phi$','FontSize',10);
ylabel('Power $\log_{10}$($^\circ$C$^2$ year)');
title('(a) Temperature Anomalies');
grid on;
% $$$ title('ACCESS-CM2 PI-control Percentile-Averaged Spectra');
set(gca,'Position',poss(1,:));

% Spectral average (budget terms):
colors = {'m','b','k','r',[0 0.5 0],'c','y',[0.5 0 0],[0 0 0.5],[0.6 0.6 0.6],[0.2 0.2 0.2]};
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
box on;
set(gca,'xtick',1./xlab);
if (ti>=2)
set(gca,'xticklabel',xlab);
xlabel('Period (years)');
else
    set(gca,'xticklabel',[]);
end
legend(names{ti},'FontSize',12,'Location','SouthWest');
if (TorH)
% $$$     ylim([-21 -17.25]);
% $$$     ylim([-19 -16.25]);
    ylim([-20 -17.25]);
% $$$     ylim([-21 -18.25]);
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
title([letlab{ti} ' ' tnames{ti} ' Budget'],'Color',tcols{ti});
end

%%% Lag correlations:

[tmp ind1] = min(abs(P-10));[tmp ind2] = min(abs(P-30));
ts1 = mean(ZvP.Tp(ind1:ind2,:),1)';
ts1 = filter_field(ts1,12*10+1,'-t');
ts1(isnan(ts1)) = 0;

[tmp ind1] = min(abs(P-90));[tmp ind2] = min(abs(P-100));
ts2 = mean(TvP.Tp(ind1:ind2,:),1)';
ts2 = filter_field(ts2,12*10+1,'-t'); % Filtering
ts2(isnan(ts2)) = 0;

lags = [-100:1:100];
cor = zeros(size(lags));

for li = 1:length(lags)
    lag = lags(li)*12;

    ts_lagged = zeros(size(ts1));
    if (lag<0)
        ts_lagged(1:(end+lag)) = ts1((-lag+1):end);
        ts_lagged((end+lag+1):end) = 0;
    elseif (lag == 0)
        ts_lagged = ts1;
    else
        ts_lagged(1:(lag)) = 0;
        ts_lagged((lag+1):end) = ts1(1:(end-lag));
    end

    cor(li) = corr(ts_lagged,ts2);
end

[max_cor, max_lag] = max(cor);
[min_cor, min_lag] = min(cor);
max_lag = lags(max_lag);
min_lag = lags(min_lag);

ts_max = zeros(size(ts1));
lag = max_lag*12;
if (lag<0)
    ts_max(1:(end+lag)) = ts1((-lag+1):end);
    ts_max((end+lag+1):end) = 0;
elseif (lag == 0)
    ts_max = ts1;
else
    ts_max(1:(lag)) = 0;
    ts_max((lag+1):end) = ts1(1:(end-lag));
end

ts_min = zeros(size(ts1));
lag = min_lag*12;
if (lag<0)
    ts_min(1:(end+lag)) = ts1((-lag+1):end);
    ts_min((end+lag+1):end) = 0;
elseif (lag == 0)
    ts_min = ts1;
else
    ts_min(1:(lag)) = 0;
    ts_min((lag+1):end) = ts1(1:(end-lag));
end

figure;
subplot(2,1,1);
plot(time,ts1);
hold on;
plot(time,ts2,'-r');
plot(time,ts_max,'--k');
plot(time,ts_min,':k');
legend('$\Theta_z(10<p<30,t)$','$\Theta_\Theta(90<p<1000,t)$', ...
       ['$\Theta_z(10<p<30,t-%03d)$',max_lag), ...
       sprintf('$\Theta_z(10<p<30,t-%03d)$',min_lag));
subplot(2,1,2);
plot(lags,cor);
xlabel('Lag (years)');
ylabel('Correlation');


