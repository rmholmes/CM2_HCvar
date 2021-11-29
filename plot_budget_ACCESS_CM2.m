% $$$ 
% $$$ 
% $$$ %%% Linear percentile plotting:
% $$$ colors = {'m','b','k','r',[0 0.5 0],'c'};     
% $$$     DER = 1;
% $$$     figure;
% $$$     set(gcf,'Position',[1921           1        1920        1005]);
% $$$     subplot(1,3,1);
% $$$     bvars = {'TEN_c','ADV_c','FOR_c','RMIX_c','VMIX_c'};
% $$$     for gi=1:length(bvars)
% $$$         eval(['var = mean(ZvP.' bvars{gi} ',2);']);
% $$$         if (DER)
% $$$             var = diff(var,[],1)./dP*100/rho0/Cp/mean(Vtot);
% $$$             plot(var/1e-9,P,'-','color',colors{gi},'linewidth',2);        
% $$$         else 
% $$$             plot(var/1e15,Pe,'-','color',colors{gi},'linewidth',2);
% $$$         end
% $$$         hold on;
% $$$     end
% $$$     plot([0 0],[0 100],'--k');
% $$$     ylabel('Percentile');
% $$$     set(gca,'ydir','reverse');
% $$$     ylim([0 100]);
% $$$     if (~DER)
% $$$         xlim([-2 2]);
% $$$         xlabel('Upward Vertical heat transport (PW)');
% $$$         legend('$\partial\mathcal{H}_z/\partial t$','$\mathcal{A}_z$','$\mathcal{F}_z$',...
% $$$                '$\mathcal{M}_z^{neutral}$',['$\' ...
% $$$                             'mathcal{M}_z^{vertical}$']);
% $$$     else
% $$$         xlabel('Temperature Tendency ($10^{-9}$ $^\circ$Cs $^{-1}$)');
% $$$         legend('$\partial\mathcal{\Theta}_z/\partial t$', ...
% $$$                '$\frac{100}{\rho_0C_p\mathcal{V}_T}\partial\mathcal{A}_z/\partial p_z$',...
% $$$                '$\frac{100}{\rho_0C_p\mathcal{V}_T}\partial\mathcal{F}_z/\partial p_z$',...
% $$$                '$\frac{100}{\rho_0C_p\mathcal{V}_T}\partial\mathcal{M}_z^{neutral}/\partial p_z$',...
% $$$                '$\frac{100}{\rho_0C_p\mathcal{V}_T}\partial\mathcal{M}_z^{vertical}/\partial p_z$');
% $$$     end
% $$$     set(gca,'Position',[0.06   0.1400    0.2580    0.8150]);
% $$$ % $$$ $100(\rho_0 C_p \mathcal{V}_T)^{-1} \partial \mathcal{F}_z/\partial p_z$    
% $$$     subplot(1,3,2);
% $$$ % $$$     TvP.FOR_c = TvP.FOR_c+2*TvP.JSH_c;
% $$$ % $$$     TvP.ADV_c = TvP.ADV_c+2*TvP.JSH_c;
% $$$ % $$$     Tv.FOR_c = Tv.FOR_c-TvP.JSH_c;
% $$$ % $$$     Tv.ADV_c = Tv.ADV_c-TvP.JSH_c;
% $$$     bvars = {'TEN_c','ADV_c','FOR_c','RMIX_c','VMIX_c'};
% $$$     for gi=1:length(bvars)
% $$$         eval(['var = mean(TvP.' bvars{gi} ',2);']);
% $$$         if (DER)
% $$$             var = diff(var,[],1)./dP*100/rho0/Cp/mean(Vtot);
% $$$             plot(var/1e-9,P,'-','color',colors{gi},'linewidth',2);        
% $$$         else 
% $$$             plot(var/1e15,Pe,'-','color',colors{gi},'linewidth',2);
% $$$         end
% $$$ % $$$         eval(['var = mean(Tv.' bvars{gi} ',2);']);
% $$$ % $$$         plot(var/1e15,Te,':','color',colors{gi},'linewidth',2);
% $$$         hold on;
% $$$     end
% $$$     plot([0 0],[0 100],'--k');
% $$$     set(gca,'ydir','reverse');
% $$$     ylim([0 100]);
% $$$     set(gca,'yticklabel',[]);
% $$$     if (~DER)
% $$$         xlabel('Cold-to-warm Diathermal heat transport (PW)');
% $$$         legend('$\partial\mathcal{H}_\Theta/\partial t$','$\mathcal{M}_\Theta^{numerical}$','$\mathcal{F}_\Theta$',...
% $$$            '$\mathcal{M}_\Theta^{neutral}$','$\mathcal{M}_\Theta^{vertical}$');
% $$$     else
% $$$         xlabel('Temperature Tendency ($10^{-9}$ $^\circ$Cs $^{-1}$)');
% $$$         legend('$\partial\mathcal{\Theta}_\Theta/\partial t$', ...
% $$$                '$\frac{100}{\rho_0C_p\mathcal{V}_T}\partial\mathcal{M}_\Theta^{numerical}/\partial p_\Theta$',...
% $$$                '$\frac{100}{\rho_0C_p\mathcal{V}_T}\partial\mathcal{F}_\Theta/\partial p_\Theta$',...
% $$$                '$\frac{100}{\rho_0C_p\mathcal{V}_T}\partial\mathcal{M}_\Theta^{neutral}/\partial p_\Theta$',...
% $$$                '$\frac{100}{\rho_0C_p\mathcal{V}_T}\partial\mathcal{M}_\Theta^{vertical}/\partial p_\Theta$');
% $$$     end
% $$$     set(gca,'Position',[0.3661    0.1400    0.2580    0.8150]);
% $$$ 
% $$$     subplot(1,3,3);
% $$$     bvars = {'TEN_c','ADV_c','FOR_c','RMIX_c'};
% $$$     for gi=1:length(bvars)
% $$$         eval(['var = mean(YvP.' bvars{gi} ',2);']);
% $$$         if (DER)
% $$$             var = diff(var,[],1)./dP*100/rho0/Cp/mean(Vtot);
% $$$             plot(var/1e-9,P,'-','color',colors{gi},'linewidth',2);        
% $$$         else 
% $$$             plot(var/1e15,Pe,'-','color',colors{gi},'linewidth',2);
% $$$         end
% $$$         hold on;
% $$$     end
% $$$     plot([0 0],[0 100],'--k');
% $$$     ylim([0 100]);
% $$$     if (~DER)
% $$$         xlabel('Southward Meridional heat transport (PW)');
% $$$         legend('$\partial\mathcal{H}_\phi/\partial t$','$\mathcal{A}_\phi^{advective}$','$\mathcal{F}_\phi$',...
% $$$                '$\mathcal{A}_\phi^{diffusive}$');
% $$$     else
% $$$         xlabel('Temperature Tendency ($10^{-9}$ $^\circ$Cs $^{-1}$)');
% $$$         legend('$\partial\mathcal{\Theta}_\phi/\partial t$', ...
% $$$                '$\frac{100}{\rho_0C_p\mathcal{V}_T}\partial\mathcal{A}_\phi^{advective}/\partial p_\phi$',...
% $$$                '$\frac{100}{\rho_0C_p\mathcal{V}_T}\partial\mathcal{F}_\phi/\partial p_\phi$',...
% $$$                '$\frac{100}{\rho_0C_p\mathcal{V}_T}\partial\mathcal{M}_\phi^{diffusive}/\partial p_\phi$');
% $$$     end
% $$$     set(gca,'Position',[0.7    0.1400    0.2580    0.8150]);
% $$$ 
    
    %%% Linear true variable plotting:
    
    DER = 0;
    degCc = 86400*365;
    
    colors = {'m','b','k','r',[0 0.5 0],'c'};     
    figure;
    set(gcf,'Position',[1921           1        1920        1005]);
    subplot(1,3,1);
    bvars = {'TEN_c','ADV_c','FOR_c','RMIX_c','VMIX_c'};
    for gi=1:length(bvars)
        eval(['var = mean(Zv.' bvars{gi} ',2);']);
        if (DER)
            var = -diff(var,[],1)/rho0/Cp./mean(Zv.V,2);
            plot(var*degCc,Z,'-','color',colors{gi},'linewidth',2);        
        else
            plot(var/1e15,Ze,'-','color',colors{gi},'linewidth',2);
        end
        hold on;
    end
    plot([0 0],[0 10000],'--k');
    ylabel('Depth (m)');
    set(gca,'ydir','reverse');
    ylim([0 500]);
    xlim([-2 2]);
    xlabel('Upward Vertical heat transport (PW)');
    legend('$\partial\mathcal{H}_z/\partial t$','$\mathcal{A}_z$','$\mathcal{F}_z$',...
           '$\mathcal{M}_z^{neutral}$',['$\' ...
                        'mathcal{M}_z^{vertical}$']);
    set(gca,'Position',[0.06   0.1400    0.2580    0.8150]);
    
    % Ticks:
    % 0-500m panel:
    yyaxis right;
    ylabel('Percentile');
    set(gca,'ydir','reverse');
    ylim([0 500]);
    Pticks = [0:2.5:12.5];
    Zticks = zeros(length(Pticks),1);
    for ii=1:length(Zticks)
        Zticks(ii) = interp1(mean(Zv.P,2),Ze, Pticks(ii),'linear');
    end
    set(gca,'ytick',Zticks);
    set(gca,'yticklabel',Pticks);

% $$$     % 500-5000m panel:
% $$$     yyaxis right;
% $$$     set(gca,'ydir','reverse');
% $$$     ylabel('Percentile');
% $$$     ylim([500 5000]);
% $$$     Pticks = [10:10:100];
% $$$     Zticks = zeros(length(Pticks),1);
% $$$     for ii=1:length(Zticks)
% $$$         Zticks(ii) = interp1(mean(Zv.P,2),Ze, Pticks(ii),'linear');
% $$$     end
% $$$     set(gca,'ytick',Zticks);
% $$$     set(gca,'yticklabel',Pticks);


    subplot(1,3,2);
    bvars = {'TEN_c','ADV_c','FOR_c','RMIX_c','VMIX_c'};
    for gi=1:length(bvars)
        eval(['var = mean(Tv.' bvars{gi} ',2);']);
        if (DER)
            var = -diff(var,[],1)/rho0/Cp./mean(Tv.V,2);
            plot(var*degCc,T,'-','color',colors{gi},'linewidth',2);        
        else
            plot(var/1e15,Te,'-','color',colors{gi},'linewidth',2);
        end
        hold on;
    end
    plot([0 0],[-2 34],'--k');
    ylim([-2 34]);
    xlabel('Cold-to-warm Diathermal heat transport (PW)');
    legend('$\partial\mathcal{H}_\Theta/\partial t$','$\mathcal{M}_\Theta^{numerical}$','$\mathcal{F}_\Theta$',...
           '$\mathcal{M}_\Theta^{neutral}$','$\mathcal{M}_\Theta^{vertical}$');

    yyaxis right;
    ylabel('Percentile');
    set(gca,'ydir','normal');
    ylim([-2 34]);
    Pticks = [0:2.5:10 20:10:50 80 95 99];
    Tticks = zeros(length(Pticks),1);
    for ii=1:length(Tticks)
        [vec inds] = unique(mean(Tv.P,2));
        Tticks(ii) = interp1(vec,Te(inds), Pticks(ii),'linear');
    end
    set(gca,'ytick',flipud(Tticks));
    set(gca,'yticklabel',fliplr(Pticks));


    
    subplot(1,3,3);
    bvars = {'TEN_c','ADV_c','FOR_c','RMIX_c'};
    for gi=1:length(bvars)
        eval(['var = mean(Yv.' bvars{gi} ',2);']);
        if (DER)
            var = -diff(var,[],1)/rho0/Cp./mean(Yv.V,2);
            plot(var*degCc,latv,'-','color',colors{gi},'linewidth',2);        
        else
            plot(var/1e15,latv_edges,'-','color',colors{gi},'linewidth',2);
        end
        hold on;
    end
    plot([0 0],[-90 90],'--k');
    ylim([-80 80]);
    xlim([-2 2]);
    ylabel('Latitude ($^\circ$N)');
    xlabel('Southward Meridional heat transport (PW)');
    legend('$\partial\mathcal{H}_\phi/\partial t$','$\mathcal{A}_\phi^{advective}$','$\mathcal{F}_\phi$',...
               '$\mathcal{A}_\phi^{diffusive}$');

    yyaxis right;
    ylabel('Percentile');
    set(gca,'ydir','normal');
    ylim([-80 80]);
    Pticks = [0 5 10:10:90 95];
    Yticks = zeros(length(Pticks),1);
    for ii=1:length(Yticks)
        [x ind] = unique(mean(Yv.P,2));
        Yticks(ii) = interp1(x,latv_edges(ind), Pticks(ii),'linear');
    end
    set(gca,'ytick',Yticks);
    set(gca,'yticklabel',Pticks);
