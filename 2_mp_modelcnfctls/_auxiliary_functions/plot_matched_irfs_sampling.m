function plot_matched_irfs_sampling(Pi_m_fit_struct,Y_m_fit_struct,R_n_m_fit_struct,n_fig,save_name)

global T_use Y_m_lb_emp Y_m_ub_emp Y_m_emp  Pi_m_lb_emp Pi_m_ub_emp Pi_m_emp  R_n_m_lb_emp R_n_m_ub_emp R_n_m_emp save_fig
%----------------------------------------------------------------
% General Settings
%----------------------------------------------------------------

% plot horizon

IRF_hor_plot = T_use;

% results scaling

scale = 1;

% go to results folder

%----------------------------------------------------------------
% Figure
%----------------------------------------------------------------

settings.colors.black  = [0 0 0];
settings.colors.grey   = [230/255 230/255 230/255];
settings.colors.dgrey  = [130/255 130/255 130/255];
settings.colors.blue   = [116/255 158/255 178/255];
settings.colors.lblue  = 0.25 * settings.colors.blue + 0.75 * [1 1 1];
settings.colors.red   = [255, 0, 0]/255;
settings.colors.lred  = 0.25 *[1,0,0] + 0.75 * [1 1 1];

plotwidth = 0.27;
gapsize = 0.05;
gapsize_edges = (1-3*plotwidth-2*gapsize)/2;
left_pos = [gapsize_edges, gapsize_edges + gapsize + plotwidth, gapsize_edges + 2 * gapsize + 2 * plotwidth];


figure(n_fig)
set(gcf, 'units','normalized','outerposition',[0.0 0.0 0.9 0.9]);

subplot(1,3,1)
pos = get(gca, 'Position');
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
pos(1) = left_pos(1);
pos(3) = plotwidth;
set(gca,'Position', pos)
jbfill(1:1:IRF_hor_plot,scale*Y_m_lb_emp(1:IRF_hor_plot)',...
       scale*Y_m_ub_emp(1:IRF_hor_plot)',settings.colors.grey,settings.colors.grey,0,1);
hold on

jbfill(1:1:IRF_hor_plot,scale*Y_m_fit_struct.lb(1:IRF_hor_plot)',...
       scale*Y_m_fit_struct.ub(1:IRF_hor_plot)',settings.colors.blue,settings.colors.blue,0,0.5);
hold on

plot(1:1:IRF_hor_plot,scale*Y_m_emp(1:IRF_hor_plot),'Color',settings.colors.dgrey,'LineWidth',4)
hold on
plot(1:1:IRF_hor_plot,scale*Y_m_fit_struct.median(1:IRF_hor_plot),':','Color',settings.colors.red,'LineWidth',4)
hold on
xlim([1 IRF_hor_plot])
%ylim([-2 2])
%yticks([-2 -1 0 1 2])
set(gcf,'color','w')
title('Output','interpreter','latex','fontsize',24)
xlabel('Horizon','interpreter','latex','FontSize',20)
ylabel('\% Deviation','interpreter','latex','FontSize',20)
grid on
hold off


subplot(1,3,2)
pos = get(gca, 'Position');
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
pos(1) = left_pos(2);
pos(3) = plotwidth;
set(gca,'Position', pos)
jbfill(1:1:IRF_hor_plot,scale*Pi_m_lb_emp(1:IRF_hor_plot)',...
       scale*Pi_m_ub_emp(1:IRF_hor_plot)',settings.colors.grey,settings.colors.grey,0,1);
hold on

jbfill(1:1:IRF_hor_plot,scale*4*Pi_m_fit_struct.lb(1:IRF_hor_plot)',...
       scale*4*Pi_m_fit_struct.ub(1:IRF_hor_plot)',settings.colors.blue,settings.colors.blue,0,0.5);
hold on

plot(1:1:IRF_hor_plot,scale*Pi_m_emp(1:IRF_hor_plot),'Color',settings.colors.dgrey,'LineWidth',4)
hold on
plot(1:1:IRF_hor_plot,scale*4*Pi_m_fit_struct.median(1:IRF_hor_plot),':','Color',settings.colors.red,'LineWidth',4)
hold on
xlim([1 IRF_hor_plot])
%ylim([-2 2])
%yticks([-2 -1 0 1 2])
set(gcf,'color','w')
title('Inflation','interpreter','latex','fontsize',24)
xlabel('Horizon','interpreter','latex','FontSize',20)
grid on
hold off

subplot(1,3,3)
pos = get(gca, 'Position');
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
pos(1) = left_pos(3);
pos(3) = plotwidth;
set(gca,'Position', pos)
jbfill(1:1:IRF_hor_plot,scale*R_n_m_lb_emp(1:IRF_hor_plot)',...
       scale*R_n_m_ub_emp(1:IRF_hor_plot)',settings.colors.grey,settings.colors.grey,0,1);
hold on

jbfill(1:1:IRF_hor_plot,scale*4*R_n_m_fit_struct.lb(1:IRF_hor_plot)',...
       scale*4*R_n_m_fit_struct.ub(1:IRF_hor_plot)',settings.colors.blue,settings.colors.blue,0,0.5);
hold on

plot(1:1:IRF_hor_plot,scale*R_n_m_emp(1:IRF_hor_plot),'Color',settings.colors.dgrey,'LineWidth',4)
hold on
plot(1:1:IRF_hor_plot,scale*4*R_n_m_fit_struct.median(1:IRF_hor_plot),':','Color',settings.colors.red,'LineWidth',4)
hold on
xlim([1 IRF_hor_plot])
%ylim([-3 1])
%yticks([-3 -2 -1 0 1])
set(gcf,'color','w')
title('Interest Rate','interpreter','latex','fontsize',24)
xlabel('Horizon','interpreter','latex','FontSize',20)
legend({'','','Empirical', 'Median', 'Mode'},'Location','Northeast','fontsize',18,'interpreter','latex')
grid on
hold off

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 2.25*pos(3) 1*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');

if save_fig == 1
    print(save_name, '-dpng')
end
end