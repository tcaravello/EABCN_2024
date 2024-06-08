function plot_inferred_shocks(m_fit_struct)

global save_fig

% plot horizon

IRF_hor_plot = 25;

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
settings.colors.dblue   = [116/255 158/255 178/255]*0.6;
settings.colors.lblue  = 0.25 * settings.colors.blue + 0.75 * [1 1 1];
settings.colors.red   = [255, 0, 0]/255;
settings.colors.lred  = 0.25 *[1,0,0] + 0.75 * [1 1 1];


figure(20)
pos = get(gca, 'Position');
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
jbfill(1:1:IRF_hor_plot,scale*m_fit_struct.lb(1:IRF_hor_plot)',...
       scale*m_fit_struct.ub(1:IRF_hor_plot)',settings.colors.blue,settings.colors.blue,0,0.5);
hold on
plot(1:1:IRF_hor_plot,scale*m_fit_struct.median(1:IRF_hor_plot),':','Color',settings.colors.red,'LineWidth',4)
grid on
title('Implicit news shocks')
legend({'','Median'},'Location','best','fontsize',14,'interpreter','latex')

if save_fig == 1
%print('models_posterior_inferred_news_shocks', '-dpng')
end
end