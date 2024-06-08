function plot_matched_irfs(model, Y_m_fit_collector_f, Pi_m_fit_collector_f, R_n_m_fit_collector_f,m_f_use, m_f_use_index, label_ind)

global T_use Y_m_lb_emp Y_m_ub_emp Y_m_emp  Pi_m_lb_emp Pi_m_ub_emp Pi_m_emp  R_n_m_lb_emp R_n_m_ub_emp R_n_m_emp save_fig

% plot horizon

m_f_use_base = m_f_use;
m_f_use = m_f_use_base(m_f_use_index);


if label_ind == 0
    behav = '_m_d';
else
    behav = '_m_f';
end


IRF_hor_plot = T_use;

% results scaling

scale = 1;

if strcmp(model, '/rank')
    model_legend = 'RANK';
elseif strcmp(model, '/hybrid_olg')
       model_legend = 'OLG';
elseif strcmp(model, '/hank')
    model_legend = 'HANK';
end

beta_legend_use = cell(length(m_f_use)+2,1);
beta_legend_use{1} = '';
beta_legend_use{2} = 'Empirical';
for j = 1:length(m_f_use)
    if m_f_use(j)>0.999
        beta_legend_use{j+2} = [model_legend ', $m = $' num2str(1)];
    else
        beta_legend_use{j+2} = [model_legend ', $m = $' num2str(m_f_use(j))];
    end
end

%----------------------------------------------------------------
% Figure
%----------------------------------------------------------------

settings.colors.black  = [0 0 0];
settings.colors.grey   = [230/255 230/255 230/255];
settings.colors.dgrey  = [130/255 130/255 130/255];
settings.colors.blue   = [116/255 158/255 178/255];
settings.colors.lblue  = 0.25 * settings.colors.blue + 0.75 * [1 1 1];

color_mat = [       0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840];


plotwidth = 0.27;
gapsize = 0.05;
gapsize_edges = (1-3*plotwidth-2*gapsize)/2 +0.02;
left_pos = [gapsize_edges, gapsize_edges + gapsize + plotwidth, gapsize_edges + 2 * gapsize + 2 * plotwidth];


close all

figure(1)
set(gcf, 'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);

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
plot(1:1:IRF_hor_plot,scale*Y_m_emp(1:IRF_hor_plot),':','Color',settings.colors.dgrey,'LineWidth',4)
hold on

for j = 1:length(m_f_use)
plot(1:1:IRF_hor_plot,scale*Y_m_fit_collector_f(1:IRF_hor_plot,m_f_use_index(j)),'LineWidth',4,'Color',color_mat(j,:))
hold on
end

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
plot(1:1:IRF_hor_plot,scale*Pi_m_emp(1:IRF_hor_plot),':','Color',settings.colors.dgrey,'LineWidth',4)
hold on
for j = 1:length(m_f_use)
plot(1:1:IRF_hor_plot,scale*4*Pi_m_fit_collector_f(1:IRF_hor_plot,m_f_use_index(j)),'LineWidth',4, 'Color',color_mat(j,:))
hold on
end
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
plot(1:1:IRF_hor_plot,scale*R_n_m_emp(1:IRF_hor_plot),':','Color',settings.colors.dgrey,'LineWidth',4)
hold on
for j = 1:length(m_f_use)
plot(1:1:IRF_hor_plot,scale*4*R_n_m_fit_collector_f(1:IRF_hor_plot,m_f_use_index(j)),'LineWidth',4, 'Color',color_mat(j,:))
hold on
end
xlim([1 IRF_hor_plot])
set(gcf,'color','w')
title('Interest Rate','interpreter','latex','fontsize',24)
xlabel('Horizon','interpreter','latex','FontSize',20)
legend(beta_legend_use,'Location','Northeast','fontsize',14,'interpreter','latex')
grid on
hold off

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 2.25*pos(3) 1*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
if save_fig == 1
    print([model_legend '_IRFs_many_behav'  behav], '-dpng')
end

end