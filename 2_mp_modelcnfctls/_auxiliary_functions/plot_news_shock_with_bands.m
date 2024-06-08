function plot_news_shock_with_bands(model,Y_m_struct,Pi_m_struct,R_n_m_struct,n_fig)

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
settings.colors.lblue  = 0.25 * settings.colors.blue + 0.75 * [1 1 1];
settings.colors.red   = [255, 0, 0]/255;
settings.colors.lred  = 0.25 *[1,0,0] + 0.75 * [1 1 1];

plotwidth = 0.27;
gapsize = 0.05;
gapsize_edges = (1-3*plotwidth-2*gapsize)/2;
left_pos = [gapsize_edges, gapsize_edges + gapsize + plotwidth, gapsize_edges + 2 * gapsize + 2 * plotwidth];

% plot several news shocks.
news_shocks_plot = [0;4;9;14];
n_shocks = length(news_shocks_plot);

figure(n_fig)
set(gcf, 'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
for j=1:n_shocks
    h_plot(j) = subplot(ceil(n_shocks/2),2,j);
    jbfill(1:1:IRF_hor_plot,scale*Y_m_struct.lb(1:IRF_hor_plot, news_shocks_plot(j)+1)',...
           scale*Y_m_struct.ub(1:IRF_hor_plot,news_shocks_plot(j)+1)',settings.colors.blue,settings.colors.blue,0,0.5);
    hold on
    plot(1:1:IRF_hor_plot,scale*Y_m_struct.median(1:IRF_hor_plot, news_shocks_plot(j)+1),':','Color',settings.colors.red,'LineWidth',4)
    hold on
    %plot(1:1:IRF_hor_plot,scale*Y_m_mat(1:IRF_hor_plot, news_shocks_plot(j)+1),'Color',settings.colors.blue,'LineWidth',4)
    %hold on
    hold off;
    title(strcat('hor =',num2str(news_shocks_plot(j)))) 
    ylabel('%');
    xlim([0,IRF_hor_plot]);
    grid on
end
pause(0.001)
legend({'','Median','Mode'},'Location','best','fontsize',14,'interpreter','latex')
h_plot=axes('Position',[0.25,0,.5,.5],'Xlim',[0 1],'Ylim',[0 1]);
sgtitle('Output response to news shocks')
set(h_plot,'Visible','off');

if save_fig == 1
    %print([model(2:end) '_Y_news_shocks'], '-dpng')
end


figure(n_fig+1)
set(gcf, 'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
for j=1:n_shocks
    h_plot(j) = subplot(ceil(n_shocks/2),2,j);
    jbfill(1:1:IRF_hor_plot,scale*4*Pi_m_struct.lb(1:IRF_hor_plot, news_shocks_plot(j)+1)',...
           scale*4*Pi_m_struct.ub(1:IRF_hor_plot,news_shocks_plot(j)+1)',settings.colors.blue,settings.colors.blue,0,0.5);
    hold on
    plot(1:1:IRF_hor_plot,scale*4*Pi_m_struct.median(1:IRF_hor_plot, news_shocks_plot(j)+1),':','Color',settings.colors.red,'LineWidth',4)
    hold on
    %plot(1:1:IRF_hor_plot,scale*4*Pi_m_mat(1:IRF_hor_plot, news_shocks_plot(j)+1),'Color',settings.colors.blue,'LineWidth',4)
    %hold on
    hold off;
    title(strcat('hor =',num2str(news_shocks_plot(j)))) 
    ylabel('%');
    xlim([0,IRF_hor_plot]);
    grid on
end
pause(0.001)
legend({'','Median','Mode'},'Location','best','fontsize',14,'interpreter','latex')
h_plot=axes('Position',[0.25,0,.5,.5],'Xlim',[0 1],'Ylim',[0 1]);
sgtitle('Inflation response to news shocks')
set(h_plot,'Visible','off');

if save_fig == 1
%print([model(2:end) '_Pi_news_shocks'], '-dpng')
end



figure(n_fig+2)
set(gcf, 'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
for j=1:n_shocks
    h_plot(j) = subplot(ceil(n_shocks/2),2,j);
    jbfill(1:1:IRF_hor_plot,scale*4*R_n_m_struct.lb(1:IRF_hor_plot, news_shocks_plot(j)+1)',...
           scale*4*R_n_m_struct.ub(1:IRF_hor_plot,news_shocks_plot(j)+1)',settings.colors.blue,settings.colors.blue,0,0.5);
    hold on
    plot(1:1:IRF_hor_plot,scale*4*R_n_m_struct.median(1:IRF_hor_plot, news_shocks_plot(j)+1),':','Color',settings.colors.red,'LineWidth',4)
    hold on
    %plot(1:1:IRF_hor_plot,scale*4*R_n_m_mat(1:IRF_hor_plot, news_shocks_plot(j)+1),'Color',settings.colors.blue,'LineWidth',4)
    %hold on
    hold off;
    title(strcat('hor =',num2str(news_shocks_plot(j)))) 
    ylabel('%');
    xlim([0,IRF_hor_plot]);
    grid on
end
pause(0.001)
legend({'','Median','Mode'},'Location','best','fontsize',14,'interpreter','latex')
h_plot=axes('Position',[0.25,0,.5,.5],'Xlim',[0 1],'Ylim',[0 1]);
sgtitle('Nominal rate response to news shocks')
set(h_plot,'Visible','off');

if save_fig == 1
    %print([model(2:end) '_R_n_news_shocks'], '-dpng')
end
end