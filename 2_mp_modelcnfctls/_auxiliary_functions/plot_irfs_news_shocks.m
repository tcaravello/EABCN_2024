function plot_irfs_news_shocks(model, Y_m_collector_f, Pi_m_collector_f, R_n_m_collector_f, m_f_use, m_f_use_index, label_ind,varargin)


global save_fig 

m_f_use_base = m_f_use;
m_f_use = m_f_use_base(m_f_use_index);

% results scaling

if label_ind == 0
    behav = '_m_d';
else
    behav = '_m_f';
end

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

%%
%----------------------------------------------------------------
% Figure
%----------------------------------------------------------------

linewidth_use = 3;
settings.colors.black  = [0 0 0];
settings.colors.grey   = [230/255 230/255 230/255];
settings.colors.dgrey  = [130/255 130/255 130/255];
settings.colors.blue   = [116/255 158/255 178/255];
settings.colors.lblue  = 0.25 * settings.colors.blue + 0.75 * [1 1 1];

IRF_hor_plot = 50;
% plot several news shocks.
news_shocks_plot = [0;5;10;20];
n_shocks = length(news_shocks_plot) ;
figure(2)
set(gcf, 'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
pos = get(gca, 'Position');
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex')
for j=1:n_shocks
    h_plot(j) = subplot(ceil(n_shocks/2),2,j);
    for i = 1:length(m_f_use)
    plot(1:1:IRF_hor_plot,scale*Y_m_collector_f(1:IRF_hor_plot, news_shocks_plot(j) +1,m_f_use_index(i)),'LineWidth',linewidth_use)
    hold on
    end
    hold off;
     title(strcat('Date \hspace{1mm}',num2str(news_shocks_plot(j))),'interpreter','latex','fontsize',16)  
    ylabel('%');
    xlim([0,IRF_hor_plot]);
        if j>=3
        xlabel('Horizon','interpreter','latex','FontSize',20)
    end

    if ceil(j/2)~=j/2
        ylabel('\% Deviation','interpreter','latex','FontSize',20)
    end
end
pause(0.001)
legend(beta_legend_use{3:end},'Location','best','fontsize',14,'interpreter','latex')
h_plot=axes('Position',[0.25,0,.5,.5],'Xlim',[0 1],'Ylim',[0 1]);
sgtitle('Output response to news shocks','interpreter','latex','fontsize',24)
set(h_plot,'Visible','off');

if save_fig == 1
 print([model_legend '_Y_news_shocks_many_behav' behav ], '-dpng')
end

figure(3)
set(gcf, 'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
pos = get(gca, 'Position');
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex')
for j=1:n_shocks
    h_plot(j) = subplot(ceil(n_shocks/2),2,j);
    for i = 1:length(m_f_use)
    plot(1:1:IRF_hor_plot,scale*4*Pi_m_collector_f(1:IRF_hor_plot, news_shocks_plot(j) +1,m_f_use_index(i)),'LineWidth',linewidth_use)
    hold on
    end
    hold off;
    title(strcat('Date \hspace{1mm}',num2str(news_shocks_plot(j))),'interpreter','latex','fontsize',16)  
    ylabel('%');
    xlim([0,IRF_hor_plot]);
        if j>=3
        xlabel('Horizon','interpreter','latex','FontSize',20)
    end

    if ceil(j/2)~=j/2
        ylabel('\% Deviation','interpreter','latex','FontSize',20)
    end
end
pause(0.001)
legend(beta_legend_use{3:end},'Location','best','fontsize',10,'interpreter','latex')
h_plot=axes('Position',[0.25,0,.5,.5],'Xlim',[0 1],'Ylim',[0 1]);
sgtitle('Inflation response to news shocks','interpreter','latex','fontsize',24)
set(h_plot,'Visible','off');

if save_fig == 1
 print([model_legend '_Pi_news_shocks_many_behav' behav ], '-dpng')
end


figure(4)
set(gcf, 'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
pos = get(gca, 'Position');
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex')
for j=1:n_shocks
    h_plot(j) = subplot(ceil(n_shocks/2),2,j);
    for i = 1:length(m_f_use)
    plot(1:1:IRF_hor_plot,scale*4*R_n_m_collector_f(1:IRF_hor_plot, news_shocks_plot(j) +1,m_f_use_index(i)),'LineWidth',linewidth_use)
    hold on
    end
    hold off;
     title(strcat('Date \hspace{1mm}',num2str(news_shocks_plot(j))),'interpreter','latex','fontsize',16) 
    ylabel('%');
    xlim([0,IRF_hor_plot]);
        if j>=3
        xlabel('Horizon','interpreter','latex','FontSize',20)
    end

    if ceil(j/2)~=j/2
        ylabel('\% Deviation','interpreter','latex','FontSize',20)
    end
end
pause(0.001)
legend(beta_legend_use{3:end},'Location','best','fontsize',10,'interpreter','latex')
h_plot=axes('Position',[0.25,0,.5,.5],'Xlim',[0 1],'Ylim',[0 1]);
sgtitle('Interest rate response to news shocks','interpreter','latex','fontsize',24)
set(h_plot,'Visible','off');

if save_fig == 1
 print([model_legend '_R_n_news_shocks_many_behav' behav], '-dpng')
end



if length(varargin) == 2
C_m_collector = varargin{1};
Inv_m_collector = varargin{2};


figure(9)
set(gcf, 'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
pos = get(gca, 'Position');
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex')
for j=1:n_shocks
    h_plot(j) = subplot(ceil(n_shocks/2),2,j);
    for i = 1:length(m_f_use)
    plot(1:1:IRF_hor_plot,scale*4*C_m_collector(1:IRF_hor_plot, news_shocks_plot(j) +1,m_f_use_index(i)),'LineWidth',linewidth_use)
    hold on
    end
    hold off;
     title(strcat('Date \hspace{1mm}',num2str(news_shocks_plot(j))),'interpreter','latex','fontsize',16)  
    ylabel('%');
    xlim([0,IRF_hor_plot]);
        if j>=3
        xlabel('Horizon','interpreter','latex','FontSize',20)
    end

    if ceil(j/2)~=j/2
        ylabel('\% Deviation','interpreter','latex','FontSize',20)
    end
end
pause(0.001)
legend(beta_legend_use{3:end},'Location','best','fontsize',10,'interpreter','latex')
h_plot=axes('Position',[0.25,0,.5,.5],'Xlim',[0 1],'Ylim',[0 1]);
sgtitle('Consumption response to news shocks')
set(h_plot,'Visible','off');

if save_fig == 1
 print([model_legend '_C_news_shocks_many_behav' behav], '-dpng')
end

figure(10)
set(gcf, 'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
pos = get(gca, 'Position');
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex')
for j=1:n_shocks
    h_plot(j) = subplot(ceil(n_shocks/2),2,j);
    for i = 1:length(m_f_use)
    plot(1:1:IRF_hor_plot,scale*4*Inv_m_collector(1:IRF_hor_plot, news_shocks_plot(j) +1,m_f_use_index(i)),'LineWidth',linewidth_use)
    hold on
    end
    hold off;
     title(strcat('Date \hspace{1mm}',num2str(news_shocks_plot(j))),'interpreter','latex','fontsize',16)  
    ylabel('%');
    xlim([0,IRF_hor_plot]);
        if j>=3
        xlabel('Horizon','interpreter','latex','FontSize',20)
    end

    if ceil(j/2)~=j/2
        ylabel('\% Deviation','interpreter','latex','FontSize',20)
    end
end
pause(0.001)
legend(beta_legend_use{3:end},'Location','best','fontsize',10,'interpreter','latex')
h_plot=axes('Position',[0.25,0,.5,.5],'Xlim',[0 1],'Ylim',[0 1]);
sgtitle('Investment response to news shocks')
set(h_plot,'Visible','off');

if save_fig == 1
 print([model_legend '_Inv_news_shocks_many_behav' behav], '-dpng')
end



end


end

