function plot_irfs_coefs_sup_t(irfs_plot,horizon,vars_plot,varNames_paper,non_lin_name,save_figs)

irf_coefs = irfs_plot.irf_coefs;
irf_se = irfs_plot.irf_se;
irf_t_stats = irfs_plot.irf_t_stats;
sup_t_cv_linear = irfs_plot.sup_t_cv_linear;
sup_t_cv_non_linear = irfs_plot.sup_t_cv_non_linear;


font_size_title = 18;
font_size_axis = 16;

if strcmp(non_lin_name,'linear')
% get the t_stats of the non-linear terms
irf_coefs_linear = irf_coefs';
irf_se_linear = irf_se';
t_stats_linear = irf_t_stats';
% use 5% confidence values. 
irf_se_linear_upper_1 = irf_coefs_linear+(ones(horizon+1,1)*sup_t_cv_linear(:,1)').*irf_se_linear;
irf_se_linear_lower_1 = irf_coefs_linear-(ones(horizon+1,1)*sup_t_cv_linear(:,1)').*irf_se_linear;
irf_se_linear_upper_2 = irf_coefs_linear+(ones(horizon+1,1)*sup_t_cv_linear(:,2)').*irf_se_linear;
irf_se_linear_lower_2 = irf_coefs_linear-(ones(horizon+1,1)*sup_t_cv_linear(:,2)').*irf_se_linear;
else
% get the t_stats of the non-linear terms
irf_coefs_linear = irf_coefs(1:2:end,:)';
irf_coefs_non_linear = irf_coefs(2:2:end,:)';
irf_se_linear = irf_se(1:2:end,:)';
irf_se_non_linear = irf_se(2:2:end,:)';
t_stats_linear = irf_t_stats(1:2:end,:)';
t_stats_non_linear = irf_t_stats(2:2:end,:)';
% use 5% confidence values.
irf_se_linear_upper_1 = irf_coefs_linear+(ones(horizon+1,1)*sup_t_cv_linear(:,1)').*irf_se_linear;
irf_se_linear_lower_1 = irf_coefs_linear-(ones(horizon+1,1)*sup_t_cv_linear(:,1)').*irf_se_linear;
irf_se_non_linear_upper_1 = irf_coefs_non_linear+(ones(horizon+1,1)*sup_t_cv_non_linear(:,1)').*irf_se_non_linear;
irf_se_non_linear_lower_1 = irf_coefs_non_linear-(ones(horizon+1,1)*sup_t_cv_non_linear(:,1)').*irf_se_non_linear;
irf_se_linear_upper_2 = irf_coefs_linear+(ones(horizon+1,1)*sup_t_cv_linear(:,2)').*irf_se_linear;
irf_se_linear_lower_2 = irf_coefs_linear-(ones(horizon+1,1)*sup_t_cv_linear(:,2)').*irf_se_linear;
irf_se_non_linear_upper_2 = irf_coefs_non_linear+(ones(horizon+1,1)*sup_t_cv_non_linear(:,2)').*irf_se_non_linear;
irf_se_non_linear_lower_2 = irf_coefs_non_linear-(ones(horizon+1,1)*sup_t_cv_non_linear(:,2)').*irf_se_non_linear;
end
dataFrequency = 'M';

time = (0:horizon)';    % time horizon for IRFs
nvar = length(vars_plot);

if strcmp(non_lin_name,'linear')

for j=1:nvar
    jj = vars_plot(j);
    h_plot(j) = subplot(ceil(nvar/2),2,j);
    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[irf_se_linear_upper_1(1,jj); irf_se_linear_lower_1(1:end,jj); flipud([irf_se_linear_upper_1(1:end,jj); irf_se_linear_lower_1(end,jj)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.4);
    set(hh,'edgecolor','none'); 
    hold on
    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[irf_se_linear_upper_2(1,jj); irf_se_linear_lower_2(1:end,jj); flipud([irf_se_linear_upper_2(1:end,jj); irf_se_linear_lower_2(end,jj)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.2);
    set(hh,'edgecolor','none'); 
    hold on
    p1=plot(time, irf_coefs_linear(:,jj),'k', 'Linewidth', 1.5); 
    if ~ismember(0,get(gca,'ylim'))
        line(get(gca,'xlim'),[0 0],'Color','k')
    end

    box on
    grid on ;hold off;
    title(varNames_paper{jj},'FontSize',font_size_title) 
    ylabel('\%');
    xlim([0,horizon]);
    set(gca,'FontSize',font_size_axis)
    if j>ceil(nvar/2)
    if dataFrequency == 'M'
        xlabel('Months');
        xticks(0:10:horizon)
    elseif dataFrequency == 'Q'
        xlabel('Quarters');
        xticks(0:5:horizon)
    end   
    end
end
pause(0.001)
h_plot=axes('Position',[0.25,0,.5,.5],'Xlim',[0 1],'Ylim',[0 1]);
set(h_plot,'Visible','off');

if save_figs == 1
    saveas(gcf,strcat('_results/irf_linear_coef_',non_lin_name,'_sup_t'),'png')
end



else

figure('Position',[100 100 1000 600],'PaperPositionMode','Auto','DefaultAxesFontSize',13);
for j=1:nvar 
    jj = vars_plot(j);
    h_plot(j) = subplot(ceil(nvar/2),2,j);
    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[irf_se_linear_upper_1(1,jj); irf_se_linear_lower_1(1:end,jj); flipud([irf_se_linear_upper_1(1:end,jj); irf_se_linear_lower_1(end,jj)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.4);
    set(hh,'edgecolor','none'); 
    hold on
    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[irf_se_linear_upper_2(1,jj); irf_se_linear_lower_2(1:end,jj); flipud([irf_se_linear_upper_2(1:end,jj); irf_se_linear_lower_2(end,jj)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.2);
    set(hh,'edgecolor','none'); 
    hold on
    p1=plot(time, irf_coefs_linear(:,jj),'k', 'Linewidth', 1.5); 
    if ~ismember(0,get(gca,'ylim'))
        line(get(gca,'xlim'),[0 0],'Color','k')
    end

    box on
    grid on ;hold off;
    title(strcat(varNames_paper{jj},', linear'), 'FontSize',font_size_title) 
    ylabel('\%');
    xlim([0,horizon]);
    set(gca,'FontSize',font_size_axis)
    if j>ceil(nvar/2)
    if dataFrequency == 'M'
        xlabel('Months');
        xticks(0:10:horizon)
    elseif dataFrequency == 'Q'
        xlabel('Quarters');
        xticks(0:5:horizon)
    end 
    end
end
pause(0.001)
h_plot=axes('Position',[0.25,0,.5,.5],'Xlim',[0 1],'Ylim',[0 1]);
set(h_plot,'Visible','off');
if save_figs ==1
saveas(gcf,strcat('_results/irf_linear_coef_',non_lin_name,'_sup_t'),'png')
end

figure('Position',[100 100 1000 600],'PaperPositionMode','Auto','DefaultAxesFontSize',13);
for j=1:nvar 
    jj = vars_plot(j);
    h_plot(j) = subplot(ceil(nvar/2),2,j);
    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[irf_se_non_linear_upper_1(1,jj); irf_se_non_linear_lower_1(1:end,jj); flipud([irf_se_non_linear_upper_1(1:end,jj); irf_se_non_linear_lower_1(end,jj)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.4);
    set(hh,'edgecolor','none'); 
    hold on
    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[irf_se_non_linear_upper_2(1,jj); irf_se_non_linear_lower_2(1:end,jj); flipud([irf_se_non_linear_upper_2(1:end,jj); irf_se_non_linear_lower_2(end,jj)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.2);
    set(hh,'edgecolor','none'); 
    hold on
    p1=plot(time, irf_coefs_non_linear(:,jj),'k', 'Linewidth', 1.5); 
    if ~ismember(0,get(gca,'ylim'))
        line(get(gca,'xlim'),[0 0],'Color','k')
    end

    box on
    grid on ;hold off;
    title(strcat(varNames_paper{jj},', non-linear'),'interpreter','latex', 'FontSize',font_size_title) 
    ylabel('\%');
    xlim([0,horizon]);
    set(gca,'FontSize',font_size_axis)
    if j>ceil(nvar/2)
    if dataFrequency == 'M'
        xlabel('Months');
        xticks(0:10:horizon)
    elseif dataFrequency == 'Q'
        xlabel('Quarters');
        xticks(0:5:horizon)
    end   
    end
end
pause(0.001)
h_plot=axes('Position',[0.25,0,.5,.5],'Xlim',[0 1],'Ylim',[0 1]);
set(h_plot,'Visible','off');
if save_figs ==1
saveas(gcf,strcat('_results/irf_non_linear_coef_',non_lin_name,'_sup_t'),'png')
end


end