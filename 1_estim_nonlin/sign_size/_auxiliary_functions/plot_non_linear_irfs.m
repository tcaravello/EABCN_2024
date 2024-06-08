function plot_non_linear_irfs(irfs_plot, irfs_plot_linear,horizon,vars_plot,varNames_paper,non_lin_name, benchmark_linear,save_figs)

irf_coefs = irfs_plot.irf_coefs;
irf_se = irfs_plot.irf_se;
irf_t_stats = irfs_plot.irf_t_stats;
significance_factor_1 = irfs_plot.significance_factor_1;
significance_factor_2 = irfs_plot.significance_factor_2;




% get the t_stats of the non-linear terms
irf_coefs_linear = irf_coefs(1:2:end,:)';
irf_coefs_non_linear = irf_coefs(2:2:end,:)';
irf_se_linear = irf_se(1:2:end,:)';
irf_se_non_linear = irf_se(2:2:end,:)';
irf_var_mat = irfs_plot.irf_var_mat;
dataFrequency = 'M';
font_size_title = 18;
font_size_axis = 16;
nvar = length(vars_plot);
nvar_2 = size(irf_coefs_linear,2);
time = (0:horizon)';    % time horizon for IRFs

if strcmp(non_lin_name,'sign')
    str1 = 'positive shock';
    str2 = 'negative shock';
    mat_quad = [1 1;
                1 -1];
elseif strcmp(non_lin_name,'size')
        str1 = 'Small shock';
        str2 = 'Large shock';
        mat_quad = [1 0;
                1 1];
        if benchmark_linear == 1
            str1 = 'Linear';
        end
end


irf_se_max = zeros(horizon+1,nvar_2);
irf_se_min = zeros(horizon+1,nvar_2);



for j = 1:nvar_2
    for i = 1:horizon+1
        V_aux = mat_quad*irf_var_mat(:,:,i,j)*mat_quad';
        se_aux = sqrt(diag(V_aux));
        irf_se_max(i,j) = se_aux(1);
        irf_se_min(i,j) = se_aux(2);
    end
end

if benchmark_linear == 1

irf_coefs_linear_3 =  irfs_plot_linear.irf_coefs';
irf_se_linear_3 = irfs_plot_linear.irf_se';
% use 5% confidence values. 
irf_max = irf_coefs_linear_3;
irf_max_upper  = irf_coefs_linear_3+significance_factor_2*irf_se_linear_3;
irf_max_lower = irf_coefs_linear_3-significance_factor_2*irf_se_linear_3;
irf_min = irf_coefs_linear*mat_quad(2,1)+irf_coefs_non_linear*mat_quad(2,2);
irf_min_upper = irf_min+significance_factor_2*irf_se_min;
irf_min_lower = irf_min-significance_factor_2*irf_se_min;

else
irf_max = irf_coefs_linear*mat_quad(1,1)+irf_coefs_non_linear*mat_quad(1,2);
irf_min = irf_coefs_linear*mat_quad(2,1)+irf_coefs_non_linear*mat_quad(2,2);
irf_max_upper = irf_max+significance_factor_2*irf_se_max;
irf_max_lower = irf_max-significance_factor_2*irf_se_max;
irf_min_upper = irf_min+significance_factor_2*irf_se_min;
irf_min_lower = irf_min-significance_factor_2*irf_se_min;
end

%close all
figure('Position',[100 100 1000 600],'PaperPositionMode','Auto','DefaultAxesFontSize',13);
for j=1:nvar 
    jj = vars_plot(j);
    h_plot(j) = subplot(ceil(nvar/2),2,j);
    hh = fill([time(1); time(1:end); flipud([time(1:end);time(end)])],[irf_max_upper(1,jj); irf_max_lower(1:end,jj);flipud([irf_max_upper(1:end,jj); irf_max_lower(end,jj)])],[0.1, 0.4470,0.7410]); set(hh,'facealpha',.2); set(hh,'edgecolor','none');
    hold on
    hh_2 = fill([time(1); time(1:end); flipud([time(1:end);time(end)])] ...
        ,[irf_min_upper(1,jj); irf_min_lower(1:end,jj);flipud([irf_min_upper(1:end,jj); irf_min_lower(end,jj)])] ...
        ,[0.7410,0.4470, 0.1]); set(hh_2,'facealpha',.2); set(hh_2,'edgecolor','none');
    hold on
    p1= plot(time, irf_max(:,jj),'k', 'Linewidth', 1.5);
    hold on
    plot(time, irf_min(:,jj),'r', 'Linewidth', 1.5);
    if ~ismember(0,get(gca,'ylim'))
        line(get(gca,'xlim'),[0 0],'Color','k')
    end
    box on
    grid on ;hold off;
    title(varNames_paper{jj}, 'FontSize',font_size_title) 
    ylabel('\%');
    xlim([0,horizon]);
    set(gca,'FontSize',font_size_axis)
    if ceil(nvar/2)<j
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
sgtitle(strcat('Black line:',{' '},str1,'. Red Line:',{' '},str2),'FontSize',font_size_title)
if save_figs ==1
saveas(gcf,strcat('_results/irf_comparison_',non_lin_name),'png')
end


end