%% predictability of the proxy, and weights
% plot weights (regress both the linear and quadratic term on residuals,
% then compute weights)

% compute weights
[a_min, a_max] = bounds(z,'all');
size_grid = 200;
a_grid = linspace(-10,10,size_grid)';
nrep_boot_weights = 1000;
qe_boot_weights = [.1;.05];
[~,~,~,~,residualized_z] = run_local_projections_predictability(z_vec,data,lags,1,1,opt);
residualized_z = z_vec;
if non_linearity ~=3
%[irf_coefs_proxy,irf_se_proxy,irf_t_stats_proxy,irf_r_sq,residualized_z] = run_local_projections_predictability(z_vec,data,lags,1,1,opt);
% compute multiple regression coefficient, we'll need it for the weights.
cov_aux = cov(residualized_z);
beta_mr_weights = [cov_aux(1,2)/cov_aux(2,2);cov_aux(1,2)/cov_aux(1,1)];
% using those betas, define the random variables to compute covariances
% with
f_z_1 = residualized_z(:,1)-beta_mr_weights(1)*residualized_z(:,2);
f_z_2 = residualized_z(:,2)-beta_mr_weights(2)*residualized_z(:,1);
end


% [weights_1, ~] = compute_weights(a_grid,residualized_z(:,1),f_z_1,var(f_z_1));
% [weights_2, ~]  = compute_weights(a_grid,residualized_z(:,1),f_z_2,var(f_z_2));
% % decompose both in odd and even part
% weights_1_even = (weights_1+flip(weights_1))/2;
% weights_1_odd = (weights_1-flip(weights_1))/2;
% weights_2_even = (weights_2+flip(weights_2))/2;
% weights_2_odd = (weights_2-flip(weights_2))/2;

% bootstrap SE

if non_linearity==0
linear_comb_mat = [0.5 0.5;
                   0.5 -0.5];
elseif non_linearity==1
linear_comb_mat = mat_quad;
elseif non_linearity==2
   if size_non_linearity == 1
                    % compare linear with large shock
%     linear_comb_mat = [1 beta_mr_weights(2);
%                 1 1];
if benchmark_linear == 1
% linear_comb_mat = [1 1];
linear_comb_mat = [1 beta_mr_weights(2); % OK this worked!! maybe I can use it above to avoid running things twice.
                 1 1];
else
linear_comb_mat = mat_quad;
end
    elseif size_non_linearity == 2
            % compare +1% shock with +6% shock
    linear_comb_mat = [1 1;
                1 6];     
    elseif size_non_linearity == 3
            % compare +1% shock with +6% shock
    linear_comb_mat = [1 1;
                1 36];    
   end

elseif non_linearity==3
linear_comb_mat = 1;
end

if non_linearity==3
    results_weights_all = bootstrap_weights(residualized_z(:,1),residualized_z(:,1),nrep_boot_weights,a_grid,qe_boot_weights,linear_comb_mat);
else
    results_weights_all = bootstrap_weights(residualized_z(:,1),[f_z_1 f_z_2],nrep_boot_weights,a_grid,qe_boot_weights,linear_comb_mat);
end

% if benchmark_linear == 1
%     % compute the weights of the linear case and use them here.
%     [~,~,~,residualized_z_lin] = run_local_projections_predictability(z_vec_3,data,lags,1,1,opt);
% results_weights_all_linear = bootstrap_weights(residualized_z(:,1),residualized_z(:,1),nrep_boot_weights,a_grid,qe_boot_weights,1);
% end
%%
% plot weights
font_size = 7;
put_legend = 0; %1 if you want the legend in each of the top subplots
significance_level_weights = 10; %5 means 5 level, 10 means ten per cent
if significance_level_weights == 5
    low_ci_ind = 3;
    up_ci_ind = 4;
else
    low_ci_ind = 1;
    up_ci_ind = 2;
end


if non_linearity == 3
% weights and CIs
weights_1 = results_weights_all.weights(:,1);
%weights_1_low_ci = results_weights_all.ci_total(:,low_ci_ind,1);
%weights_1_up_ci = results_weights_all.ci_total(:,up_ci_ind,1);
weights_1_ci = results_weights_all.ci_total(:,:,1);
weights_1_even = results_weights_all.weights_even(:,1);
weights_1_odd = results_weights_all.weights_odd(:,1);
else
    % weights and CIs
weights_1 = results_weights_all.weights(:,1);
%weights_1_low_ci = results_weights_all.ci_total(:,low_ci_ind,1);
%weights_1_up_ci = results_weights_all.ci_total(:,up_ci_ind,1);
weights_1_ci = results_weights_all.ci_total(:,:,1);
weights_2 = results_weights_all.weights(:,2);
weights_2_ci = results_weights_all.ci_total(:,:,2);
%weights_2_low_ci = results_weights_all.ci_total(:,low_ci_ind,2);
%weights_2_up_ci = results_weights_all.ci_total(:,up_ci_ind,2);

% even and odd part
weights_1_even = results_weights_all.weights_even(:,1);
weights_2_even = results_weights_all.weights_even(:,2);
weights_1_odd = results_weights_all.weights_odd(:,1);
weights_2_odd = results_weights_all.weights_odd(:,2);

% weights and CIs: transformation
weights_1_trans = results_weights_all.weights_trans(:,1);
weights_1_ci_trans = results_weights_all.ci_total_trans(:,:,1);
%weights_1_low_ci_trans = results_weights_all.ci_total_trans(:,low_ci_ind,1);
%weights_1_up_ci_trans = results_weights_all.ci_total_trans(:,up_ci_ind,1);

weights_2_trans = results_weights_all.weights_trans(:,2);
weights_2_ci_trans = results_weights_all.ci_total_trans(:,:,2);
%weights_2_up_ci_trans = results_weights_all.ci_total_trans(:,up_ci_ind,2);
% even and odd: transformation
weights_1_even_trans = results_weights_all.weights_even_trans(:,1);
weights_2_even_trans = results_weights_all.weights_even_trans(:,2);
weights_1_odd_trans = results_weights_all.weights_odd_trans(:,1);
weights_2_odd_trans = results_weights_all.weights_odd_trans(:,2);
end

%close all 

% titles for the subplots
if non_linearity==0
    str_title_trans_1 = ', Positive';
    str_title_trans_2 = ', Negative';
elseif non_linearity==1 
    str_title_trans_1 = ", +"+num2str(mat_quad(1,2))+"\% shock";
    str_title_trans_2 = ", -"+num2str(mat_quad(2,2))+"\% shock";
elseif non_linearity==2
    if size_non_linearity == 1
if benchmark_linear == 1
    str_title_trans_1 = ', benchmark';
    str_title_trans_2 = ', large shock';
else
    str_title_trans_1 = ', small shock';
    str_title_trans_2 = ', large shock';
end

    elseif size_non_linearity == 2
    str_title_trans_1 = ", +" + num2str(mat_quad(1,2))+"\% shock";
    str_title_trans_2 = ", -" + num2str(mat_quad(2,2))+"\% shock";
    elseif size_non_linearity == 3
    str_title_trans_1 = ", +"+num2str(sqrt(mat_quad(1,2)))+" shock";
    str_title_trans_2 = ", -"+num2str(sqrt(mat_quad(2,2)))+" shock";
    end
elseif non_linearity == 3
end

if non_linearity == 3
figure;
subplot(2,1,1)
hh=fill([a_grid(1); a_grid(1:end); flipud([a_grid(1:end); a_grid(end)])],[weights_1_ci(1,low_ci_ind); weights_1_ci(1:end,low_ci_ind); flipud([weights_1_ci(1:end,up_ci_ind); weights_1_ci(end,up_ci_ind)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.4);
set(hh,'edgecolor','none'); 
hold on
plot(a_grid,weights_1, 'Color','blue','LineWidth',1.5)
yline(0,'k-')
if put_legend ==1
legend('95\% CI','weights','location','best','FontSize',font_size)
end
ylim([min(weights_1_ci(:,low_ci_ind))-0.05,max(weights_1_ci(:,up_ci_ind))+0.05])
title('Weights, linear regressor')
hold off


subplot(2,1,2)
plot(a_grid,weights_1, 'Color','blue','LineWidth',1.5)
hold on
plot(a_grid,weights_1_even, 'Color','black','LineWidth',1.1,'LineStyle','--')
hold on
plot(a_grid,weights_1_odd, 'Color','red','LineWidth',1.1,'LineStyle','--')
yline(0,'k-')
legend('total','even','odd','','location','best','FontSize',font_size)
ylim([min(weights_1)-0.05,max(weights_1)+0.05])
title('Weights, linear regressor')
hold off
if save_figs ==1
saveas(gcf,strcat('Figures_ours/weights_linear_',str_fig_save),'png')
end

else

figure;
subplot(2,2,1)
hh=fill([a_grid(1); a_grid(1:end); flipud([a_grid(1:end); a_grid(end)])],[weights_1_ci(1,low_ci_ind); weights_1_ci(1:end,low_ci_ind); flipud([weights_1_ci(1:end,up_ci_ind); weights_1_ci(end,up_ci_ind)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.4);
set(hh,'edgecolor','none'); 
hold on
plot(a_grid,weights_1, 'Color','blue','LineWidth',1.5)
yline(0,'k-')
if put_legend ==1
legend('95\% CI','weights','location','best','FontSize',font_size)
end
ylim([min(weights_1_ci(:,low_ci_ind))-0.05,max(weights_1_ci(:,up_ci_ind))+0.05])
title('Weights, linear regressor')
hold off

subplot(2,2,2)
hh=fill([a_grid(1); a_grid(1:end); flipud([a_grid(1:end); a_grid(end)])],[weights_2_ci(1,low_ci_ind); weights_2_ci(1:end,low_ci_ind); flipud([weights_2_ci(1:end,up_ci_ind); weights_2_ci(end,up_ci_ind)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.4);
set(hh,'edgecolor','none'); 
hold on
plot(a_grid,weights_2, 'Color','blue','LineWidth',1.5)
yline(0,'k-')
if put_legend ==1
legend('95\% CI','weights','location','best','FontSize',font_size)
end
ylim([min(weights_2_ci(:,low_ci_ind))-0.05,max(weights_2_ci(:,up_ci_ind))+0.05])
title(strcat('Weights ',str_title),'interpreter','Latex')
hold off

subplot(2,2,3)
plot(a_grid,weights_1, 'Color','blue','LineWidth',1.5)
hold on
plot(a_grid,weights_1_even, 'Color','black','LineWidth',1.1,'LineStyle','--')
hold on
plot(a_grid,weights_1_odd, 'Color','red','LineWidth',1.1,'LineStyle','--')
yline(0,'k-')
legend('total','even','odd','','location','best','FontSize',font_size)
ylim([min(weights_1)-0.05,max(weights_1)+0.05])
title('Weights, linear regressor')
hold off

subplot(2,2,4)
plot(a_grid,weights_2, 'Color','blue','LineWidth',1.5)
hold on
plot(a_grid,weights_2_even, 'Color','black','LineWidth',1.1,'LineStyle','--')
hold on
plot(a_grid,weights_2_odd, 'Color','red','LineWidth',1.1,'LineStyle','--')
ylim([min(weights_2)-0.05,max(weights_2)+0.05])
yline(0,'k--')
legend('total','even','odd','','location','best','FontSize',font_size)
title(strcat('Weights ',str_title),'interpreter','Latex')
hold off

if save_figs ==1
saveas(gcf,strcat('Figures_ours/weights_non_linear_original_',str_fig_save),'png')
end

% plot weights on the plotted IRFs

figure;
subplot(2,2,1)
hh=fill([a_grid(1); a_grid(1:end); flipud([a_grid(1:end); a_grid(end)])],[weights_1_ci_trans(1,low_ci_ind); weights_1_ci_trans(1:end,low_ci_ind); flipud([weights_1_ci_trans(1:end,up_ci_ind); weights_1_ci_trans(end,up_ci_ind)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.4);
set(hh,'edgecolor','none'); 
hold on
plot(a_grid,weights_1_trans, 'Color','blue','LineWidth',1.5)
yline(0,'k-')
if put_legend ==1
legend('95\% CI','weights','location','best','FontSize',font_size)
end
ylim([min(weights_1_ci_trans(:,low_ci_ind))-0.05,max(weights_1_ci_trans(:,up_ci_ind))+0.05])
title(strcat('Weights',str_title_trans_1),'interpreter','Latex')
hold off

subplot(2,2,2)
hh=fill([a_grid(1); a_grid(1:end); flipud([a_grid(1:end); a_grid(end)])],[weights_2_ci_trans(1,low_ci_ind); weights_2_ci_trans(1:end,low_ci_ind); flipud([weights_2_ci_trans(1:end,up_ci_ind); weights_2_ci_trans(end,up_ci_ind)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.4);
set(hh,'edgecolor','none'); 
hold on
plot(a_grid,weights_2_trans, 'Color','blue','LineWidth',1.5)
yline(0,'k-')
if put_legend ==1
legend('95\% CI','weights','location','best','FontSize',font_size)
end
ylim([min(weights_2_ci_trans(:,low_ci_ind))-0.05,max(weights_2_ci_trans(:,up_ci_ind))+0.05])
title(strcat('Weights',str_title_trans_2),'interpreter','Latex')
hold off



subplot(2,2,3)
plot(a_grid,weights_1_trans, 'Color','blue','LineWidth',1.5)
hold on
plot(a_grid,weights_1_even_trans, 'Color','black','LineWidth',1.1,'LineStyle','--')
hold on
plot(a_grid,weights_1_odd_trans, 'Color','red','LineWidth',1.1,'LineStyle','--')
yline(0,'k-')
legend('total','even','odd','','location','best','FontSize',font_size)
ylim([min([weights_1_trans,weights_1_even_trans,weights_1_odd_trans],[],'all')-0.05,max([weights_1_trans,weights_1_even_trans,weights_1_odd_trans],[],'all')+0.05])
title(strcat('Weights',str_title_trans_1),'interpreter','Latex')
hold off

subplot(2,2,4)
plot(a_grid,weights_2_trans, 'Color','blue','LineWidth',1.5)
hold on
plot(a_grid,weights_2_even_trans, 'Color','black','LineWidth',1.1,'LineStyle','--')
hold on
plot(a_grid,weights_2_odd_trans, 'Color','red','LineWidth',1.1,'LineStyle','--')
ylim([min([weights_2_trans,weights_2_even_trans,weights_2_odd_trans],[],'all')-0.05,max([weights_2_trans,weights_2_even_trans,weights_2_odd_trans],[],'all')+0.05])
yline(0,'k--')
legend('total','even','odd','','location','best','FontSize',font_size)
title(strcat('Weights',str_title_trans_2),'interpreter','Latex')
hold off
if save_figs ==1
saveas(gcf,strcat('Figures_ours/weights_non_linear_trans_',str_fig_save),'png')
end

end