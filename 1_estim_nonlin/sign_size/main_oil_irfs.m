%% COMPUTE THE MAIN NON-LINEAR IRFS

% Tomas E. Caravello and Pedro Martinez Bruera

% this vesion June 7, 2024

%% HOUSEKEEPING
clear all
close all
clc

warning('off','MATLAB:dispatcher:nameConflict')
path = '/Users/tomyc/Dropbox (MIT)/EACBN_2024/code';
vintage = '';
task = '/1_estim_nonlin/sign_size';

addpath([path vintage task '/_auxiliary_functions'])
addpath([path vintage task '/_results'])
addpath([path '/_data'])
addpath(genpath([path '/_data/Kanzig']))
cd([path vintage task])

% initialize random number generator to be able to replicate results exactly
rng default

% Set text interpreter for figures to latex
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
%% Estimation settings
horizon = 48; %horizon for the local projection

lags = 18; %number of lags to use for controls

levels = 1; %1 if run the LP in levels, 0 if wanna use the long difference

non_linearity = 2; % 0 uses |x|, 2 is size non-linearity, 3 uses just the linear term

transformation = 1; %order of the polynomial to use.

use_own_lags_only = 0; %1 if use own lags only, 0 if you wanna use lags for all outcome variables in all regressions (like a VAR)
% if you chose 1, then pick in the next line which variables to use as
% controls in all equations (if you put zero then it is all)

original_data = 1; %1 loads Kanzig's original data

% {'Real oil price','Oil production','Oil stocks', 'World Industrial Production', 'Industrial Production', 'Consumer Price Index','Unemployment','Fed Funds Rate', 'Core CPI'};

vars_plot = [1,5,6,8]; %just plot real oil price, industrial production, CPI, and fed funds rate

trend = 0; %add a linear trend as Tenreyro and Thwaites (2016)

do_weights = 1; % 1 if you wanna compute and plot the weights

do_triple_test = 0; % 1 if you wanna compute and plot the weights

compute_sup_t = 0; % 1 if you wanna compute the sup-t bands in Montiel-Olea and Plagborg-Moeller (2019)

save_figs = 0; %0 if don't wanna save figs

% create non-linear function for size. Compute quantiles of z
opt = 0; % 0 means Newey-West with bandwidth=i; 1 means optimal bandwidth (takes much longer to run)'
har = 1; % 1 means Newey-West, 0 means just use heteroskedasticiy constitent as suggested by Plagborg-Moeller and Montiel Olea (ECTA, 2022) 
% note: all the output is ordered like: each column is a different horizon.
% And variables are order by two rows, that is: first two rows is the two
% coefficients in oil price, second two is the ones in oil production and
% so on.

alpha_1 = 0.32;
alpha_2 = 0.05;

significance_factor_1 = norminv(1-alpha_1/2,0,1); % 1 standard deviation

significance_factor_2 = norminv(1-alpha_2/2,0,1); %2 standard deviations 

benchmark_linear = 1; %1 if you want the size plot to use the linear (no non-linear regressor) as benchmark
%% Read in data
load_data_oil;

%% Histogram of the shocks
get_histogram_shocks;

%% Test symmetry
clc
if do_triple_test == 1
[h,p,s] = triplestest(z, 0.05);
% p is p-value, h = 0 is no reject, h = 1 is rejecting the null.
end
%%
opt_use = struct();
opt_use.order = [1,2,3,4]; %variables to use as controls in all equations
opt_use.level = levels;
opt_use.trend = trend;
opt_use.har = har;
opt_use.transformation = transformation;

if non_linearity==0
    z_vec = [z,abs(z)];
    non_lin_name = 'sign';
    %fun_sign = @(x) ((x>bound).*(x-bound)+(x<-bound).*(-x-bound));
    %z_vec = [z,fun_sign(z)];
elseif non_linearity==1 
    non_lin_name = 'sign';
    z_vec = [z,z.*z];
elseif non_linearity==2
    non_lin_name = 'size';
    z_vec_quantiles = quantile(z(z~=0),[0.2,0.8]);
    bound = max(abs(z_vec_quantiles));
    %bound = 0.5;
    fun_size = @(x) ((x>bound).*(x-bound)+(x<-bound).*(x+bound));
    if transformation == 1
        z_vec = [z,fun_size(z)];
    elseif transformation == 2
        z_vec = [z,z.*abs(z)];
    elseif transformation == 3
        z_vec = [z,z.^3];
    end
elseif non_linearity == 3
    non_lin_name = 'linear';
    z_vec = z;
end


irfs_plot_base = struct();
irfs_plot_linear= struct();

if use_own_lags_only == 1
    [irf_coefs,irf_se,irf_t_stats,irf_var_mat,extra_stuff] = run_non_linear_local_projections(data,z,[],lags,horizon+1,non_lin_name);
    [irfs_plot_linear.irf_coefs,irfs_plot_linear.irf_se,irfs_plot_linear.irf_t_stats,irfs_plot_linear.irf_var_mat,extra_stuff] = run_non_linear_local_projections(data,z,[],lags,horizon+1,'linear');
else
    [irf_coefs,irf_se,irf_t_stats,irf_var_mat,extra_stuff] = run_non_linear_local_projections(data,z,data,lags,horizon+1,non_lin_name,opt_use);
    [irfs_plot_linear.irf_coefs,irfs_plot_linear.irf_se,irfs_plot_linear.irf_t_stats,irfs_plot_linear.irf_var_mat,extra_stuff] = run_non_linear_local_projections(data,z,data,lags,horizon+1,'linear',opt_use);
end
 
irfs_plot_base.irf_coefs = irf_coefs;
irfs_plot_base.irf_se= irf_se;
irfs_plot_base.irf_t_stats = irf_t_stats;
irfs_plot_base.irf_var_mat = irf_var_mat;
irfs_plot_base.significance_factor_1 = significance_factor_1;
irfs_plot_base.significance_factor_2 = significance_factor_2;

%% Plot IRF with CI
close all
plot_irfs_coefs(irfs_plot_base,horizon,vars_plot,varNames_paper,non_lin_name,save_figs)

%% Non-linear IRFs
if ~strcmp(non_lin_name, 'linear')
plot_non_linear_irfs(irfs_plot_base, irfs_plot_linear,horizon,vars_plot,varNames_paper,non_lin_name, benchmark_linear,save_figs)
end
%% Compute weights

% re-write this as a function that takes some inputs and gives you some
% outputs.
if do_weights == 1
compute_and_plot_weights;
end

%% Sup-t bands
if compute_sup_t == 1
% Compute the variance covariance matrix as in Barnichon and Brownless (2019)
% This will give you the variance-covariance matrix of all parameters.
z_use = z_vec(lags+1:end,:);
level_collector = [];
se_collector = [];
sig_levels = [.9;.95;.99];
nvar_sup_t = length(vars_plot);
sup_t_cv_linear = zeros(size(data,2),size(sig_levels,1));
sup_t_cv_non_linear = zeros(size(data,2),size(sig_levels,1)); 
% 
for jj = 1:length(vars_plot)
    j = vars_plot(jj);
y_use = data(lags+1:end,j); %use the variable of interest
x_aux = lagmatrix(data(:,opt_use.order),1:lags);
if sum(opt_use.order==j)==0 %if the variable is not used as control for all, add lags of the same variable.
        y_lags_to_use = lagmatrix(data,1:lags);
        x_use_2 = [x_aux(lags+1:end,:), y_lags_to_use(lags+1:end,j)];
    else
        x_use_2 = x_aux(lags+1:end,:); %do not add lags for the same variable, because it is already included in x.
end 
% if use_own_lags_only == 1
%        x_use = lagmatrix(data(:,order),1:lags);
%        y_lags_to_use = lagmatrix(data(1:end,j),1:lags);
%        x_use_2 = [x_use(lags+1:end,:) y_lags_to_use(lags+1:end,:)];
% else
% x_mat = lagmatrix(data,1:lags);
% x_use_2 = x_mat(lags+1:end,:);
% end
irfs_joint_lp = locproj_2(y_use,z_use,x_use_2,0,horizon,'reg');
intervals_lp = locproj_conf_2(irfs_joint_lp,horizon,0);
% use the sup_t critical value
for i=1:size(sig_levels,1)
sup_t_cv_linear(j,i) =suptcritval_plugin(sig_levels(i),intervals_lp.var_cov_mat(1:horizon+1,1:horizon+1),1e5);
if non_linearity ~=3
sup_t_cv_non_linear(j,i) =suptcritval_plugin(sig_levels(i),intervals_lp.var_cov_mat(horizon+2:2*(horizon+1),horizon+2:2*(horizon+1)),1e5);
end
end
level_collector = [level_collector,irfs_joint_lp];
se_collector = [se_collector,intervals_lp];
end
%
irfs_plot_base_sup_t = struct();
irfs_plot_base_sup_t = irfs_plot_base;
irfs_plot_base_sup_t.sup_t_cv_linear = sup_t_cv_linear; 
irfs_plot_base_sup_t.sup_t_cv_non_linear = sup_t_cv_non_linear; 

plot_irfs_coefs_sup_t(irfs_plot_base_sup_t,horizon,vars_plot,varNames_paper,non_lin_name,save_figs)

end
%% Robustness on z
% number of observations "big shocks"
z_big_index = (z_vec(14:end,2)~=0);
z_small_index = logical((z_vec(14:end,1)~=0).*(z_vec(14:end,2)==0));
z_no_shock_index = (z_vec(14:end,1)==0);

% obtain interanual real oil price change/ IP/ inflation, nominal wages
data_control = data(13:end-1,:)-data(1:end-13,:);
data_control(:,7) = data(13:end-1,7);
data_control(:,8) = data(13:end-1,8);

%perform balance table.
fprintf('Variable & Large Shocks & Small Shocks & No Shock & Total\\\\\n')
for jj=1:length(vars_plot)
j = vars_plot(jj);
fprintf('%s & $%4.2f$ & $%4.2f$ & $%4.2f$ & $%4.2f$ \\\\ \n',varNames_paper{j},...
    mean(data_control(z_big_index,j)),...
    mean(data_control(z_small_index,j)),...
    mean(data_control(z_no_shock_index,j)),...
    mean(data_control(:,j)))
% standard errors
fprintf(' & $(%4.2f)$ & $(%4.2f)$ & $(%4.2f)$ & $(%4.2f)$ \\\\ \n',...
    std(data_control(z_big_index,j)),...
    std(data_control(z_small_index,j)),...
    std(data_control(z_no_shock_index,j)),...
    std(data_control(:,j)))
end
fprintf('N & $%4.0f$ & $%4.0f$ & $%4.0f$ & $%4.0f$ \\\\', sum(z_big_index),sum(z_small_index),sum(z_no_shock_index),size(z_vec,1))

