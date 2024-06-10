%% ESTIMATE MONETARY IRFs UNDER DIFFERENT ESTIMATORS

% This code estimates monetary IRFs using the estimators considered in Li,
% Plagborg-Møller and Wolf (2024, Journal of Econometrics)

% Tomas Caravello

% This version: June 7, 2024

%% HOUSEKEEPING
 
clc
clear all
close all

warning('off','MATLAB:dispatcher:nameConflict')
path = '/Users/tomyc/Dropbox (MIT)/EACBN_2024/code';
vintage = '';
task = '/1_estim_nonlin/estim';

addpath([path vintage task '/_aux_fun'])
addpath(genpath([path vintage task '/_aux_fun/Estimation_Routines']))
addpath([path '/_data'])

save_fig = 0;
%% Load Data
load_var_data;

startdate = find(date == 1969);
enddate = find(date == 2006.75);

extended_dataset = 0; %1 if you want to add more variables.

if extended_dataset == 1
    vardata_aux = [ad_shock gdp_tot  unemp cons inv lab cpi_inf tb3 lpcom];
    series_names = {'shock','Output','Unemployment','Consumption','Investment','Hours','Inflation','Interest rate','Commodity Prices'};
    str_fig = '_extended';
else
    vardata_aux = [ad_shock gdp_tot unemp cpi_inf tb3];
    series_names = {'shock','Output','Unemployment','Inflation','Interest rate'};
    str_fig = '_base';
end

outcomes_interest = {'Output','Unemployment','Inflation','Interest rate'};

vardata_aux = vardata_aux(startdate:enddate,:);
disp('Note: I am setting missing values of the shock to 0.')
vardata_aux(isnan(vardata_aux)) = 0;
data_y = vardata_aux;

data_sim.data_y = data_y;
%% Settings
settings = struct;

% horizon
settings.est.IRF_hor    = 21; % contemporaneous response (time 0) + 20 quarters in the future.

% identification method
settings.est.with_shock = 0; %this is: set the shock as perfect measure.
settings.est.recursive_shock = 1; %use recursive identification, the rr_shock is ordered third.
settings.est.with_IV = 0; % this is: set the shock as IV for SVAR.

shocks_pos = find(ismember(series_names,'shock')); 

settings.est.recursive_shock_pos = shocks_pos; % what position is the shock in?
settings.est.normalize_with_shock_std_dev = 1; % normalize shock to have std of 1.

% outcomes of interest
n_vars = length(outcomes_interest);
vars_irf = zeros(length(outcomes_interest),1);
for i = 1:n_vars
    vars_irf(i)= find(ismember(series_names,outcomes_interest{i})); %where the variables of interest are in terms of the IRF
end

% variables of interest
settings.est.IRF_response_var_pos = vars_irf; %variables of interest for the IRF
settings.est.est_normalize_var_pos = vars_irf(end); %normalize the response wrt the interest rate

%lags
settings.est.est_n_lag  = 0; %do not estimate the lag length
settings.est.n_lags_fix = 4; %fixed at 4 lags.
% THE FOLLOWING TWO LINES ONLY APPLY IF WE ESTIMATE THE LAG LENGTH
settings.est.n_lags_max = 8; %max number of lags if we choose it we the BIC
settings.est.est_n_lag_BIC  = 1; %1 uses the BIC, 0 uses AIC

%specific procedures:

% BVAR prior

settings.est.prior.towards_random_walk = 1; % prior shrinking towards random walk? otherwise towards zero
settings.est.bvar_glp                  = 1; % use Giannone, Lenza & Primiceri (2015) BVAR procedure?
                                            % otherwise use basic BVAR with default MN prior (see settings below)

% this sets priors for paramters, only used if bvar_glp == 0, otherwise the
% GLP procedure selects the prior.

settings.est.prior.tight_overall       = 0.04;
settings.est.prior.tight_nonown_lag    = 0.25;
settings.est.prior.decay_power         = 2;
settings.est.prior.tight_exogenous     = 1e5;

% BVAR posterior draws

settings.est.posterior_ndraw = 10000; % number of posterior draws (if set to 0, only use posterior mean of VAR coefficient to compute posterior mean of IRF)

% Penalized LP:

% we will try several lambdas, estimate in one subsample and then forecast
% in the rest. We pick the lambda to minimize the average forecast error.

settings.est.lambdaRange   = [0.001:0.005:0.021, 0.05:0.1:1.05, 2:1:19, 20:20:100, 200:200:2000]; % cross validation grid, scaled up by T
settings.est.irfLimitOrder = 2; % shrink towards polynomial of that order
settings.est.CV_folds      = 5; % Number of folds used for cross validation

% % VAR model averaging
% 
% settings.est.average_store_weight       = [2, 5, 9, 15, 21]; % store model weights at which horizon
% settings.est.average_store_submodel_irf = 0; % store IRF of each submodel? Only store if want to comput oracle weight
% settings.est.average_max_lags           = 1; % include lags up to n_lags_max? otherwise up to estimated lags
% settings.est.average_options            = optimoptions('quadprog','Display','off');

%% Estimate IRFs

% standard VAR

[IRF_VAR,~] ...
    = SVAR_est(data_sim,settings,0);

% Bias-corrected VAR using the Pope (1990) correction
[IRF_VAR_BC,~] ...
    = SVAR_est(data_sim,settings,1);

% Bayesian VAR
[IRF_BVAR,n_lags_est,~] = BVAR_est(data_sim,settings);
%%
% standard LP 

[IRF_LP,~] ...
    = LP_est(data_sim,settings,0);

% Bias-Corrected LP

[IRF_LP_BC,~] ...
    = LP_est(data_sim,settings,1);

% Penalized LP
[IRF_LP_pen,~, ~]...
= LP_shrink_est(data_sim,settings);
%% Plot Results

cd([path vintage task '/_results'])

close all

f = figure(position=[100 100 1200 800]);
time_vec = (0:(settings.est.IRF_hor-1))'; 

colors = struct();
colors.red    = [1 0 0]; %red
colors.lred    = 0.25 * colors.red  + 0.5*[1 1 1]; %light red
colors.dred    = 0.25 * colors.red  + 0.75*[0 0 0]; %dark red

colors.blue   = [116/255 158/255 178/255];
colors.lblue  = 0.25 * colors.blue + 0.5 * [1 1 1];
colors.dblue  = 0.25 * colors.blue + 0.75 * [0 0 0];

legend_vec = {'','VAR', 'BC VAR', 'BVAR', 'LP','BC LP', 'Pen LP'};

for jj = 1:length(outcomes_interest)
subplot(2,2,jj)
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
yline(0,'color','black','LineWidth',1.5,'LineStyle','--')
hold on
plot(time_vec,IRF_VAR(:,jj),'color',colors.red,'LineWidth',3)
hold on
plot(time_vec,IRF_VAR_BC(:,jj),'color',colors.red,'LineWidth',3,'LineStyle','--')
hold on
plot(time_vec,IRF_BVAR(:,jj),'color',colors.red,'LineWidth',3,'LineStyle',':')
hold on
plot(time_vec,IRF_LP(:,jj),'color',colors.blue,'LineWidth',3)
hold on
plot(time_vec,IRF_LP_BC(:,jj),'color',colors.blue,'LineWidth',3,'LineStyle','--')
hold on
plot(time_vec,IRF_LP_pen(:,jj),'color',colors.blue,'LineWidth',3,'LineStyle',':')
hold on
title(outcomes_interest{jj},'interpreter','latex','fontsize',20)
if jj == 4
legend(legend_vec,'location','best','FontSize',14,NumColumns=2)
end
end
if save_fig == 1
    print(['irfs_estim' str_fig],'-dpng');
end
