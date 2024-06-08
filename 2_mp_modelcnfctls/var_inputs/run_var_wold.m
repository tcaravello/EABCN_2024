%% WOLD SHOCKS IRFs
% Tomas Caravello, Alisdair McKay, Christian Wolf

% this version: June 7, 2024

%% HOUSEKEEPING
 
clc
clear all
close all

warning('off','MATLAB:dispatcher:nameConflict')

path = '/Users/tomyc/Dropbox (MIT)/EACBN_2024/code';
session = '/2_mp_modelcnfctls';
task = '/var_inputs';

addpath([path session '/_auxiliary_functions'])
addpath([path '/_data/cmw'])

save_results = 1;

plot_wold_irfs = 0; % 1 if you want to plot the IRFs, 0 if you don't, then the code just produces the results.

cd([path session task]);


%% DATA

% import data

data_table = readtable('_data_cmw.csv',detectImportOptions('_data_cmw.csv'));
data = table2array(data_table);
date = data(:,1);

series = data_table.Properties.VariableNames;
series = series(2:end);

% raw macro outcomes

gdp       = data(:,2);
unemp     = data(:,3);
ffr       = 4 * data(:,4);
infl      = 100 * 4 * data(:,5);
inv       = data(:,6);
cons      = data(:,7);
lab       = data(:,8);
lab_share = data(:,9);
lab_prod  = data(:,10);
tfp       = data(:,11);

% macro outcome transformations: Hamilton Filter the desired variables.

gdp       = 100 * stat_transform(gdp,1);
inv       = 100 * stat_transform(inv,1);
cons      = 100 * stat_transform(cons,1);
lab       = 100 * stat_transform(lab,1);
lab_share = 100 * stat_transform(lab_share,1);
lab_prod  = 100 * stat_transform(lab_prod,1);
tfp       = 100 * stat_transform(tfp,1);

% dates

startdate = find(date == 1960);
enddate = find(date == 2019.75);

% collect VAR inputs

vardata = [unemp gdp inv cons lab tfp lab_prod lab_share infl ffr]; 
series_names = {'Unemployment', 'Output', 'Investment', 'Consumption','Hours',...
    'TFP', 'Labor Productivity', 'Labor Share', 'Inflation', 'FFR'};

vardata = vardata(startdate:enddate,:);

% de-trend: this just demeans the data for inflation, FFR and unemployment.

vardata = detrend(vardata);

%% SETTINGS

n_lags     = 4;                    % number of lags

constant   = 0;                    % constant? 0: no constant, 1: just constant, 2: constant and linear trend.
% no need to add the constant again since we already detrended the data.

IRF_hor    = 250;
n_draws    = 50; %use small number of draws just for illustration purposes.
n_y        = size(vardata,2);

%% VAR ESTIMATION

%----------------------------------------------------------------
% Estimate Reduced-Form VAR
%----------------------------------------------------------------

% Prior is very loose so this is essentially equivalent to standard VAR.

T = size(vardata,1) - n_lags;
[B_draws,Sigma_draws,B_OLS,Sigma_OLS] = bvar_fn(vardata,n_lags,constant,n_draws);

%----------------------------------------------------------------
% OLS Wold IRFs
%----------------------------------------------------------------

% extract VAR inputs
    
Sigma_u   = Sigma_OLS;
B         = B_OLS;

% benchmark rotation: since we need an arbitrary rotation, we just use
% cholesky, but any other will do. Ordering does not matter here.

bench_rot = chol(Sigma_u,'lower');

% Wold IRFs

IRF_Wold = zeros(n_y,n_y,IRF_hor); % row is variable, column is shock
IRF_Wold(:,:,1) = eye(n_y);

for l = 1:IRF_hor
    
    if l < IRF_hor
        for j=1:min(l,n_lags)
            IRF_Wold(:,:,l+1) = IRF_Wold(:,:,l+1) + B(1+(j-1)*n_y:j*n_y,:)'*IRF_Wold(:,:,l-j+1);
        end
    end
    
end

W = bench_rot;

% get IRFs

IRF_OLS = NaN(n_y,n_y,IRF_hor);
for i_hor = 1:IRF_hor
    IRF_OLS(:,:,i_hor) = IRF_Wold(:,:,i_hor) * W;
end

% collect results

IS.Theta_OLS = squeeze(IRF_OLS);

%----------------------------------------------------------------
% Wold IRFs
%----------------------------------------------------------------

IS.Theta     = NaN(n_y,n_y,IRF_hor,n_draws);
IS.Theta_med = NaN(n_y,n_y,IRF_hor);
IS.Theta_lb  = NaN(n_y,n_y,IRF_hor);
IS.Theta_ub  = NaN(n_y,n_y,IRF_hor);

% do the same for each posterior draw.

for i_draw = 1:n_draws

% extract VAR inputs
    
Sigma_u   = Sigma_draws(:,:,i_draw);
B         = B_draws(:,:,i_draw);

% benchmark rotation

bench_rot = chol(Sigma_u,'lower');

% Wold IRFs

IRF_Wold = zeros(n_y,n_y,IRF_hor); % row is variable, column is shock
IRF_Wold(:,:,1) = eye(n_y);

for l = 1:IRF_hor
    
    if l < IRF_hor
        for j=1:min(l,n_lags)
            IRF_Wold(:,:,l+1) = IRF_Wold(:,:,l+1) + B(1+(j-1)*n_y:j*n_y,:)'*IRF_Wold(:,:,l-j+1);
        end
    end
    
end

W = bench_rot;

% get IRFs

IRF_idraw = NaN(n_y,n_y,IRF_hor);
for i_hor = 1:IRF_hor
    IRF_idraw(:,:,i_hor) = IRF_Wold(:,:,i_hor) * W;
end

% collect results

IS.Theta(:,:,:,i_draw) = squeeze(IRF_idraw);

end

% compute percentiles

for ii=1:n_y
    for jj=1:n_y
        for kk = 1:IRF_hor
            IS.Theta_med(ii,jj,kk) = quantile(IS.Theta(ii,jj,kk,:),0.5);
            IS.Theta_lb(ii,jj,kk) = quantile(IS.Theta(ii,jj,kk,:),0.16);
            IS.Theta_ub(ii,jj,kk) = quantile(IS.Theta(ii,jj,kk,:),0.84);
        end
    end
end

% re-shuffle ordering for plots

IS.Theta_med = permute(IS.Theta_med,[3 1 2]); % order: horizon, variable, shock
IS.Theta_lb  = permute(IS.Theta_lb,[3 1 2]);
IS.Theta_ub  = permute(IS.Theta_ub,[3 1 2]);

%% COMPUTE SECOND-MOMENT RESULTS

% VMA-implied variance-covariance matrix

IS.cov = zeros(n_y,n_y);
for i_hor = 1:IRF_hor
    IS.cov = IS.cov + IS.Theta_OLS(:,:,i_hor) * IS.Theta_OLS(:,:,i_hor)';
end

% correlations

IS.corr = zeros(n_y,n_y);
for i_y = 1:n_y
    for ii_y = 1:n_y
        IS.corr(i_y,ii_y) = IS.cov(i_y,ii_y)/sqrt(IS.cov(i_y,i_y) * IS.cov(ii_y,ii_y));
    end
end

% frequency bands

omega_1 = (2*pi)/32;
omega_2 = (2*pi)/6;

IS.freq_var = diag(freq_var_fn(IS.Theta_OLS,omega_1,omega_2));

%% PLOT IRFs

%----------------------------------------------------------------
% General Settings
%----------------------------------------------------------------

if plot_wold_irfs == 1

% plot horizon

IRF_hor_plot = 28;

% color settings

settings.colors.black  = [0 0 0];
settings.colors.grey   = [230/255 230/255 230/255];
settings.colors.orange = [204/255 102/255 0/255];
settings.colors.green = [37/255 152/255 14/255];
settings.colors.navyblue = [0/255 0/255 50/255];
settings.colors.purple = [160/255 32/255 240/255];
settings.colors.list = [settings.colors.black;settings.colors.grey;settings.colors.orange];

% figure spacing

plotwidth = 0.15;
gapsize = 0.04;
gapsize_edges = (1-5*plotwidth-4*gapsize)/2;
left_pos = [gapsize_edges, gapsize_edges + gapsize + plotwidth, gapsize_edges + 2 * gapsize + 2 * plotwidth, ...
    gapsize_edges + 3 * gapsize + 3 * plotwidth, gapsize_edges + 4 * gapsize + 4 * plotwidth];

%----------------------------------------------------------------
% Plot Wold IRFs
%----------------------------------------------------------------

for ii_y = 1:n_y % shock

    figure(ii_y)

    for i_y = 1:n_y % variable
        subplot(2,5,i_y)
        pos = get(gca, 'Position');
        set(gca,'FontSize',15)
        set(gca,'Position', pos)
        set(gca,'TickLabelInterpreter','latex')
        hold on
        jbfill(0:1:IRF_hor_plot-1,IS.Theta_lb(1:IRF_hor_plot,i_y,ii_y)',...
            IS.Theta_ub(1:IRF_hor_plot,i_y,ii_y)',settings.colors.grey,settings.colors.grey,0,1);
        plot(0:1:IRF_hor_plot-1,IS.Theta_med(1:IRF_hor_plot,i_y,ii_y)','Color',settings.colors.navyblue,'LineWidth',4)
        hold on
        xlim([0 IRF_hor_plot-1])
        set(gcf,'color','w')
        title(series_names(i_y),'interpreter','latex','fontsize',24)
        xlabel('Horizon','interpreter','latex','FontSize',20)
        if i_y == 1 || i_y == 6
            ylabel('\% Deviation','interpreter','latex','FontSize',20)
        end
        grid on
        hold off
    end
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) 3*pos(3) 2*pos(4)]);
    set(gcf, 'PaperPositionMode', 'auto');

end
end
%% SAVE RESULTS
if save_results == 1

IS_wold = IS;

cd([path session task '/_results']);

save wold_results IS_wold series_names

cd([path session task]);
end