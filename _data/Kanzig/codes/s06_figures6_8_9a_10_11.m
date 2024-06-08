%% Replication files for "The macroeconomic effects of oil supply news"
% Step 6: This file creates figures 6, 8, 9a, 10, and 11 in the paper

% Diego R. Känzig
% LBS, September 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

% add tools directory
addpath(genpath('auxfiles'))

% initialize random number generator
rng default

% Set text interpreter to latex
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% change specific options in this file
% data frequency
dataFrequency = 'M';

% Estimation sample
smplStart_base = '1974M01'; 
smplEnd   = '2017M12'; 

% Instrument sample
% range has to be contained in estimation sample
smplStartProxy_base = '1975M01'; 
smplEndProxy   = '2017M12'; 

% VAR specifics
p          = 12;         % Lag order
horizon    = 50;         % Horizon for IRFs
shockType  = 'custom';   % one standard deviation 'sd' or 'custom'
shockSize  = 10;         % if custom, specify shock size here
alpha      = 0.1;        % Significance level for bands (alpha=0.1 => 90% CIs (two SD))
alpha2     = 0.32;
nsim       = 10000;      % number of simulations in bootstrap
bootType   = 'mbb1block';% Moving block bootstrap

% proxy
ncontract = 14;          % Principal component of contracts spanning first year of term structure

% switches
compFEVDs  = false;      % Compute variance decomposition
includeBase = false;

verbo = false;           % Control verbosity
saveFigs   = true;       % Save figures to disk


%% Read in data
% data for baseline model
load('../data/dataBaseM')
% data: transformed endogenous variables 
% dataExo: exogenous variables (e.g. constant, trend)
% sampleDates: sample dates (string format)
% sampleDatesNum: sample dates (numeric format, e.g. 2000 = 2000M1)
% varNames: labels of variables

% number of variables in baseline VAR
nvar_base = size(data,2);  

% names for paper
varNames_paper_base = {'Real oil price','World oil production','World oil inventories','World industrial production','U.S. industrial production','U.S. CPI'};

% select sample
smplStartInd = find(strcmp(sampleDates,smplStart_base));
smplEndInd   = find(strcmp(sampleDates,smplEnd));

data = data(smplStartInd:smplEndInd,:);
dataExo = dataExo(smplStartInd:smplEndInd,:);

data_base = data;
dataExo_base = dataExo;

% data to extend model
load('../data/dataExtM')
load('../data/appendix/dataExtStocksM')
data_extended = [data_extended data_extended_stocks];
varNames_paper_extended = [varNames_paper_extended varNames_paper_extended_stocks];


% number of extended VARs
n_extended = size(data_extended,2);

data_extended = data_extended(smplStartInd:smplEndInd,:);
sampleDatesRaw = sampleDates(smplStartInd:smplEndInd,:);
sampleDatesRawNum = sampleDatesNum(smplStartInd:smplEndInd,:);

%% Augment baseline VAR by one variable at a time and map out IRFs

IRFs_extended = nan(horizon+1,nvar_base+1,n_extended);
IRFsupper_extended = nan(horizon+1,nvar_base+1,n_extended);
IRFslower_extended = nan(horizon+1,nvar_base+1,n_extended);
IRFsupper2_extended = nan(horizon+1,nvar_base+1,n_extended);
IRFslower2_extended = nan(horizon+1,nvar_base+1,n_extended);

for id = 1:n_extended 
    
    disp(id)
    
    % correct start date if necessary (if additional series only spans a
    % shorter sample)
    startInd_extended = find(~isnan(data_extended(:,id)),1);
    
    if startInd_extended>find(strcmp(sampleDatesRaw,smplStartProxy_base))-p 
        smplStartProxy = sampleDatesRaw(startInd_extended+p); 
    else
        smplStartProxy = smplStartProxy_base; 
    end
    
    %% Data manipulations
    data = [data_base data_extended(:,id)];     
    nvar = size(data,2);
    
    varNames_paper = [varNames_paper_base varNames_paper_extended(id)];
    varNames = varNames_paper;
    
    data = data(startInd_extended:smplEndInd,:);
    dataExo = dataExo_base(startInd_extended:smplEndInd,:);
    sampleDates = sampleDatesRaw(startInd_extended:smplEndInd,:);
    sampleDatesNum = sampleDatesRawNum(startInd_extended:smplEndInd,:);
    
    
    %% External instruments VAR

    loadProxy;
    
    runProxyVAR; 

    % save outputs
    IRFs_extended(:,:,id) = IRFs_proxy;   % unit shock IRFs to pos oil price shock
    IRFsupper_extended(:,:,id) = IRFsupper_proxy;
    IRFslower_extended(:,:,id) = IRFslower_proxy;
    IRFsupper2_extended(:,:,id) = IRFsupper2_proxy;
    IRFslower2_extended(:,:,id) = IRFslower2_proxy;
end

% save results
save('../results/IRFsext','IRFs_extended','IRFsupper_extended','IRFslower_extended','IRFsupper2_extended','IRFslower2_extended')

%% Generate figures
time = (0:horizon)';

% get baseline responses
load('../results/IRFsbench')

% expectations
j = 12;
signIRFs = 1;
figure('Position',[100 100 645 255],'PaperPositionMode','Auto','DefaultAxesFontSize',13);
subplot(1,2,1); 
hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper_extended(1,end,j); IRFslower_extended(1:end,end,j); flipud([IRFsupper_extended(1:end,end,j); IRFslower_extended(end,end,j)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.2);
set(hh,'edgecolor','none');
hold on;
hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper2_extended(1,end,j); IRFslower2_extended(1:end,end,j); flipud([IRFsupper2_extended(1:end,end,j); IRFslower2_extended(end,end,j)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.4);
set(hh,'edgecolor','none');
plot(time, signIRFs*IRFs_extended(:,end,j),'k', 'Linewidth', 1.5); hold on;
ylabel('\%');
title('Oil price expectations') 
if ~ismember(0,get(gca,'ylim'))
    line(get(gca,'xlim'),[0 0],'Color','k')
end
grid on ;hold off;
xlabel('Months');
xlim([0,horizon]);
xticks([0:10:horizon]);
ax = gca;
ax.YRuler.Exponent = 0;
j = j+1;
subplot(1,2,2); 
hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper_extended(1,end,j); IRFslower_extended(1:end,end,j); flipud([IRFsupper_extended(1:end,end,j); IRFslower_extended(end,end,j)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.2);
set(hh,'edgecolor','none');
hold on;
hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper2_extended(1,end,j); IRFslower2_extended(1:end,end,j); flipud([IRFsupper2_extended(1:end,end,j); IRFslower2_extended(end,end,j)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.4);
set(hh,'edgecolor','none');
plot(time, signIRFs*IRFs_extended(:,end,j),'k', 'Linewidth', 1.5); hold on;
ylabel('ppt');
title('Inflation expectations') 
if ~ismember(0,get(gca,'ylim'))
    line(get(gca,'xlim'),[0 0],'Color','k')
end
grid on ;hold off;
xlabel('Months');
xlim([0,horizon]);
xticks([0:10:horizon]);
ax = gca;
ax.YRuler.Exponent = 0;
tightfig;
if saveFigs
    print('-dpdf', gcf, '../results/figure6a');  
end

% uncertainty
j = j+1;
signIRFs = 1;
figure('Position',[100 100 645 255],'PaperPositionMode','Auto','DefaultAxesFontSize',13);
subplot(1,2,1); 
hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper_extended(1,end,j); IRFslower_extended(1:end,end,j); flipud([IRFsupper_extended(1:end,end,j); IRFslower_extended(end,end,j)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.2);
set(hh,'edgecolor','none');
hold on;
hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper2_extended(1,end,j); IRFslower2_extended(1:end,end,j); flipud([IRFsupper2_extended(1:end,end,j); IRFslower2_extended(end,end,j)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.4);
set(hh,'edgecolor','none');
plot(time, signIRFs*IRFs_extended(:,end,j),'k', 'Linewidth', 1.5); hold on;
ylabel('\%');
title('Financial uncertainty') 
if ~ismember(0,get(gca,'ylim'))
    line(get(gca,'xlim'),[0 0],'Color','k')
end
grid on ;hold off;
xlabel('Months');
xlim([0,horizon]);
xticks([0:10:horizon]);
ax = gca;
ax.YRuler.Exponent = 0;
j = j+1;
subplot(1,2,2); 
hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper_extended(1,end,j); IRFslower_extended(1:end,end,j); flipud([IRFsupper_extended(1:end,end,j); IRFslower_extended(end,end,j)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.2);
set(hh,'edgecolor','none');
hold on;
hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper2_extended(1,end,j); IRFslower2_extended(1:end,end,j); flipud([IRFsupper2_extended(1:end,end,j); IRFslower2_extended(end,end,j)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.4);
set(hh,'edgecolor','none');
plot(time, signIRFs*IRFs_extended(:,end,j),'k', 'Linewidth', 1.5); hold on;
ylabel('\%');
title('Geopolitical risk') 
if ~ismember(0,get(gca,'ylim'))
    line(get(gca,'xlim'),[0 0],'Color','k')
end
grid on ;hold off;
xlabel('Months');
xlim([0,horizon]);
xticks([0:10:horizon]);
ax = gca;
ax.YRuler.Exponent = 0;
tightfig;
if saveFigs
    print('-dpdf', gcf, '../results/figure6b');  
end

% consumer prices
% headline
names={'Headline','Energy','Non-durables','Core','Durables','Services'};
figure('Position',[100 100 1000 600],'PaperPositionMode','Auto','DefaultAxesFontSize',13);
signIRFs = 1;
subplot(2,3,1)
hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper_base(1,6); IRFslower_base(1:end,6); flipud([IRFsupper_base(1:end,6); IRFslower_base(end,6)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.2);
set(hh,'edgecolor','none');
hold on;
hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper2_base(1,6); IRFslower2_base(1:end,6); flipud([IRFsupper2_base(1:end,6); IRFslower2_base(end,6)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.4);
set(hh,'edgecolor','none');
plot(time, signIRFs*IRFs_base(:,6),'k', 'Linewidth', 1.5); hold on;
if ~ismember(0,get(gca,'ylim'))
    line(get(gca,'xlim'),[0 0],'Color','k')
end
grid on ;hold off;
title(names{1});
xlabel('Months');
ylabel('\%');
xlim([0,horizon]);
xticks([0:10:horizon]);
% other prices
inds = [2 3 1 4 5];
ii = 2;
for j=inds
    subplot(2,3,ii); 
    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper_extended(1,end,j); IRFslower_extended(1:end,end,j); flipud([IRFsupper_extended(1:end,end,j); IRFslower_extended(end,end,j)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.2);
    set(hh,'edgecolor','none');
    hold on;
    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper2_extended(1,end,j); IRFslower2_extended(1:end,end,j); flipud([IRFsupper2_extended(1:end,end,j); IRFslower2_extended(end,end,j)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.4);
    set(hh,'edgecolor','none');
    plot(time, signIRFs*IRFs_extended(:,end,j),'k', 'Linewidth', 1.5); hold on;
    if ~ismember(0,get(gca,'ylim'))
        line(get(gca,'xlim'),[0 0],'Color','k')
    end
    grid on ;hold off;
    title(names{ii}); 
    xlabel('Months');
    ylabel('\%');
    xlim([0,horizon]);
    xticks([0:10:horizon]);
    if ii==6
        ylim([-0.4 0.6])
    end
    ii = ii+1;
end
tightfig;
if saveFigs
    print('-dpdf', gcf, '../results/figure8');  
end

% activity
figure('Position',[100 100 1000 255],'PaperPositionMode','Auto','DefaultAxesFontSize',13);
signIRFs = 1;
subplot(1,3,1);
hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper_base(1,5); IRFslower_base(1:end,5); flipud([IRFsupper_base(1:end,5); IRFslower_base(end,5)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.2);
set(hh,'edgecolor','none');
hold on;
hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper2_base(1,5); IRFslower2_base(1:end,5); flipud([IRFsupper2_base(1:end,5); IRFslower2_base(end,5)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.4);
set(hh,'edgecolor','none');
plot(time, signIRFs*IRFs_base(:,5),'k', 'Linewidth', 1.5); hold on;
if ~ismember(0,get(gca,'ylim'))
    line(get(gca,'xlim'),[0 0],'Color','k')
end
grid on ;hold off;
title('Industrial production')
xlabel('Months');
ylabel('\%');
xlim([0,horizon]);
xticks([0:10:horizon]);
ax = gca;
ax.YRuler.Exponent = 0;
subplot(1,3,2); 
j = 6;
hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper_extended(1,end,j); IRFslower_extended(1:end,end,j); flipud([IRFsupper_extended(1:end,end,j); IRFslower_extended(end,end,j)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.2);
set(hh,'edgecolor','none');
hold on;
hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper2_extended(1,end,j); IRFslower2_extended(1:end,end,j); flipud([IRFsupper2_extended(1:end,end,j); IRFslower2_extended(end,end,j)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.4);
set(hh,'edgecolor','none');
plot(time, signIRFs*IRFs_extended(:,end,j),'k', 'Linewidth', 1.5); hold on;
if ~ismember(0,get(gca,'ylim'))
    line(get(gca,'xlim'),[0 0],'Color','k')
end
grid on ;hold off;
title(varNames_paper_extended{j}) 
xlabel('Months');
ylabel('ppt');
xlim([0,horizon]);
xticks([0:10:horizon]);
ax = gca;
ax.YRuler.Exponent = 0;
subplot(1,3,3); 
j = j+1;
hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper_extended(1,end,j); IRFslower_extended(1:end,end,j); flipud([IRFsupper_extended(1:end,end,j); IRFslower_extended(end,end,j)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.2);
set(hh,'edgecolor','none');
hold on;
hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper2_extended(1,end,j); IRFslower2_extended(1:end,end,j); flipud([IRFsupper2_extended(1:end,end,j); IRFslower2_extended(end,end,j)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.4);
set(hh,'edgecolor','none');
plot(time, signIRFs*IRFs_extended(:,end,j),'k', 'Linewidth', 1.5); hold on;
if ~ismember(0,get(gca,'ylim'))
    line(get(gca,'xlim'),[0 0],'Color','k')
end
grid on ;hold off;
title(varNames_paper_extended{j}) 
xlabel('Months');
ylabel('\%');
xlim([0,horizon]);
xticks([0:10:horizon]);
yticks([-0.6:0.2:0.2]);
ylim([-0.75 0.25])
ax = gca;
ax.YRuler.Exponent = 0;
tightfig;
if saveFigs
    print('-dpdf', gcf, '../results/figure9a');  
end

% monetary policy
figure('Position',[100 100 1000 255],'PaperPositionMode','Auto','DefaultAxesFontSize',13);
ii = 1;
for j=16:18
    signIRFs = 1;
    subplot(1,3,ii); 
    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper_extended(1,end,j); IRFslower_extended(1:end,end,j); flipud([IRFsupper_extended(1:end,end,j); IRFslower_extended(end,end,j)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.2);
    set(hh,'edgecolor','none');
    hold on;
    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper2_extended(1,end,j); IRFslower2_extended(1:end,end,j); flipud([IRFsupper2_extended(1:end,end,j); IRFslower2_extended(end,end,j)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.4);
    set(hh,'edgecolor','none');
    plot(time, signIRFs*IRFs_extended(:,end,j),'k', 'Linewidth', 1.5); hold on;
    if ~ismember(0,get(gca,'ylim'))
        line(get(gca,'xlim'),[0 0],'Color','k')
    end
    grid on ;hold off;
    title(varNames_paper_extended{j}) 
    xlabel('Months');
    ylabel('\%');
    if ii<3
        ylabel('ppt')
    end
    xlim([0,horizon]);
    xticks([0:10:horizon]);
    ax = gca;
    ax.YRuler.Exponent = 0;
    ii = ii+1;
end
tightfig;
if saveFigs
    print('-dpdf', gcf, '../results/figure10');  
end

% NEER
figure('Position',[100 100 645 265],'PaperPositionMode','Auto','DefaultAxesFontSize',13);
signIRFs = 1;
subplot(1,2,1);
j = 20;
hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper_extended(1,end,j); IRFslower_extended(1:end,end,j); flipud([IRFsupper_extended(1:end,end,j); IRFslower_extended(end,end,j)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.2);
set(hh,'edgecolor','none');
hold on;
hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper2_extended(1,end,j); IRFslower2_extended(1:end,end,j); flipud([IRFsupper2_extended(1:end,end,j); IRFslower2_extended(end,end,j)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.4);
set(hh,'edgecolor','none');
plot(time, signIRFs*IRFs_extended(:,end,j),'k', 'Linewidth', 1.5); hold on;
ylabel('\%');

title('Narrow index') 
if ~ismember(0,get(gca,'ylim'))
    line(get(gca,'xlim'),[0 0],'Color','k')
end
grid on ;hold off;
xlabel('Months');
xlim([0,horizon]);
xticks([0:10:horizon]);
ylim([-3 2])
ax = gca;
ax.YRuler.Exponent = 0;
subplot(1,2,2); 
j = j-1;
hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper_extended(1,end,j); IRFslower_extended(1:end,end,j); flipud([IRFsupper_extended(1:end,end,j); IRFslower_extended(end,end,j)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.2);
set(hh,'edgecolor','none');
hold on;
hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper2_extended(1,end,j); IRFslower2_extended(1:end,end,j); flipud([IRFsupper2_extended(1:end,end,j); IRFslower2_extended(end,end,j)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.4);
set(hh,'edgecolor','none');
plot(time, signIRFs*IRFs_extended(:,end,j),'k', 'Linewidth', 1.5); hold on;
ylabel('\%');
title('Broad index') 
if ~ismember(0,get(gca,'ylim'))
    line(get(gca,'xlim'),[0 0],'Color','k')
end
grid on ;hold off;
xlabel('Months');
xlim([0,horizon]);
xticks([0:10:horizon]);
ax = gca;
ax.YRuler.Exponent = 0;
ylim([-3 0.5])
tightfig;
if saveFigs
    print('-dpdf', gcf, '../results/figure11a');  
end

% Billateral exchange rates (nominal)
netExpClassification = [2 1 1 2 1 2 2 2 2 1 1 1]; % 1 net importer, 2: net exporter
C = linspecer(length(netExpClassification)); 
figure('Position',[100 100 645 460],'PaperPositionMode','Auto','DefaultAxesFontSize',13);
signIRFs = 1;
ii = 1;
ind = 21;
for j=ind-1+find(netExpClassification==1)
    subplot(1,2,1); 
    plot(time, signIRFs*IRFs_extended(:,end,j), 'Linewidth', 1.5,'color',C(ii,:)); hold on;
    ii = ii+1;
end
ylim([-9 2])
line(get(gca,'xlim'),[0 0],'Color','k')
grid on 
title('Oil importer currencies') 
xlabel('Months');
ylabel('\%');
xlim([0,horizon]);
xticks([0:10:horizon]);
ax = gca;
ax.YRuler.Exponent = 0;
legend(varNames_paper_extended{ind-1+find(netExpClassification==1)},'Location', 'southoutside','fontsize',10)

for j=ind-1+find(netExpClassification==2)
    subplot(1,2,2); 
    plot(time, signIRFs*IRFs_extended(:,end,j), 'Linewidth', 1.5,'color',C(ii,:)); hold on;
    ii = ii+1;
end
ylim([-9 2])
line(get(gca,'xlim'),[0 0],'Color','k')
grid on 
title('Oil exporter currencies') 
xlabel('Months');
ylabel('\%');
xlim([0,horizon]);
xticks([0:10:horizon]);
ax = gca;
ax.YRuler.Exponent = 0;
legend(varNames_paper_extended{ind-1+find(netExpClassification==2)},'Location', 'southoutside','fontsize',10)
tightfig;
if saveFigs
    saveas(gcf,'../results/figure11b','epsc'); 
end


%% Figures for appendix

% expenditures (components)
names={'Energy','Durables','Non-durables','Services'};
figure('Position',[100 100 645 600],'PaperPositionMode','Auto','DefaultAxesFontSize',13);
ii = 1;
signIRFs = 1;
for j=[8 10 9 11]
    subplot(2,2,ii); %note the index
    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper_extended(1,end,j); IRFslower_extended(1:end,end,j); flipud([IRFsupper_extended(1:end,end,j); IRFslower_extended(end,end,j)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.2);
    set(hh,'edgecolor','none');
    hold on;
    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper2_extended(1,end,j); IRFslower2_extended(1:end,end,j); flipud([IRFsupper2_extended(1:end,end,j); IRFslower2_extended(end,end,j)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.4);
    set(hh,'edgecolor','none');
    plot(time, signIRFs*IRFs_extended(:,end,j),'k', 'Linewidth', 1.5); 
    if ~ismember(0,get(gca,'ylim'))
        line(get(gca,'xlim'),[0 0],'Color','k')
    end
    grid on ;hold off;
    title(names{ii}); 
    xlabel('Months');
    ylabel('\%');

    xlim([0,horizon]);
    xticks([0:10:horizon]);
    ii = ii+1;
end
tightfig;
if saveFigs
    print('-dpdf', gcf, '../results/appendix/figurea7');  
end

% stock price indices
figure('Position',[100 100 1000 600],'PaperPositionMode','Auto','DefaultAxesFontSize',13);
ii = 1;
signIRFs = 1;
for j=34:39
    subplot(2,3,ii); 
    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper_extended(1,end,j); IRFslower_extended(1:end,end,j); flipud([IRFsupper_extended(1:end,end,j); IRFslower_extended(end,end,j)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.2);
    set(hh,'edgecolor','none');
    hold on;
    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper2_extended(1,end,j); IRFslower2_extended(1:end,end,j); flipud([IRFsupper2_extended(1:end,end,j); IRFslower2_extended(end,end,j)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.4);
    set(hh,'edgecolor','none');
    plot(time, signIRFs*IRFs_extended(:,end,j),'k', 'Linewidth', 1.5); hold on;
    if ~ismember(0,get(gca,'ylim'))
        line(get(gca,'xlim'),[0 0],'Color','k')
    end
    grid on ;hold off;
    title(varNames_paper_extended{j}) 
    xlabel('Months');
    ylabel('\%');
    xlim([0,horizon]);
    xticks([0:10:horizon]);
    ii = ii+1;
end
tightfig;
if saveFigs
    print('-dpdf', gcf, '../results/appendix/figurea8');  
end


% Robustness
varNames_paper = {'Real oil price','World oil production','World oil inventories','World industrial production','U.S. industrial production','U.S. CPI'};
figure('Position',[100 100 1000 600],'PaperPositionMode','Auto','DefaultAxesFontSize',13);
signIRFs = 1;
for ii = 1:n_extended
    % exclude series with shorter sample to allow for better comparison
    if strcmp(varNames_paper_extended(ii),'Oil price expectations') || strcmp(varNames_paper_extended(ii),'Inflation expectations (Michigan)') || strcmp(varNames_paper_extended(ii),'Geopolitical risk') || strcmp(varNames_paper_extended(ii),'Russia')
        continue
    end
    for j=1:nvar_base 
        subplot(2,3,j); 
        hold on
        plot(time, signIRFs*IRFs_extended(:,j,ii),'k');  
        if ~ismember(0,get(gca,'ylim'))
            line(get(gca,'xlim'),[0 0],'Color','k')
        end
        grid on; box on
        title(varNames_paper{j}) 
        xlabel('Months');
        ylabel('\%');

        xlim([0,horizon]);
        xticks([0:10:horizon]);
    end
end
tightfig;
if saveFigs
    print('-dpdf', gcf, '../results/appendix/figurea15'); 
end