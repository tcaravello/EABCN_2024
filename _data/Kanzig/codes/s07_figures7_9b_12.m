%% Replication files for "The macroeconomic effects of oil supply news"
% Step 7: This file creates figures 7, 9b and 12 in the paper

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
dataFrequency = 'Q';

% Estimation sample
smplStart_base = '1974Q1'; 
smplEnd   = '2017Q4'; 

% Instrument sample
% range has to be contained in estimation sample
smplStartProxy_base = '1975Q1'; 
smplEndProxy   = '2017Q4'; 

% VAR specifics
p          = 4;          % Lag order
horizon    = 20;         % Horizon for IRFs
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
includeBase = false;     % Include baseline response in plots

verbo = false;           % Control verbosity
saveFigs   = true;       % Save figures to disk


%% Read in data
% data for baseline model
load('../data/dataBaseQ')
% data: transformed endogenous variables 
% dataExo: exogenous variables (e.g. constant, trend)
% sampleDates: sample dates (string format)
% sampleDatesNum: sample dates (numeric format, e.g. 2000 = 2000Q1)
% varNames: labels of variables

% number of variables in baseline VAR
nvar_base = size(data,2);  

% names for paper
varNames_paper_base = {'Real oil price','World oil production','World oil inventories','World industrial production','U.S. GDP','U.S. CPI'};

% select sample
smplStartInd = find(strcmp(sampleDates,smplStart_base));
smplEndInd   = find(strcmp(sampleDates,smplEnd));

data = data(smplStartInd:smplEndInd,:);
dataExo = dataExo(smplStartInd:smplEndInd,:);

data_base = data;
dataExo_base = dataExo;

% data to extend model
load('../data/dataExtQ')

% number of extended VARs
n_extended = size(data_extended,2);

data_extended = data_extended(smplStartInd:smplEndInd,:);
sampleDatesRaw = sampleDates(smplStartInd:smplEndInd,:);
sampleDatesRawNum = sampleDatesNum(smplStartInd:smplEndInd,:);

% read in results from monthly VAR
load('../results/IRFsext')
IRFs_extendedM = IRFs_extended;
IRFsupper_extendedM = IRFsupper_extended;
IRFslower_extendedM = IRFslower_extended;
IRFsupper2_extendedM = IRFsupper2_extended;
IRFslower2_extendedM = IRFslower2_extended;

%% Augment quarterly VAR by one variable at a time and map out IRFs

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

%% Run baseline VAR at quarterly frequency
rng default

smplStartProxy = smplStartProxy_base; 
data = data_base;
nvar = nvar_base;

loadProxy;
    
runProxyVAR; 

%% Plotting
time = (0:horizon)';

% Inflation expectations
figure('Position',[100 100 645 255],'PaperPositionMode','Auto','DefaultAxesFontSize',13);
signIRFs = 1;
subplot(1,2,1);
j = 2;
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
title('Households (Michigan survey)') 
xlabel('Quarters');
ylabel('ppt');
xlim([0,horizon]);
ylim([-0.1 0.4])
subplot(1,2,2)
j = 1;
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
title('Forecasters (SPF)')
xlabel('Quarters');
ylabel('ppt');
xlim([0,horizon]);
ylim([-0.1 0.4])
tightfig;
if saveFigs
    print('-dpdf', gcf, '../results/figure7');  
end

% activity indicators
figure('Position',[100 100 1000 255],'PaperPositionMode','Auto','DefaultAxesFontSize',13);
signIRFs = 1;
subplot(1,3,1);
hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper_proxy(1,5); IRFslower_proxy(1:end,5); flipud([IRFsupper_proxy(1:end,5); IRFslower_proxy(end,5)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.2);
set(hh,'edgecolor','none');
hold on;
hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper2_proxy(1,5); IRFslower2_proxy(1:end,5); flipud([IRFsupper2_proxy(1:end,5); IRFslower2_proxy(end,5)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.4);
set(hh,'edgecolor','none');
plot(time, signIRFs*IRFs_proxy(:,5),'k', 'Linewidth', 1.5); hold on;
if ~ismember(0,get(gca,'ylim'))
    line(get(gca,'xlim'),[0 0],'Color','k')
end
grid on ;hold off;
title('Real GDP') 
xlabel('Quarters');
ylabel('\%');
xlim([0,horizon]);
xticks([0:5:horizon]);
subplot(1,3,3)
j = 3;
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
title(varNames_paper_extended{3}) 
xlabel('Quarters');
ylabel('\%');
xlim([0,horizon]);
xticks([0:5:horizon]);
subplot(1,3,2)
j = 4;
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
title(varNames_paper_extended{4}) 
xlabel('Quarters');
ylabel('\%');
xlim([0,horizon]);
xticks([0:5:horizon]);
tightfig;
if saveFigs
    print('-dpdf', gcf, '../results/figure9b');  
end

% trade
timeM = (0:50)';
figure('Position',[100 100 645 255],'PaperPositionMode','Auto','DefaultAxesFontSize',13);
signIRFs = 1;
j = 33;
subplot(1,2,1);
hh=fill([timeM(1); timeM(1:end); flipud([timeM(1:end); timeM(end)])],[IRFsupper_extendedM(1,end,j); IRFslower_extendedM(1:end,end,j); flipud([IRFsupper_extendedM(1:end,end,j); IRFslower_extendedM(end,end,j)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.2);
set(hh,'edgecolor','none');
hold on;
hh=fill([timeM(1); timeM(1:end); flipud([timeM(1:end); timeM(end)])],[IRFsupper2_extendedM(1,end,j); IRFslower2_extendedM(1:end,end,j); flipud([IRFsupper2_extendedM(1:end,end,j); IRFslower2_extendedM(end,end,j)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.4);
set(hh,'edgecolor','none');
plot(timeM, signIRFs*IRFs_extendedM(:,end,j),'k', 'Linewidth', 1.5); hold on;
if ~ismember(0,get(gca,'ylim'))
    line(get(gca,'xlim'),[0 0],'Color','k')
end
grid on ;hold off;
title('Terms of trade') 
xlabel('Months');
ylabel('\%');
xlim([0,50]);
j = 5;
subplot(1,2,2)
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
title('Trade balance')
xlabel('Quarters');
ylabel('ppt');
xlim([0,20]);
yticks(-0.1:0.02:0.02)
tightfig;
if saveFigs
    print('-dpdf', gcf, '../results/figure12');  
end

% baseline responses for appendix
figure('Position',[100 100 1000 600],'PaperPositionMode','Auto','DefaultAxesFontSize',13);
signIRFs = 1;
for j=1:nvar %variable
    h(j) = subplot(ceil(nvar/3),3,j); %note the index

    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper_proxy(1,j); IRFslower_proxy(1:end,j); flipud([IRFsupper_proxy(1:end,j); IRFslower_proxy(end,j)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.2);
    set(hh,'edgecolor','none'); 

    hold on;
    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper2_proxy(1,j); IRFslower2_proxy(1:end,j); flipud([IRFsupper2_proxy(1:end,j); IRFslower2_proxy(end,j)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.4);
    set(hh,'edgecolor','none');

    if includeBase
        load('../results/IRFsbench')
        p2=plot(time, signIRFs*IRFs_base(:,j), 'Linewidth', 2,'Color',[0.8500, 0.3250, 0.0980],'LineStyle',':');
    end

    p1=plot(time, signIRFs*IRFs_proxy(:,j),'k', 'Linewidth', 1.5); 
    if ~ismember(0,get(gca,'ylim'))
        line(get(gca,'xlim'),[0 0],'Color','k')
    end
    box on
    grid on ;hold off;
    title(varNames_paper_base{j}) 
    ylabel('\%');
    xlim([0,horizon]);
    xlabel('Quarters');
    xticks(0:5:horizon)
    if j==1 && includeBase
        legend([p1 p2],{legendText,'Baseline'})
    end
end
pause(0.001)
h=axes('Position',[0.25,0,.5,.5],'Xlim',[0 1],'Ylim',[0 1]);
set(h,'Visible','off');
if k==1
    string_1stage = ['First stage regression: F: ',num2str(olsEst.F,' %2.2f'),', robust F: ',num2str(olsEst.Frobust,' %2.2f'),', $R^2$: ',num2str(olsEst.R2*100,' %1.2f'),'\%, Adjusted $R^2$: ',num2str(olsEst.R2adj*100,' %1.2f'),'\%'];
    text('Position',[-0.16 -0.002],'string',string_1stage,'FontSize',14);
end
tightfig;
if saveFigs
    print('-dpdf', gcf, '../results/appendix/figurea31');  
end

