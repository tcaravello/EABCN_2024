%% Replication files for "The macroeconomic effects of oil supply news"
% This file creates figure A.9 in the appendix

% Diego R. Känzig
% LBS, September 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all
clc

% add tools directory
addpath(genpath('../auxfiles'))

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
smplStart = '1974M01'; 
smplEnd   = '2017M12'; 

% Instrument sample
% range has to be contained in estimation sample
smplStartProxy = '1975M01'; 
smplEndProxy   = '2017M12'; 

% VAR specifics
p          = 12;       
horizon    = 50;
shockType  = 'custom';   % one standard deviation 'sd' or 'custom'
shockSize  = 10;         % if custom, specify shock size here
alpha      = 0.1;        % Significance level for bands (alpha=0.1 => 90% CIs (two SD))
alpha2     = 0.32;
nsim       = 10000;      % number of simulations in bootstrap
bootType   = 'mbb1block'; 

% proxy
ncontract = 14; 

% switches
compFEVDs  = false;
includeBase = true;

verbo = false;
saveFigs   = true;


%% Read in data
load('../../data/dataBaseM')
% data: transformed endogenous variables 
% dataExo: exogenous variables (e.g. constant, trend)
% sampleDates: sample dates (string format)
% sampleDatesNum: sample dates (numeric format, e.g. 2000 = 2000M1)
% varNames: labels of variables

% augment by FFR
load('../../data/dataExtM')

data = [data data_extended(:,16)];

% number of variables in VAR
nvar = size(data,2);  

% names for paper
varNames_paper = {'Real oil price','World oil production','World oil inventories','World industrial production','U.S. industrial production','U.S. CPI','Fed funds rate'};
varNames = varNames_paper;

% select sample
smplStartInd = find(strcmp(sampleDates,smplStart));
smplEndInd   = find(strcmp(sampleDates,smplEnd));

data = data(smplStartInd:smplEndInd,:);
dataExo = dataExo(smplStartInd:smplEndInd,:);
sampleDates = sampleDates(smplStartInd:smplEndInd,:);
sampleDatesNum = sampleDatesNum(smplStartInd:smplEndInd,:);


%% External instruments VAR

% load the proxy
loadProxy;

% run proxy VAR
runProxyVAR; 

%% Figure
time = (0:horizon)';    

figure('Position',[100 100 1050 500],'PaperPositionMode','Auto','DefaultAxesFontSize',13);
signIRFs = 1;
for j=1:nvar % variable
    h(j) = subplot(ceil(nvar/4),4,j);
    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper_proxy(1,j); IRFslower_proxy(1:end,j); flipud([IRFsupper_proxy(1:end,j); IRFslower_proxy(end,j)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.2);
    set(hh,'edgecolor','none');

    hold on;
    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper2_proxy(1,j); IRFslower2_proxy(1:end,j); flipud([IRFsupper2_proxy(1:end,j); IRFslower2_proxy(end,j)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.4);
    set(hh,'edgecolor','none');
    if includeBase && j<7
        load('../../results/IRFsbench')
        p2=plot(time, signIRFs*IRFs_base(:,j), 'Linewidth', 1.5,'Color',[0.8500, 0.3250, 0.0980],'LineStyle',':');
    end
    p1=plot(time, signIRFs*IRFs_proxy(:,j),'k', 'Linewidth', 1.5); hold on;
    if ~ismember(0,get(gca,'ylim'))
        line(get(gca,'xlim'),[0 0],'Color','k')
    end
    grid on ;hold off;
    title(varNames_paper{j}) 
    if dataFrequency == 'M'
        xlabel('Months');
    elseif dataFrequency == 'Q'
        xlabel('Quarters');
    end   
    ylabel('\%');
    xlim([0,horizon]);
    xticks([0:10:horizon]);
    if j==1
        legend([p1 p2],{'FFR Model','Baseline'})
    end
end
pause(0.001)
pos = get(h,'Position');
new = mean(cellfun(@(v)v(1),pos(1:2)));
set(h(j-2),'Position',[(pos{1}(1)+pos{2}(1))/2 pos{end}(2:end)])
set(h(j-1),'Position',[(pos{2}(1)+pos{3}(1))/2 pos{end}(2:end)])
set(h(j),'Position',[(pos{3}(1)+pos{4}(1))/2 pos{end}(2:end)])

h=axes('Position',[0.25,0,.675,.5],'Xlim',[0 1],'Ylim',[0 1]);
set(h,'Visible','off');
if k==1
    string_1stage = ['First stage regression: F: ',num2str(olsEst.F,' %2.2f'),', robust F: ',num2str(olsEst.Frobust,' %2.2f'),', $R^2$: ',num2str(olsEst.R2*100,' %1.2f'),'\%, Adjusted $R^2$: ',num2str(olsEst.R2adj*100,' %1.2f'),'\%'];
    text('Position',[-.07 -0.05],'string',string_1stage,'FontSize',14);
end
tightfig;
if saveFigs
    print('-dpdf', gcf, '../../results/appendix/figurea9');  
end
