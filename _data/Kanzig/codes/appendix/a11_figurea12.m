%% Replication files for "The macroeconomic effects of oil supply news"
% This file generates figure A.12

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

% number of variables in VAR
nvar = size(data,2);  

% names for paper
varNames_paper = {'Real oil price','World oil production','World oil inventories','World industrial production','U.S. industrial production','U.S. CPI'};
varNames_paperVD = {'Real oil price','Oil production','Oil inventories','World IP','U.S. IP','U.S. CPI'};

% select sample
smplStartInd = find(strcmp(sampleDates,smplStart));
smplEndInd   = find(strcmp(sampleDates,smplEnd));

data = data(smplStartInd:smplEndInd,:);
dataExo = dataExo(smplStartInd:smplEndInd,:);
sampleDates = sampleDates(smplStartInd:smplEndInd,:);
sampleDatesNum = sampleDatesNum(smplStartInd:smplEndInd,:);

%% External instruments VAR

% load the proxy
loadProxyOrdinary;

runProxyVAR; 

%% Figure
time = (0:horizon)';    % time horizon for IRFs and FEVDs

figure('Position',[100 100 1000 600],'PaperPositionMode','Auto','DefaultAxesFontSize',13);
signIRFs = 1;
for j=1:nvar % variable
    h(j) = subplot(ceil(nvar/3),3,j); 

    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper_proxy(1,j); IRFslower_proxy(1:end,j); flipud([IRFsupper_proxy(1:end,j); IRFslower_proxy(end,j)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.2);
    set(hh,'edgecolor','none'); 

    hold on;
    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper2_proxy(1,j); IRFslower2_proxy(1:end,j); flipud([IRFsupper2_proxy(1:end,j); IRFslower2_proxy(end,j)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.4);
    set(hh,'edgecolor','none');

    if includeBase
        load('../../results/IRFsbench')
        p2=plot(time, signIRFs*IRFs_base(:,j), 'Linewidth', 2,'Color',[0.8500, 0.3250, 0.0980],'LineStyle',':');
    end

    p1=plot(time, signIRFs*IRFs_proxy(:,j),'k', 'Linewidth', 1.5); 
    if ~ismember(0,get(gca,'ylim'))
        line(get(gca,'xlim'),[0 0],'Color','k')
    end
    box on
    grid on ;hold off;
    title(varNames_paper{j}) 
    ylabel('\%');
    xlim([0,horizon]);
    xlabel('Months');
    xticks(0:10:horizon)
    if j==1 && includeBase
        legend([p1 p2],{'Ordinary ann.','Baseline'})
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
    print('-dpdf', gcf, '../../results/appendix/figurea12');  
end

