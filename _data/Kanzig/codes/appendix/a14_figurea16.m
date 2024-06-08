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
nsim       = 1000;      % number of simulations in bootstrap
bootType   = 'mbb1block'; 

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

contracts = [2 4 7 9 13 14];
maxcontract = length(contracts);
IRFs_extended = nan(horizon+1,nvar,maxcontract);
IRFsupper_extended = nan(horizon+1,nvar,maxcontract);
IRFslower_extended = nan(horizon+1,nvar,maxcontract);

isave = 1;
for ncontract = contracts
    
    % load the proxy
    loadProxy;
    
    runProxyVAR; 

    % save outputs
    IRFs_extended(:,:,isave) = IRFs_proxy;   % unit shock IRFs to pos oil price shock
    IRFsupper_extended(:,:,isave) = IRFsupper_proxy;
    IRFslower_extended(:,:,isave) = IRFslower_proxy;

    isave = isave+1;
end

%% Figure
time = (0:horizon)';    % time horizon for IRFs and FEVDs

colorz = linspecer(20);
colorz(6,:) = [0 0 0];

figure('Position',[100 100 1000 600],'PaperPositionMode','Auto','DefaultAxesFontSize',13);
signIRFs = 1;
for j=1:nvar % variable
    h(j) = subplot(ceil(nvar/3),3,j);
    for kk = [1 2 3 4 5]
    p(kk) = plot(time, signIRFs*IRFs_extended(:,j,kk),'Color',colorz(kk,:), 'Linewidth', 1.5,'LineStyle',':'); hold on;
    end
    kk=6;
    p(kk) = plot(time, signIRFs*IRFs_extended(:,j,kk),'Color',colorz(kk,:), 'Linewidth', 2)
    if ~ismember(0,get(gca,'ylim'))
        line(get(gca,'xlim'),[0 0],'Color','k')
    end
    if j==1
        legend([p(6) p(1) p(2) p(3) p(4) p(5) ],{'COMP','1M','3M','6M','9M','12M'})
    end
    grid on ;hold off;
    title(varNames_paper{j}) 
    xlabel('Months');
    xticks([0:10:horizon])
    ylabel('\%');
    xlim([0,horizon]);
end
pause(0.001)
if mod(nvar,2)~=0
    pos = get(h,'Position');
    new = mean(cellfun(@(v)v(1),pos(1:2)));
    set(h(j-1),'Position',[(pos{1}(1)+pos{2}(1))/2 pos{end}(2:end)])
    set(h(j),'Position',[(pos{2}(1)+pos{3}(1))/2 pos{end}(2:end)])
end
h=axes('Position',[0,0,1,1],'Xlim',[0 1],'Ylim',[0 1]);
set(h,'Visible','off');
if saveFigs
    saveas(gcf,'../../results/appendix/figurea16','epsc2')
end

