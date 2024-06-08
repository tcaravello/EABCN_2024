%% Replication files for "The macroeconomic effects of oil supply news"
% This file generates figure A.11

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

% preallocate
IRFs_collect = nan(horizon+1,nvar,2);
IRFsupper_collect  = nan(horizon+1,nvar,2);
IRFslower_collect  = nan(horizon+1,nvar,2);
IRFsupper2_collect = nan(horizon+1,nvar,2);
IRFslower2_collect = nan(horizon+1,nvar,2);
   
% run with raw and refined instrument
ii = 1;
for instRefine = [false true]
    % load the proxy
    loadProxyCleaned;

    runProxyVAR; 

    % save outputs
    IRFs_collect(:,:,ii) = IRFs_proxy;  
    IRFsupper_collect(:,:,ii) = IRFsupper_proxy;
    IRFslower_collect(:,:,ii) = IRFslower_proxy;
    IRFsupper2_collect(:,:,ii) = IRFsupper2_proxy;
    IRFslower2_collect(:,:,ii) = IRFslower2_proxy;
    ii= ii+1;
end


%% Figure 
time = (0:horizon)';    % time horizon for IRFs and FEVDs

figure('Position',[100 100 1000 600],'PaperPositionMode','Auto','DefaultAxesFontSize',13);
signIRFs = 1;
for j=1:nvar % variable
    h(j) = subplot(ceil(nvar/3),3,j); 
    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper_collect(1,j,2); IRFslower_collect(1:end,j,2); flipud([IRFsupper_collect(1:end,j,2); IRFslower_collect(end,j,2)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.2);
    set(hh,'edgecolor','none'); 

    hold on;
    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper2_collect(1,j,2); IRFslower2_collect(1:end,j,2); flipud([IRFsupper2_collect(1:end,j,2); IRFslower2_collect(end,j,2)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.4);
    set(hh,'edgecolor','none');
    
    plot(time, signIRFs*IRFsupper_collect(:,j,1),'Color',[0.8500, 0.3250, 0.0980], 'Linestyle','--');
    plot(time, signIRFs*IRFslower_collect(:,j,1),'Color',[0.8500, 0.3250, 0.0980], 'Linestyle','--');
    plot(time, signIRFs*IRFsupper2_collect(:,j,1),'Color',[0.8500, 0.3250, 0.0980], 'Linestyle',':');
    plot(time, signIRFs*IRFslower2_collect(:,j,1),'Color',[0.8500, 0.3250, 0.0980], 'Linestyle',':');
        
    p1=plot(time, signIRFs*IRFs_collect(:,j,1), 'Linewidth', 1.5,'Color',[0.8500, 0.3250, 0.0980]); hold on;
    p2=plot(time, signIRFs*IRFs_collect(:,j,2),'k', 'Linewidth', 1.5); hold on;
    
    if ~ismember(0,get(gca,'ylim'))
        line(get(gca,'xlim'),[0 0],'Color','k')
    end
    grid on ;hold off;
    title(varNames_paper{j}) 
    xlabel('Months');
    ylabel('\%');
    xlim([0,horizon]);
    xticks([0:10:horizon]);
    if j==1
        legend([p1 p2],{'Raw','Refined'})
    end
    if j==4
        ylim([-1.5 1.5])
    end
    if j==6
        ylim([-0.5 1.2])
    end
end
tightfig;
if saveFigs
    print('-dpdf', gcf, '../../results/appendix/figurea11');  
end

