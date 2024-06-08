%% Replication files for "The macroeconomic effects of oil supply news"
% Step 3: This file creates figures 3 and 5 in the paper

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
smplStartProxy = '1974M01'; 
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
cumProxy = false;
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

%% Internal instruments VAR

% load the proxy
loadProxy;

runProxyVARPM; 

%% Figure

% plot the results
time = (0:horizon)';    % time horizon for IRFs and FEVDs

figure('Position',[100 100 1000 600],'PaperPositionMode','Auto','DefaultAxesFontSize',13);
if IRFs_proxy(1,1)>0
    signIRFs = 1;
else
    signIRFs = -1;
end
for j=1:nvar-1 % variable
    subplot(2,3,j); 
    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper_proxy(1,j+1); IRFslower_proxy(1:end,j+1); flipud([IRFsupper_proxy(1:end,j+1); IRFslower_proxy(end,j+1)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.2);
    set(hh,'edgecolor','none');

    hold on;
    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper2_proxy(1,j+1); IRFslower2_proxy(1:end,j+1); flipud([IRFsupper2_proxy(1:end,j+1); IRFslower2_proxy(end,j+1)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.4);
    set(hh,'edgecolor','none');

    if includeBase
        load('../../results/IRFsbench')
        p2=plot(time, signIRFs*IRFs_base(:,j), 'Linewidth', 1.5,'Color',[0.8500, 0.3250, 0.0980],'LineStyle','-');
        plot(time, signIRFs*IRFsupper_base(:,j), 'Color',[0.8500, 0.3250, 0.0980],'LineStyle','--');
        plot(time, signIRFs*IRFslower_base(:,j),'Color',[0.8500, 0.3250, 0.0980],'LineStyle','--');
        plot(time, signIRFs*IRFsupper2_base(:,j),'Color',[0.8500, 0.3250, 0.0980],'LineStyle',':');
        plot(time, signIRFs*IRFslower2_base(:,j),'Color',[0.8500, 0.3250, 0.0980],'LineStyle',':');
    end
    p1=plot(time, signIRFs*IRFs_proxy(:,j+1),'k', 'Linewidth', 1.5); hold on;

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
    if j==1
        ylim([-10 30])
    end
    if j==1
        ylim([-10 30])
    end
    if j==6
        ylim([-1 1])
    end
    xlim([0,horizon]);
    xticks([0:10:horizon]);
    if j==1 && includeBase
        legend([p1 p2],{'Internal inst.','External inst.'})
    end

end
pause(0.001)
h=axes('Position',[0.25,0,.5,.5],'Xlim',[0 1],'Ylim',[0 1]);
set(h,'Visible','off');
tightfig;
if saveFigs
    print('-dpdf', gcf, '../../results/appendix/figurea14');  
end

