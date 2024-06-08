%% Replication files for "The macroeconomic effects of oil supply news"
% This figure generates Figure A.5

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
smplStart = '1974M01'; 
smplEnd   = '2017M12'; 

% Instrument sample
% range has to be contained in estimation sample
smplStartProxy = '1983M04'; 
smplEndProxy   = '2017M12'; 

% VAR specifics
p          = 18;       
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
includeBase = true;

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


%% Load instrument and control series

% load the proxy
loadProxy;

% load placebo
load('../../instrument/OilSurprisesMLogControl')
proxyControlRaw = [oilPlaceboWTIM(:,ncontract)]; 

%% Run heteroskedasticity-based LP

% settings
pLP = p;
horizonLP = horizon;

variableControls = true;
colval = [0.8500, 0.3250, 0.0980]; 

IRFs_LP = zeros(horizonLP+1,nvar);
IRFsupper_LP = zeros(horizonLP+1,nvar);
IRFslower_LP = zeros(horizonLP+1,nvar);
IRFsupper2_LP = zeros(horizonLP+1,nvar);
IRFslower2_LP = zeros(horizonLP+1,nvar);
for hh = 0:horizonLP
    for ii = 1:nvar

        % first clean for lags on longer sample
        yi_pre = data(pLP+1+hh:end-horizonLP+hh,ii);

        Xr_pre = [];
        for jj = 1:pLP
            if variableControls
                Xr_pre = [Xr_pre data(pLP+1-jj:end-horizonLP-jj,[1:4])];
                if ii>4
                    Xr_pre = [Xr_pre data(pLP+1-jj:end-horizonLP-jj,ii)];
                end
            else
                Xr_pre = [Xr_pre data(pLP+1-jj:end-horizonLP-jj,:)];
            end
        end
        Xr_pre = [Xr_pre dataExo(pLP+1:end-horizonLP,1:end)];

        % construct shorter sample for LP on shock
        yi_all = [nan(pLP,1); yi_pre - Xr_pre*(Xr_pre\yi_pre); nan(horizonLP,1)];
        yi = yi_all(smplStartProxyVARInd:smplEndProxyVARInd-horizonLP);
        Xr = proxyRaw(smplStartProxyInd:smplEndProxyInd-horizonLP);                      
        XrC = proxyControlRaw(smplStartProxyInd:smplEndProxyInd-horizonLP); 

        % regimes
        indsR1 = logical(statementMind(smplStartProxyInd:smplEndProxyInd-horizonLP));
        indsR2 = logical(placeboMind(smplStartProxyInd:smplEndProxyInd-horizonLP));
  
        T_OPEC = size(Xr(indsR1,:),1);
        T_Control = size(XrC(indsR2,:),1);

        % IV estimator
        XrIV = [(Xr(indsR1,:)-mean(Xr(indsR1,:)))/sqrt(T_OPEC); (XrC(indsR2,:)-mean(XrC(indsR2,:)))/sqrt(T_Control)];
        ZrIV = [(Xr(indsR1,:)-mean(Xr(indsR1,:)))/sqrt(T_OPEC); -(XrC(indsR2,:)-mean(XrC(indsR2,:)))/sqrt(T_Control)];
        yiIV = [(yi(indsR1,:)-mean(yi(indsR1,:)))/sqrt(T_OPEC); (yi(indsR2,:)-mean(yi(indsR2,:)))/sqrt(T_Control)];

        gmmLP = gmmest(XrIV,ZrIV,yiIV,false,2,hh+1);
        IRFs_LP(hh+1,ii) = gmmLP.bhat(1);
        IRFsupper_LP(hh+1,ii) = gmmLP.bhat(1) + norminv(1-alpha/2,0,1)*(gmmLP.varbhat(1,1))^0.5;
        IRFslower_LP(hh+1,ii) = gmmLP.bhat(1) - norminv(1-alpha/2,0,1)*(gmmLP.varbhat(1,1))^0.5;
        IRFsupper2_LP(hh+1,ii) = gmmLP.bhat(1) + norminv(1-alpha2/2,0,1)*(gmmLP.varbhat(1,1))^0.5;
        IRFslower2_LP(hh+1,ii) = gmmLP.bhat(1) - norminv(1-alpha2/2,0,1)*(gmmLP.varbhat(1,1))^0.5;

    end
end

% standardize
estSize = IRFs_LP(1,1);
IRFs_LP = IRFs_LP./estSize*shockSize;
IRFsupper_LP = IRFsupper_LP./estSize*shockSize;
IRFslower_LP = IRFslower_LP./estSize*shockSize;
IRFsupper2_LP = IRFsupper2_LP./estSize*shockSize;
IRFslower2_LP = IRFslower2_LP./estSize*shockSize;


%% Figure
time = (0:horizonLP)';

figure('Position',[100 100 1000 600],'PaperPositionMode','Auto','DefaultAxesFontSize',13)
signIRFs = 1;
for j=1:nvar % variable
    subplot(2,ceil(nvar/2),j) 

    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper_LP(1,j); IRFslower_LP(1:end,j); flipud([IRFsupper_LP(1:end,j); IRFslower_LP(end,j)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.2);
    set(hh,'edgecolor','none');
     hold on;

    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper2_LP(1,j); IRFslower2_LP(1:end,j); flipud([IRFsupper2_LP(1:end,j); IRFslower2_LP(end,j)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.4);
    set(hh,'edgecolor','none');

    if includeBase
        load('../../results/appendix/IRFsbenchLP')
        p2=plot(time, signIRFs*IRFs_LPbase(1:horizonLP+1,j), 'Linewidth', 1.5,'Color',[0.8500, 0.3250, 0.0980]);
        plot(time, signIRFs*IRFsupper_LPbase(1:horizonLP+1,j), 'Color',[0.8500, 0.3250, 0.0980],'LineStyle','--');
        plot(time, signIRFs*IRFslower_LPbase(1:horizonLP+1,j),'Color',[0.8500, 0.3250, 0.0980],'LineStyle','--');
        plot(time, signIRFs*IRFsupper2_LPbase(1:horizonLP+1,j),'Color',[0.8500, 0.3250, 0.0980],'LineStyle',':');
        plot(time, signIRFs*IRFslower2_LPbase(1:horizonLP+1,j),'Color',[0.8500, 0.3250, 0.0980],'LineStyle',':');
    end

    p1=plot(time, signIRFs*IRFs_LP(:,j),'k', 'Linewidth', 1.5); 
    grid on ;hold off;
    if j==1
        ylim([-20 60])
        legend([p1 p2],'Rigobon-LP','LP-IV','AutoUpdate','off')
    elseif j==2 
        ylim([-3 3])
    elseif j==4
        ylim([-3 3])
    elseif j==5
        ylim([-4 2.5])
    elseif j==6
        ylim([-1 2.5])
    end
    title(varNames_paper{j}) 
    if dataFrequency == 'M'
        xlabel('Months');
    elseif dataFrequency == 'Q'
        xlabel('Quarters');
    end  
    ylabel('\%');
    if ~ismember(0,get(gca,'ylim'))
        line(get(gca,'xlim'),[0 0],'Color','k')
    end
    xlim([0,horizonLP]);
    xticks([0:10:horizonLP]);
end
tightfig;
if saveFigs
    print('-dpdf', gcf, '../../results/appendix/figurea5');  
end

