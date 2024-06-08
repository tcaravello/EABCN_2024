%% Replication files for "The macroeconomic effects of oil supply news"
% This figure generates Figure A.4 Panel B

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

%% Run LP-IV Robustness

% load the proxy
loadProxy;

% settings
horizonLP = horizon;

colval = [0.8500, 0.3250, 0.0980]; 

% lag specifications
pSet = [12 18 24];

% all lags as controls
IRFs_LPIVs = zeros(horizonLP+1,nvar,length(pSet));

kk = 1;
for pLP = pSet
    for hh = 0:horizonLP
        for ii = 1:nvar

            % first clean for lags on longer sample
            yi_pre = data(pLP+1+hh:end-horizonLP+hh,ii);
            xi_pre = data(pLP+1:end-horizonLP,1);

            Xr_pre = [];
            for jj = 1:pLP
                Xr_pre = [Xr_pre data(pLP+1-jj:end-horizonLP-jj,:)];
            end
            Xr_pre = [Xr_pre dataExo(pLP+1:end-horizonLP,1:end)];

            % construct shorter sample for LP-IV on shock
            yi_all = [nan(pLP,1); yi_pre - Xr_pre*(Xr_pre\yi_pre); nan(horizonLP,1)];
            xi_all = [nan(pLP,1); xi_pre - Xr_pre*(Xr_pre\xi_pre); nan(horizonLP,1)];
            yi = yi_all(smplStartProxyInd:smplEndProxyInd-horizonLP);
            Xr = xi_all(smplStartProxyInd:smplEndProxyInd-horizonLP);
            Zr = proxyRaw(smplStartProxyInd:smplEndProxyInd-horizonLP);                     
            
            gmmLP = gmmest(Xr,Zr,yi,true,2,hh+1);
            IRFs_LPIVs(hh+1,ii,kk) = gmmLP.bhat(1);
         end
    end
    
    % standardize
    estSize = IRFs_LPIVs(1,1,kk);
    IRFs_LPIVs(:,:,kk) = IRFs_LPIVs(:,:,kk)./estSize*shockSize;
    
    kk = kk+1;
end

% flexible lags as controls
IRFs_LPsflex = zeros(horizonLP+1,nvar,length(pSet));

kk = 1;
for pLP = pSet
    for hh = 0:horizonLP
        for ii = 1:nvar
            % first clean for lags on longer sample
            yi_pre = data(pLP+1+hh:end-horizonLP+hh,ii);
            xi_pre = data(pLP+1:end-horizonLP,1);

            Xr_pre = [];
            for jj = 1:pLP
                Xr_pre = [Xr_pre data(pLP+1-jj:end-horizonLP-jj,[1:4])];
                if ii>4
                    Xr_pre = [Xr_pre data(pLP+1-jj:end-horizonLP-jj,ii)];
                end
            end
            Xr_pre = [Xr_pre dataExo(pLP+1:end-horizonLP,1:end)];

            % construct shorter sample for LP-IV on shock
            yi_all = [nan(pLP,1); yi_pre - Xr_pre*(Xr_pre\yi_pre); nan(horizonLP,1)];
            xi_all = [nan(pLP,1); xi_pre - Xr_pre*(Xr_pre\xi_pre); nan(horizonLP,1)];
            yi = yi_all(smplStartProxyInd:smplEndProxyInd-horizonLP);
            Xr = xi_all(smplStartProxyInd:smplEndProxyInd-horizonLP);
            Zr = proxyRaw(smplStartProxyInd:smplEndProxyInd-horizonLP);                      

            gmmLP = gmmest(Xr,Zr,yi,true,2,hh+1);
            IRFs_LPsflex(hh+1,ii,kk) = gmmLP.bhat(1);
         end
    end
    
    % standardize
    estSize = IRFs_LPsflex(1,1,kk);
    IRFs_LPsflex(:,:,kk) = IRFs_LPsflex(:,:,kk)./estSize*shockSize;
    
    kk = kk+1;
end

%% Figure
colorz = linspecer(10);

time = (0:horizonLP)';

figure('Position',[100 100 1000 600],'PaperPositionMode','Auto','DefaultAxesFontSize',13)
for j=1:nvar %variable
    subplot(2,ceil(nvar/2),j) 


    hold on
    for ii = 1:size(IRFs_LPIVs,3)
        plot(time, IRFs_LPIVs(:,j,ii),'k', 'Linewidth', 1.5, 'LineStyle','--','Color',colorz(ii,:)); 
    end
    for ii = 1:size(IRFs_LPsflex,3)
        plot(time, IRFs_LPsflex(:,j,ii),'k', 'Linewidth', 1.5,'Color',colorz(ii,:)); 
    end
    grid on;box on;hold off;
    if j==1
        ylim([-20 40])
        leg = legend({'12 lags','18 lags','24 lags','12 lags, flex.','18 lags, flex.','24 lags, flex.'},'autoupdate','off','fontsize',8);
        leg.ItemTokenSize = [18,9];
    elseif j==2 || j==4 || j==5
        ylim([-2 2])
    elseif j==6
        ylim([-0.4 1])
    end
    title(varNames_paper{j}) 
    xlabel('Months');
    ylabel('\%');
    if ~ismember(0,get(gca,'ylim'))
        line(get(gca,'xlim'),[0 0],'Color','k')
    end
    xlim([0,horizonLP]);
    xticks([0:10:horizonLP]);
end
tightfig;
if saveFigs
    print('-dpdf', gcf, '../../results/appendix/figurea4b');  
end

