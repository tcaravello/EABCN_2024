%% Replication files for "The macroeconomic effects of oil supply news"
% Step 5: This file creates figure 4b in the paper

% Diego R. Känzig
% LBS, September 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
%close all
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
smplStartProxy = '1975M01'; 
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
includeBase = false;     % Include baseline response in plots

verbo = false;           % Control verbosity
saveFigs   = true;       % Save figures to disk
doLagCrit = false;

%% Read in data
load('empirics/_data/Kanzig/data/dataBaseM')
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
loadProxy;

runProxyVAR;  

% get shock
oilSupplyNewsShock = (b1'*inv2(Sigma)*U')'; 


%% Compute IRFs using LPs (Jorda, 2005)

pLP = 1;
horizonLP = horizon;

colval = [0.8500, 0.3250, 0.0980]; 

% reduced-form
IRFs_LP = zeros(horizonLP+1,nvar);

for hh = 0:horizonLP
    for ii = 1:nvar
        
        yi = data(smplStartProxyVARInd+hh:end-horizonLP+hh,ii);

        Xr = oilSupplyNewsShock(1:end-horizonLP);  
        for jj = 1:pLP
            Xr = [Xr data(smplStartProxyVARInd-jj:end-horizonLP-jj,ii)];
        end
        Xr = [Xr dataExo(smplStartProxyInd:smplEndProxyInd-horizonLP,2:end)];

        olsLP = olsest(Xr,yi,true);
        IRFs_LP(hh+1,ii) = olsLP.bhat(1);
    end
end
IRFs_LP = IRFs_LP./IRFs_LP(1,1)*shockSize;

% bands
IRFsboot_LP = zeros(horizonLP+1,nvar,nsim);

for j = 1:nsim
    for hh = 0:horizonLP
        for ii = 1:nvar
            
            yi = bootDatas(smplStartProxyVARInd+hh:end-horizonLP+hh,ii,j);

            Xr = bootShocks(1:end-horizonLP,j);  
            for jj = 1:pLP
                Xr = [Xr bootDatas(smplStartProxyVARInd-jj:end-horizonLP-jj,ii,j)];
            end
            Xr = [Xr dataExo(smplStartProxyInd:smplEndProxyInd-horizonLP,2:end)];

            olsLP = olsest(Xr,yi);
            IRFsboot_LP(hh+1,ii,j) = olsLP.bhat(1);
        end
    end
    IRFsboot_LP(:,:,j) = IRFsboot_LP(:,:,j)./IRFsboot_LP(1,1,j)*IRFs_LP(1,1);
end

% center bootstrapped IRFs around sample estimates and get quantiles
IRFsmed_LP = quantile(IRFsboot_LP, 0.5, 3);

IRFsupper_LP = quantile(IRFsboot_LP, 1-alpha/2, 3)-IRFsmed_LP+IRFs_LP;  % rescaling does not affect the ordering, thus it is fine to take quantile first
IRFslower_LP = quantile(IRFsboot_LP, alpha/2, 3)-IRFsmed_LP+IRFs_LP;

IRFsupper2_LP = quantile(IRFsboot_LP, 1-alpha2/2, 3)-IRFsmed_LP+IRFs_LP;  
IRFslower2_LP = quantile(IRFsboot_LP, alpha2/2, 3)-IRFsmed_LP+IRFs_LP;

time = (0:horizonLP)';

%% Figure 4b
figure('Position',[10 10 1000 600],'PaperPositionMode','Auto','DefaultAxesFontSize',13);
if IRFs_LP(1,1)>0
    signIRFs = 1;
else
    signIRFs = -1;
end
for j=1:nvar % variable
    subplot(2,ceil(nvar/2),j) 

    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper_LP(1,j); IRFslower_LP(1:end,j); flipud([IRFsupper_LP(1:end,j); IRFslower_LP(end,j)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.2);
    set(hh,'edgecolor','none'); 

    hold on;
    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper2_LP(1,j); IRFslower2_LP(1:end,j); flipud([IRFsupper2_LP(1:end,j); IRFslower2_LP(end,j)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.4);
    set(hh,'edgecolor','none');

    p2=plot(time, signIRFs*IRFs_proxy(1:horizonLP+1,j), 'Color',colval, 'Linewidth', 1.5,'LineStyle','-');
    p1=plot(time, signIRFs*IRFs_LP(:,j),'k', 'Linewidth', 1.5); 
    plot(time, signIRFs*IRFsupper_proxy(1:horizonLP+1,j), 'Color',colval,'Linestyle','--');
    plot(time, signIRFs*IRFslower_proxy(1:horizonLP+1,j), 'Color',colval,'Linestyle','--');
    plot(time, signIRFs*IRFsupper2_proxy(1:horizonLP+1,j), 'Color',colval,'Linestyle',':');
    plot(time, signIRFs*IRFslower2_proxy(1:horizonLP+1,j), 'Color',colval,'Linestyle',':');
    grid on ;hold off;
    if j==1
        legend('LP-IV','Proxy-VAR','AutoUpdate','off')
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
    if j==1
        legend([p1 p2],{'LP','VAR'})
    end
end
tightfig;
if saveFigs
    print('-dpdf', gcf, '../results/figure4b');  
end

