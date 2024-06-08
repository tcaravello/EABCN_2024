%% Replication files for "The macroeconomic effects of oil supply news"
% Step 3: This file creates figures 3 and 5 in the paper

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
smplStartProxy = '1975M01'; 
smplEndProxy   = '2017M12'; 
% Missing values in proxy are censored to zero

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


%% Read in data
load('../data/dataBaseM')
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

figure('Position',[10 10 1000 600],'PaperPositionMode','Auto','DefaultAxesFontSize',13); 
for ii = 1:size(data,2)
    subplot(2,ceil(nvar/2),ii);
    hold on
    plot(sampleDatesNum,data(:,ii),'LineWidth',1.5)
    pylims = get(gca,'ylim');
    if ~ismember(0,pylims) && ~(pylims(1)>0 || pylims(2)<0)
        l1 = line(get(gca,'xlim'),[0 0],'Color','k');
        uistack(l1,'bottom');
    end
    title(varNames_paper{ii})
    xlim([sampleDatesNum(1) sampleDatesNum(end)])
    grid on
    box on
end
tightfig;
if saveFigs
    print('-dpdf', gcf, '../results/appendix/figureb1');  
end
    
%% External instruments VAR

% load the proxy
loadProxy;

% run proxy VAR
runProxyVAR; 

% save results to disk
if saveFigs
    IRFs_base = IRFs_proxy; 
    IRFsupper_base = IRFsupper_proxy; 
    IRFslower_base = IRFslower_proxy; 
    IRFsupper2_base = IRFsupper2_proxy; 
    IRFslower2_base = IRFslower2_proxy; 

    save('../results/IRFsbench','IRFs_base','IRFsupper_base','IRFslower_base','IRFsupper2_base','IRFslower2_base')
end

%% Figure 3: IRFs
time = (0:horizon)';    % time horizon for IRFs and FEVDs

figure('Position',[100 100 1000 600],'PaperPositionMode','Auto','DefaultAxesFontSize',13);
signIRFs = 1;
for j=1:nvar %variable
    h(j) = subplot(ceil(nvar/3),3,j); 

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
    title(varNames_paper{j}) 
    ylabel('\%');
    xlim([0,horizon]);
    if dataFrequency == 'M'
        xlabel('Months');
        xticks(0:10:horizon)
    elseif dataFrequency == 'Q'
        xlabel('Quarters');
        xticks(0:5:horizon)
    end   
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
    print('-dpdf', gcf, '../results/figure3');  
end


%% Historical decomposition
% get shock as in Stock and Watson (2018)
if strcmp(shockType,'custom')
    % unit normalization
    oilSupplyNewsShock = (b1unit'*inv2(Sigma)*U')'*inv2(b1unit'*inv2(Sigma)*b1unit);
else
    % one sd normalization
    oilSupplyNewsShock = (b1'*inv2(Sigma)*U')'*1; 
end

% Do historical decomposition of oil price
varCompan = varxcompan(varEst);    % companion form
F = varCompan.F;                   % max(abs(eig(F))): VAR is stable at least.
V       = oilSupplyNewsShock';     % structural errors 
nvarXeq = nvar*p;                  % number of lagged endogenous per equation

% Compute historical decompositions

% Contribution of the shock
A0_big = zeros(nvarXeq,nvar);
if strcmp(shockType,'custom')
    A0_big(1:nvar,1) = b1unit;
else
    A0_big(1:nvar,1) = b1;
end
Icomp = [eye(nvar) zeros(nvar,(p-1)*nvar)];
HDshock_big = zeros(p*nvar,T+1);
HDshock_temp = zeros(nvar,T+1);
V_big = zeros(nvar,T+1); % matrix of shocks conformable with companion
V_big(1,2:end) = V;
for i = 2:T+1
    HDshock_big(:,i) = A0_big*V_big(:,i) + F*HDshock_big(:,i-1);
    HDshock_temp(:,i) =  Icomp*HDshock_big(:,i);
end
HDshock = HDshock_temp(:,2:end)';

% bootstrap it 
bootHDshock = zeros(T,nvar,nsim);
for isim = 1:nsim
    varBoot = struct;
    varBoot.B = bootBs(:,:,isim);
    varBoot.Sigma = Sigma;
    varBoot.nexo = nexo;
    varBootCompan = varxcompan(varEst);    % companion form
    bootF = varBootCompan.F;
    bootV = bootShocks(:,isim)';                         % structural errors 

    % Contribution of the shock
    bootA0_big = zeros(nvarXeq,nvar);
    bootA0_big(1:nvar,1) = bootb1s(:,isim);
    Icomp = [eye(nvar) zeros(nvar,(p-1)*nvar)];
    bootHDshock_big = zeros(p*nvar,T+1);
    bootHDshock_temp = zeros(nvar,T+1);
    bootV_big = zeros(nvar,T+1); % matrix of shocks conformable with companion
    bootV_big(1,2:end) = bootV;
    for i = 2:T+1
        bootHDshock_big(:,i) = bootA0_big*bootV_big(:,i) + bootF*bootHDshock_big(:,i-1);
        bootHDshock_temp(:,i) =  Icomp*bootHDshock_big(:,i);
    end
    bootHDshock(:,:,isim) = bootHDshock_temp(:,2:end)';
end

% get quantiles and recenter around point estimate
HDshockmed = quantile(bootHDshock, 0.5, 3);
HDshockupper = quantile(bootHDshock, 1-alpha/2, 3)-HDshockmed+HDshock;
HDshocklower = quantile(bootHDshock, alpha/2, 3)-HDshockmed+HDshock;
HDshockupper2 = quantile(bootHDshock, 1-alpha2/2, 3)-HDshockmed+HDshock;
HDshocklower2 = quantile(bootHDshock, alpha2/2, 3)-HDshockmed+HDshock;


%% Figure 5: Historical decomposition
timeHD = sampleDatesNum(smplStartProxyVARInd:smplEndProxyVARInd);
oilDates = [1978+8/12 1980+9/12 1985+11/12 1990+7/12 1997+6/12 2002+10/12 2008+8/12 2014+10/12];
figure('Position',[100 100 1000 350],'DefaultAxesFontSize',13)
hold on
hh=fill([timeHD(1); timeHD(1:end); flipud([timeHD(1:end); timeHD(end)])],[HDshockupper(1,1); HDshocklower(1:end,1); flipud([HDshockupper(1:end,1); HDshocklower(end,1)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.2);
set(hh,'edgecolor','none'); 

hh=fill([timeHD(1); timeHD(1:end); flipud([timeHD(1:end); timeHD(end)])],[HDshockupper2(1,1); HDshocklower2(1:end,1); flipud([HDshockupper2(1:end,1); HDshocklower2(end,1)])],[0.1, 0.4470, 0.7410]); 
set(hh,'facealpha',.4);
set(hh,'edgecolor','none'); 
set(hh,'edgealpha',.4);
f1=plot(sampleDatesNum(smplStartProxyVARInd:smplEndProxyVARInd),varEst.Y(:,1)-mean(varEst.Y(:,1)),'Color',[0.8500, 0.3250, 0.0980],'LineWidth',2,'LineStyle',':');
f2=plot(sampleDatesNum(smplStartProxyVARInd:smplEndProxyVARInd),HDshock(:,1)-mean(HDshock(:,1)),'Color','k','LineWidth',1.5);
for ii=1:length(oilDates)
   xline(oilDates(ii),'LineWidth',1.2); 
end
grid on
box on
xlim([1975 2017])
ylim([-150 150])
ylabel('\%')
tightfig;

fs = 10;
legend([f1 f2],'Real oil price','Contribution oil supply news','Location','Southeast','Interpreter','latex')
str={'Iranian', 'revolution'};
annotation('textbox',[0.077, 0.85, 0.1, 0.1],'String',str,'FitBoxToText','on','EdgeColor','none','Interpreter','latex','fontsize',fs);
str={'Iran-Iraq','war'};
annotation('textbox',[0.2, 0.85, 0.1, 0.1],'String',str,'FitBoxToText','on','EdgeColor','none','Interpreter','latex','fontsize',fs);
str={'OPEC','collapse'};
annotation('textbox',[0.31, 0.85, 0.1, 0.1],'String',str,'FitBoxToText','on','EdgeColor','none','Interpreter','latex','fontsize',fs);
str={'Gulf war'};
annotation('textbox',[0.414, 0.85, 0.1, 0.1],'String',str,'FitBoxToText','on','EdgeColor','none','Interpreter','latex','fontsize',fs);
str={'Asian fin.','crisis'};
annotation('textbox',[0.49, 0.85, 0.1, 0.1],'String',str,'FitBoxToText','on','EdgeColor','none','Interpreter','latex','fontsize',fs);
str={'Venezuelan','crisis'};
annotation('textbox',[0.596, 0.85, 0.1, 0.1],'String',str,'FitBoxToText','on','EdgeColor','none','Interpreter','latex','fontsize',fs);
str={'Global fin.','crisis'};
annotation('textbox',[0.73, 0.85, 0.1, 0.1],'String',str,'FitBoxToText','on','EdgeColor','none','Interpreter','latex','fontsize',fs);
str={'Oil crash','2014'};
annotation('textbox',[0.872, 0.85, 0.1, 0.1],'String',str,'FitBoxToText','on','EdgeColor','none','Interpreter','latex','fontsize',fs);

if saveFigs
    print('-dpdf', gcf, '../results/figure5');
end
