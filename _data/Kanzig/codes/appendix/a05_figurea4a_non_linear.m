%% Replication files for "The macroeconomic effects of oil supply news"
% This figure generates Figure A.4 Panel A

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
p          = 12;       
horizon    = 50;
shockType  = 'custom';   % one standard deviation 'sd' or 'custom'
shockSize  = 1;         % if custom, specify shock size here
alpha      = 0.1;        % Significance level for bands (alpha=0.1 => 90% CIs (two SD))
alpha2     = 0.32;
nsim       = 10000;      % number of simulations in bootstrap
bootType   = 'mbb1block'; 

% proxy
ncontract = 14; 

% switches
includeBase = false;

saveFigs   = false;


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
varNames_paper = {'Real oil price','World oil production','World oil inventories','World industrial production','U.S. industrial production','U.S. CPI', 'Unemployment', 'Fed Funds'};
varNames_paperVD = {'Real oil price','Oil production','Oil inventories','World IP','U.S. IP','U.S. CPI'};

% select sample
smplStartInd = find(strcmp(sampleDates,smplStart));
smplEndInd   = find(strcmp(sampleDates,smplEnd));

data = data(smplStartInd:smplEndInd,:);
dataExo = dataExo(smplStartInd:smplEndInd,:);
sampleDates = sampleDates(smplStartInd:smplEndInd,:);
sampleDatesNum = sampleDatesNum(smplStartInd:smplEndInd,:);


%% Run LP-IV 

% load the proxy
loadProxy;

% settings
pLP = p;
horizonLP = horizon;

colval = [0.8500, 0.3250, 0.0980]; 

IRFs_LPIV = zeros(horizonLP+1,nvar);
IRFs_LPIV_non_lin = zeros(horizonLP+1,nvar);
IRFsupper_LPIV_non_lin = zeros(horizonLP+1,nvar);
IRFslower_LPIV_non_lin = zeros(horizonLP+1,nvar);

IRFsupper_LPIV = zeros(horizonLP+1,nvar);
IRFslower_LPIV = zeros(horizonLP+1,nvar);
IRFsupper2_LPIV = zeros(horizonLP+1,nvar);
IRFslower2_LPIV = zeros(horizonLP+1,nvar);
for hh = 0:horizonLP
    for ii = 1:nvar

        % first clean for lags on longer sample
        yi_pre = data(pLP+1+hh:end-horizonLP+hh,ii);
        xi_pre = data(pLP+1:end-horizonLP,1);

        Xr_pre = [];
        Xr_pre2 = [];
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
        z_vec_quantiles = quantile(Zr(Zr~=0),[0.2,0.8]);
        bound = max(abs(z_vec_quantiles));
        fun_size = @(x) ((x>bound).*(x-bound)+(x<-bound).*(x+bound));
        Zr_2 = [Zr,fun_size(Zr)];

        gmmLP = gmmest(Zr_2,Zr_2,yi,true,2,hh+1);
        IRFs_LPIV(hh+1,ii) = gmmLP.bhat(1);
        IRFsupper_LPIV(hh+1,ii) = gmmLP.bhat(1) + norminv(1-alpha/2,0,1)*(gmmLP.varbhat(1,1))^0.5;
        IRFslower_LPIV(hh+1,ii) = gmmLP.bhat(1) - norminv(1-alpha/2,0,1)*(gmmLP.varbhat(1,1))^0.5;
        IRFs_LPIV_non_lin(hh+1,ii) = gmmLP.bhat(2);
        IRFsupper_LPIV_non_lin(hh+1,ii) = gmmLP.bhat(2) + norminv(1-alpha/2,0,1)*(gmmLP.varbhat(2,2))^0.5;
        IRFslower_LPIV_non_lin(hh+1,ii) = gmmLP.bhat(2) - norminv(1-alpha/2,0,1)*(gmmLP.varbhat(2,2))^0.5;
        IRFsupper2_LPIV(hh+1,ii) = gmmLP.bhat(1) + norminv(1-alpha2/2,0,1)*(gmmLP.varbhat(1,1))^0.5;
        IRFslower2_LPIV(hh+1,ii) = gmmLP.bhat(1) - norminv(1-alpha2/2,0,1)*(gmmLP.varbhat(1,1))^0.5;
    end
end

%% Figure
time = (0:horizonLP)';

% standardize
IRFs_LPIV = IRFs_LPIV*shockSize;
IRFsupper_LPIV = IRFsupper_LPIV*shockSize;
IRFslower_LPIV = IRFslower_LPIV*shockSize;
IRFs_LPIV_large = (IRFs_LPIV+IRFs_LPIV_non_lin)*shockSize;

IRFsupper2_LPIV = IRFsupper2_LPIV*shockSize;
IRFslower2_LPIV = IRFslower2_LPIV*shockSize;

figure('Position',[100 100 1000 600],'PaperPositionMode','Auto','DefaultAxesFontSize',13)
signIRFs = 1;
for j=1:nvar % variable
    subplot(2,ceil(nvar/2),j)
    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper_LPIV(1,j); IRFslower_LPIV(1:end,j); flipud([IRFsupper_LPIV(1:end,j); IRFslower_LPIV(end,j)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.2);
    set(hh,'edgecolor','none');
    hold on;

    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper2_LPIV(1,j); IRFslower2_LPIV(1:end,j); flipud([IRFsupper2_LPIV(1:end,j); IRFslower2_LPIV(end,j)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.4);
    set(hh,'edgecolor','none');

    if includeBase
        load('IRFsbench.mat')
        p2=plot(time, signIRFs*IRFs_base(1:horizonLP+1,j), 'Linewidth', 1.5,'Color',[0.8500, 0.3250, 0.0980]);
        plot(time, signIRFs*IRFsupper_base(1:horizonLP+1,j), 'Color',[0.8500, 0.3250, 0.0980],'LineStyle','--');
        plot(time, signIRFs*IRFslower_base(1:horizonLP+1,j),'Color',[0.8500, 0.3250, 0.0980],'LineStyle','--');
        plot(time, signIRFs*IRFsupper2_base(1:horizonLP+1,j),'Color',[0.8500, 0.3250, 0.0980],'LineStyle',':');
        plot(time, signIRFs*IRFslower2_base(1:horizonLP+1,j),'Color',[0.8500, 0.3250, 0.0980],'LineStyle',':');
    end

    p1=plot(time, signIRFs*IRFs_LPIV(:,j),'k', 'Linewidth', 1.5); 
    hold on
    plot(time, signIRFs*IRFs_LPIV_large(:,j),'r', 'Linewidth', 1.5);

    grid on ;hold off;
    if j==1
        ylim([-20 40])
        legend([p1],'LP-IV','Proxy-VAR','AutoUpdate','off')
    elseif  j==4 || j==5
        ylim([-3 2])
    elseif j==2
        ylim([-3 2])
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
    print('-dpdf', gcf, '../../results/appendix/figurea4a');  
end

if saveFigs
    IRFs_LPbase = IRFs_LPIV;
    IRFsupper_LPbase = IRFsupper_LPIV;
    IRFslower_LPbase = IRFslower_LPIV;
    IRFsupper2_LPbase = IRFsupper2_LPIV;
    IRFslower2_LPbase = IRFslower2_LPIV;
    save('../../results/appendix/IRFsbenchLP', 'IRFs_LPbase', 'IRFsupper_LPbase', 'IRFslower_LPbase', 'IRFsupper2_LPbase', 'IRFslower2_LPbase')
end
