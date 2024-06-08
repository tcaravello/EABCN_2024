%% Replication files for "The macroeconomic effects of oil supply news"
% This file creates figure A.2 Panel B in the appendix

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

%% Heteroskedasticity identified VAR as in Rigobon 2003

% load the proxy
loadProxy;

thresh = 7;
proxy(abs(proxy)>thresh)=0;  % drop extreme values as a check
statementMind(abs(proxy)>thresh)=0;

% External instruments VAR
runProxyVAR;

IRFs_base = IRFs_proxy;
IRFsupper_base = IRFsupper_proxy;
IRFslower_base = IRFslower_proxy;
IRFsupper2_base = IRFsupper2_proxy;
IRFslower2_base = IRFslower2_proxy;


% Heteroskedasticity-based VAR
% load the control series
load('../../instrument/OilSurprisesMLogControl')

proxyControlRaw = [oilPlaceboWTIM(:,ncontract)]; 

proxyControl = proxyControlRaw(smplStartProxyInd:smplEndProxyInd,:);   % we loose the first p values in the VAR


% run reduced-form VAR 
varEst = varxest(data,dataExo,p);

% identification using the covariance structure between proxy and
% reduced-form residuals a la Mertens and Ravn (2013)

% only use proxy sample for identification (potentially a subset of the
% estimation sample)
nexo = size(dataExo,2);
U = varEst.U(smplStartProxyVARInd-p:smplEndProxyVARInd-p,:);    % loose first p observations in estimation  

% Split data into treatment and control sample
statementMindSel = statementMind(smplStartProxyInd:smplEndProxyInd,1);
statementMindPlaceboSel = placeboMind(smplStartProxyInd:smplEndProxyInd,1);
indsR1 = logical(statementMindSel);
indsR2 = logical(statementMindPlaceboSel);

% alternatively compute using IV:
T_OPEC = size(proxy(indsR1,:),1);
T_Control = size(proxyControl(indsR2,:),1);

XrIV = [(proxy(indsR1,:)-mean(proxy(indsR1,:)))/sqrt(T_OPEC); (proxyControl(indsR2,:)-mean(proxyControl(indsR2,:)))/sqrt(T_Control)];
ZrIV = [(proxy(indsR1,:)-mean(proxy(indsR1,:)))/sqrt(T_OPEC); -(proxyControl(indsR2,:)-mean(proxyControl(indsR2,:)))/sqrt(T_Control)];
yiIV = [(U(indsR1,:)-mean(U(indsR1,:)))/sqrt(T_OPEC); (U(indsR2,:)-mean(U(indsR2,:)))/sqrt(T_Control)];

% first stage
olsEst = olsest(ZrIV,XrIV,true,true);
uhat = olsEst.yhat;
            
% second stage
b21ib11_2SLS    =   [uhat]\yiIV;  
b1 = b21ib11_2SLS';      % 2 SLS coefficients               

% compute IRFs to shock
IRFs_proxy = varirfsingle(varEst.B(:,1+nexo:end),b1,p,horizon);
if strcmp(shockType,'custom')
    IRFs_proxy = IRFs_proxy./IRFs_proxy(1,1)*shockSize;
end


% compute the confidence bands using bootstrapping
bootIRFs = nan(horizon+1,nvar,nsim);
bootb1s = nan(nvar,nsim);
bootShocks = nan(T,nsim);
bootDatas = zeros(T+p,nvar,nsim); 

ProxyCount = zeros(nsim,size(proxy,2));
T_est = varEst.T; % length of estimation sample

if strcmp(bootType,'mbb1block')
    % if identification sample is shorter that estimation sample, censor
    % unobserved values to zero
    proxyLong = zeros(T_est, size(proxy,2));
    proxyLong(smplStartProxyVARInd-p:smplEndProxyVARInd-p,:) =  proxy;

    proxyControlLong = zeros(T_est, size(proxyControl,2));
    proxyControlLong(smplStartProxyVARInd-p:smplEndProxyVARInd-p,:) =  proxyControl;
    
    indsR1Long = zeros(T_est, size(indsR1,2));
    indsR1Long(smplStartProxyVARInd-p:smplEndProxyVARInd-p,:) =  indsR1;
    indsR2Long = zeros(T_est, size(indsR2,2));
    indsR2Long(smplStartProxyVARInd-p:smplEndProxyVARInd-p,:) =  indsR2;
    
    BlockSize = round(5.03*T_est^0.25);
    nBlock = ceil(T_est/BlockSize);
    VARBlocks = zeros(BlockSize,nvar,T_est-BlockSize+1);
    ProxyBlocks = zeros(BlockSize,size(proxyLong,2),T_est-BlockSize+1);
    ProxyControlBlocks = zeros(BlockSize,size(proxyLong,2),T_est-BlockSize+1);
    indsR1Blocks = zeros(BlockSize,size(indsR1Long,2),T_est-BlockSize+1);
    indsR2Blocks = zeros(BlockSize,size(indsR2Long,2),T_est-BlockSize+1);
    for j = 1:T_est-BlockSize+1
        VARBlocks(:,:,j) = varEst.U(j:BlockSize+j-1,:);
        ProxyBlocks(:,:,j) = proxyLong(j:BlockSize+j-1,:);
        ProxyControlBlocks(:,:,j) = proxyControlLong(j:BlockSize+j-1,:);
        indsR1Blocks(:,:,j) = indsR1Long(j:BlockSize+j-1,:);
        indsR2Blocks(:,:,j) = indsR2Long(j:BlockSize+j-1,:);
    end

    % center the bootstrapped VAR errors
    VARcentering = zeros(BlockSize,nvar);
    for j = 1:BlockSize
        VARcentering(j,:) = mean(varEst.U(j:T_est-BlockSize+j,:),1);
    end
    VARcentering = repmat(VARcentering,[nBlock,1]);
    VARcentering = VARcentering(1:T_est,:);

    %center the bootstrapped proxy variables
    Proxycentering = zeros(BlockSize,size(proxyLong,2));
    ProxyControlcentering = zeros(BlockSize,size(proxyControlLong,2));
    for j = 1:BlockSize 
        subProxy = proxyLong(j:T_est-BlockSize+j,:);
        subProxyControl = proxyControlLong(j:T_est-BlockSize+j,:);
        %Proxycentering(j,:) = mean(subProxy((subProxy(:,1) ~= 0),1),1); 
        % account for non-zero mean instrument:
        Proxycentering(j,:) = mean(subProxy((subProxy(:,1) ~= 0),1),1) - mean(proxyLong((proxyLong(:,1) ~= 0),1),1);
        ProxyControlcentering(j,:) = mean(subProxyControl((subProxyControl(:,1) ~= 0),1),1) - mean(proxyControlLong((proxyControlLong(:,1) ~= 0),1),1);

    end
    Proxycentering = repmat(Proxycentering,[nBlock,1]);
    Proxycentering = Proxycentering(1:T_est,:);
    
    ProxyControlcentering = repmat(ProxyControlcentering,[nBlock,1]);
    ProxyControlcentering = ProxyControlcentering(1:T_est,:);
    
end

j = 1;
while j <= nsim
    % generate artificial data
    
    if strcmp(bootType,'mbb1block')
        % Moving block bootstrap (Lundsford and Jentsch) using one block
        % type
        
        %draw bootstrapped residuals and proxies
        index = ceil((T_est - BlockSize + 1)*rand(nBlock,1));
        bootU = zeros(nBlock*BlockSize,nvar);
        for kk = 1:nBlock
            bootU(1+BlockSize*(kk-1):BlockSize*kk,:) = VARBlocks(:,:,index(kk,1));
        end
        bootU = bootU(1:T_est,:);
        
        bootProxy = zeros(nBlock*BlockSize,size(proxy,2));
        bootProxyControl = zeros(nBlock*BlockSize,size(proxyControl,2));
        for kk = 1:nBlock
            bootProxy(1+BlockSize*(kk-1):BlockSize*kk,:) = ProxyBlocks(:,:,index(kk,1));
            bootProxyControl(1+BlockSize*(kk-1):BlockSize*kk,:) = ProxyControlBlocks(:,:,index(kk,1));
        end
        bootProxy = bootProxy(1:T_est,:);
        bootProxyControl = bootProxyControl(1:T_est,:);
        
        %center the bootstrapped residuals and proxies
        bootU = bootU - VARcentering;
        for kk = 1:size(proxy,2)
            bootProxy((bootProxy(:,kk)~=0),kk) =...
                bootProxy((bootProxy(:,kk)~=0),kk) - Proxycentering((bootProxy(:,kk)~=0),kk);
            bootProxyControl((bootProxyControl(:,kk)~=0),kk) =...
                bootProxyControl((bootProxyControl(:,kk)~=0),kk) - ProxyControlcentering((bootProxyControl(:,kk)~=0),kk);
        end
        
        % adjust for identification sample
        bootProxy = bootProxy(smplStartProxyVARInd-p:smplEndProxyVARInd-p,:); 
        bootProxyControl = bootProxyControl(smplStartProxyVARInd-p:smplEndProxyVARInd-p,:); 
        
        % treatment months
        bootindsR1 = zeros(nBlock*BlockSize,size(indsR1,2)); 
        bootindsR2 = zeros(nBlock*BlockSize,size(indsR2,2)); 
        for kk = 1:nBlock
            bootindsR1(1+BlockSize*(kk-1):BlockSize*kk,:) = indsR1Blocks(:,:,index(kk,1));
            bootindsR2(1+BlockSize*(kk-1):BlockSize*kk,:) = indsR2Blocks(:,:,index(kk,1));
        end
        bootindsR1 = bootindsR1(1:T_est,:);
        bootindsR2 = bootindsR2(1:T_est,:);
        
        % adjust for identification sample
        bootindsR1 = logical(bootindsR1(smplStartProxyVARInd-p:smplEndProxyVARInd-p,:)); 
        bootindsR2 = logical(bootindsR2(smplStartProxyVARInd-p:smplEndProxyVARInd-p,:)); 
        
        % count the number proxy variables not censored to zero
        ProxyCount(j,:) = sum(abs(bootProxy) > 0,1);
        
        if ProxyCount(j,:)<15
            continue
        end
        
        % simulate VAR starting from initial values
        Xexo = varEst.Xexo;
        
        bootData = zeros(T_est+p,nvar); 
        bootData(1:p,:) = data(1:p,:); %initial values of y, same for all j
        for i = p+1:T_est+p
            bootData(i,:)= varEst.B*[Xexo(i-p,:)'; vec(fliplr(bootData(i-p:i-1,:)'))] ...
                             + bootU(i-p,:)'; % bootstrap
        end
        
        % re-estimate the VAR
        bootvarEst = varxest(bootData,dataExo,p);
        
        % only use proxy sample for identification
        bootU = bootvarEst.U(smplStartProxyVARInd-p:smplEndProxyVARInd-p,:);    % loose first p observations in estimation  
        bootSigma = bootU'*bootU/(T-p*nvar-nexo);
    end
        
    % structural impact matrix
    % 2SLS of U1 on U2 using proxy as an instrument (include a constant for the case that proxy is not demeaned)

    bootT_OPEC = size(bootProxy(bootindsR1,:),1);
    bootT_Control = size(bootProxyControl(bootindsR2,:),1);

    bootXrIV = [(bootProxy(bootindsR1,:)-mean(bootProxy(bootindsR1,:)))/sqrt(bootT_OPEC); (bootProxyControl(bootindsR2,:)-mean(bootProxyControl(bootindsR2,:)))/sqrt(bootT_Control)];
    bootZrIV = [(bootProxy(bootindsR1,:)-mean(bootProxy(bootindsR1,:)))/sqrt(bootT_OPEC); -(bootProxyControl(bootindsR2,:)-mean(bootProxyControl(bootindsR2,:)))/sqrt(bootT_Control)];
    bootyiIV = [(bootU(bootindsR1,:)-mean(bootU(bootindsR1,:)))/sqrt(bootT_OPEC); (bootU(bootindsR2,:)-mean(bootU(bootindsR2,:)))/sqrt(bootT_Control)];

    % first stage
    bootOlsEst = olsest(bootZrIV,bootXrIV,true,true);
    bootuhat = bootOlsEst.yhat;

    % second stage
    bootb21ib11_2SLS = [bootuhat]\bootyiIV;  
    bootb1 = bootb21ib11_2SLS';      % 2 SLS coefficients

    % compute IRFs
    bootIRFs(:,:,j)  = varirfsingle(bootvarEst.B(:,1+nexo:end),bootb1,p,horizon);
    if strcmp(shockType,'custom')
        bootIRFs(:,:,j) = bootIRFs(:,:,j)./bootIRFs(1,1,j)*shockSize;
    end 
    bootb1s(:,j) = bootb1;
    bootSigma1 = bootU(bootindsR1,:)'*bootU(bootindsR1,:)/(sum(bootindsR1)); % -p*nvar-nexo
    bootShocks(:,j) = (bootb1'*inv2(bootSigma1)*bootU')'*inv2(bootb1'*inv2(bootSigma1)*bootb1);
    bootDatas(:,:,j) = bootData;
  
    j = j+1;
end

IRFsmed   = quantile(bootIRFs, 0.5, 3);
IRFslower_proxy = quantile(bootIRFs, 1-alpha/2, 3)-IRFsmed+IRFs_proxy;  % correct for small-sample bias
IRFsupper_proxy = quantile(bootIRFs, alpha/2, 3)-IRFsmed+IRFs_proxy;
IRFslower2_proxy = quantile(bootIRFs, 1-alpha2/2, 3)-IRFsmed+IRFs_proxy;  % correct for small-sample bias
IRFsupper2_proxy = quantile(bootIRFs, alpha2/2, 3)-IRFsmed+IRFs_proxy;

%% Figure
time = (0:horizon)';    % time horizon for IRFs 

signIRFs = 1/IRFs_proxy(1,1)*shockSize;
figure('Position',[100 100 1000 600],'PaperPositionMode','Auto','DefaultAxesFontSize',13);
for j=1:nvar %variable
    subplot(ceil(nvar/3),3,j); 
    hold on;

    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[signIRFs*IRFsupper_proxy(1,j); signIRFs*IRFslower_proxy(1:end,j); flipud([signIRFs*IRFsupper_proxy(1:end,j); signIRFs*IRFslower_proxy(end,j)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.2);
    set(hh,'edgecolor','none');
    hold on;
    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[signIRFs*IRFsupper2_proxy(1,j); signIRFs*IRFslower2_proxy(1:end,j); flipud([signIRFs*IRFsupper2_proxy(1:end,j); signIRFs*IRFslower2_proxy(end,j)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.4);
    set(hh,'edgecolor','none');
    if includeBase
        p1=plot(time, IRFs_base(:,j), 'Linewidth', 1.5,'Color',[0.8500, 0.3250, 0.0980]);
        plot(time, IRFsupper_base(:,j), 'Color',[0.8500, 0.3250, 0.0980],'Linestyle','--');
        plot(time, IRFslower_base(:,j), 'Color',[0.8500, 0.3250, 0.0980],'Linestyle','--');
        plot(time, IRFsupper2_base(:,j), 'Color',[0.8500, 0.3250, 0.0980],'Linestyle',':');
        plot(time, IRFslower2_base(:,j), 'Color',[0.8500, 0.3250, 0.0980],'Linestyle',':');
    end
    p2=plot(time, signIRFs*IRFs_proxy(:,j),'k', 'Linewidth', 1.5); 
    if ~ismember(0,get(gca,'ylim'))
        line(get(gca,'xlim'),[0 0],'Color','k')
    end
    box on
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
       legend([p2,p1],{'Heteroskedasticity','External instrument'})
    end
end
tightfig;
if saveFigs
    print('-dpdf', gcf, '../../results/appendix/figurea2b');  
end

