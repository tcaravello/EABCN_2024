%% Replication files for "The macroeconomic effects of oil supply news"
% This file generates figure A.25

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
smplStart = '1974M02'; 
smplEnd   = '2014M12'; 

% Instrument sample
% range has to be contained in estimation sample
smplStartProxy = '1976M02'; 
smplEndProxy   = '2014M12'; 

% VAR specifics
p          = 24;       
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
load('../../data/appendix/dataKilianM')
% data: transformed endogenous variables 
% dataExo: exogenous variables (e.g. constant, trend)
% sampleDates: sample dates (string format)
% sampleDatesNum: sample dates (numeric format, e.g. 2000 = 2000M1)
% varNames: labels of variables

% number of variables in VAR
nvar = size(data,2);  

% names for paper
varNames_paper = {'Real oil price','World oil production','World oil inventories','Global activity'};
varNames_paperVD = varNames_paper;
diffInd     = [0 1 1 0];

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

% run reduced-form VAR 
varEst = varxest(data,dataExo,p);

% only use proxy sample for identification (potentially a subset of the
% estimation sample)
nexo = size(dataExo,2);
U = varEst.U(smplStartProxyVARInd-p:smplEndProxyVARInd-p,:);    % loose first p observations in estimation  
Sigma = U'*U/(T-p*nvar-nexo);

% first stage
olsEst = olsest(proxy,U(:,k),true,true);

if verbo
    fprintf('F-stat: %4.3f, p-value: %4.3f, F-stat (robust): %4.3f, p-value: %4.3f, R^2: %4.3f, R^2 (adj): %4.3f \n',olsEst.F,olsEst.Fpval,olsEst.Frobust,olsEst.Frobustpval,olsEst.R2,olsEst.R2adj)
end
uhat = olsEst.yhat;

% second stage
b21ib11_2SLS    =   [ones(length(proxy),1) uhat]\U(:,k+1:end);  
b21ib11 = b21ib11_2SLS(2:end,:)';      % 2 SLS coefficients               
Sig11   = Sigma(1:k,1:k);
Sig21   = Sigma(k+1:nvar,1:k);
Sig22   = Sigma(k+1:nvar,k+1:nvar);
ZZp     = b21ib11*Sig11*b21ib11'-(Sig21*b21ib11'+b21ib11*Sig21')+Sig22;
b12b12p = (Sig21- b21ib11*Sig11)'*(ZZp\(Sig21- b21ib11*Sig11));
b11b11p = Sig11-b12b12p;
b11     = sqrt(b11b11p);
b1      = [b11; b21ib11*b11];
b1unit  = [1; b21ib11]*shockSize;


% compute IRFs to shock
IRFs_proxy = varirfsingle(varEst.B(:,1+nexo:end),b1,p,horizon,diffInd);
if strcmp(shockType,'custom')
    IRFs_proxy = IRFs_proxy./IRFs_proxy(1,1)*shockSize;
end

% compute FEVDs to shock
varEst.Sigma = Sigma;   
FEVDs_proxy = varxfevdsingle(varEst,b1,horizon+1); 

% compute the confidence bands using bootstrapping
bootIRFs = nan(horizon+1,nvar,nsim);
bootFEVDs = nan(horizon+1,nvar,nsim);
bootFstat = nan(nsim,2);
bootBs = nan(nvar,size(varEst.B,2),nsim);
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

    BlockSize = round(5.03*T_est^0.25);
    nBlock = ceil(T_est/BlockSize);
    VARBlocks = zeros(BlockSize,nvar,T_est-BlockSize+1);
    ProxyBlocks = zeros(BlockSize,size(proxyLong,2),T_est-BlockSize+1);
    for j = 1:T_est-BlockSize+1
        VARBlocks(:,:,j) = varEst.U(j:BlockSize+j-1,:);
        ProxyBlocks(:,:,j) = proxyLong(j:BlockSize+j-1,:);
    end

    % center the bootstrapped VAR errors
    VARcentering = zeros(BlockSize,nvar);
    for j = 1:BlockSize
        VARcentering(j,:) = mean(varEst.U(j:T_est-BlockSize+j,:),1);
    end
    VARcentering = repmat(VARcentering,[nBlock,1]);
    VARcentering = VARcentering(1:T_est,:);

    % center the bootstrapped proxy variables
    Proxycentering = zeros(BlockSize,size(proxyLong,2));
    for j = 1:BlockSize 
        subProxy = proxyLong(j:T_est-BlockSize+j,:);
        % account for non-zero mean instrument:
        Proxycentering(j,:) = mean(subProxy((subProxy(:,1) ~= 0),1),1) - mean(proxyLong((proxyLong(:,1) ~= 0),1),1);
    end
    Proxycentering = repmat(Proxycentering,[nBlock,1]);
    Proxycentering = Proxycentering(1:T_est,:);
     
end

j = 1;
while j <= nsim
    % generate artificial data
    
    if strcmp(bootType,'mbb1block')
        % Moving block bootstrap (Lundsford and Jentsch) using one block
        % type
        
        % draw bootstrapped residuals and proxies
        index = ceil((T_est - BlockSize + 1)*rand(nBlock,1));
        bootU = zeros(nBlock*BlockSize,nvar);
        for kk = 1:nBlock
            bootU(1+BlockSize*(kk-1):BlockSize*kk,:) = VARBlocks(:,:,index(kk,1));
        end
        bootU = bootU(1:T_est,:);
        
        bootProxy = zeros(nBlock*BlockSize,size(proxy,2));
        for kk = 1:nBlock
            bootProxy(1+BlockSize*(kk-1):BlockSize*kk,:) = ProxyBlocks(:,:,index(kk,1));
        end
        bootProxy = bootProxy(1:T_est,:);

        % center the bootstrapped residuals and proxies
        bootU = bootU - VARcentering;
        for kk = 1:size(proxy,2)
            bootProxy((bootProxy(:,kk)~=0),kk) =...
                bootProxy((bootProxy(:,kk)~=0),kk) - Proxycentering((bootProxy(:,kk)~=0),kk);
        end
        
        % adjust for identification sample
        bootProxy = bootProxy(smplStartProxyVARInd-p:smplEndProxyVARInd-p,:); 
        
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

    % first stage
    bootOlsEst = olsest(bootProxy,bootU(:,1),true,true);
    bootuhat = bootOlsEst.yhat;

    bootFstat(j,:) = [bootOlsEst.F bootOlsEst.Frobust];
    
    % second stage
    bootb21ib11_2SLS    =   [ones(length(proxy),1) bootuhat]\bootU(:,k+1:end);  
    bootb21ib11 = bootb21ib11_2SLS(2:end,:)';      % 2 SLS coefficients    
    
    bootSig11   = bootSigma(1:k,1:k);
    bootSig21   = bootSigma(k+1:nvar,1:k);
    bootSig22   = bootSigma(k+1:nvar,k+1:nvar);
    bootZZp     = bootb21ib11*bootSig11*bootb21ib11'-(bootSig21*bootb21ib11'+bootb21ib11*bootSig21')+bootSig22;
    bootb12b12p = (bootSig21- bootb21ib11*bootSig11)'*(bootZZp\(bootSig21- bootb21ib11*bootSig11));
    bootb11b11p = bootSig11-bootb12b12p;
    bootb11     = sqrt(bootb11b11p);
    bootb1      = [bootb11; bootb21ib11*bootb11];

    % compute IRFs
    bootIRFs(:,:,j)  = varirfsingle(bootvarEst.B(:,1+nexo:end),bootb1,p,horizon,diffInd);
    if strcmp(shockType,'custom')
        bootIRFs(:,:,j) = bootIRFs(:,:,j)./bootIRFs(1,1,j)*shockSize;
    end 
    if compFEVDs
        bootvarEst.Sigma = bootSigma;
        bootFEVDs(:,:,j) = varxfevdsingle(bootvarEst,bootb1,horizon+1); 
    end
    bootBs(:,:,j) = bootvarEst.B;
    bootb1s(:,j) = bootb1;
    bootShocks(:,j) = (bootb1'*inv2(bootSigma)*bootU')';
    bootDatas(:,:,j) = bootData;
    
    j = j+1;
end


% compute percentile interval and center around point estimate
IRFsmed_proxy = quantile(bootIRFs, 0.5, 3);

IRFsupper_proxy = quantile(bootIRFs, 1-alpha/2, 3)-IRFsmed_proxy+IRFs_proxy;  
IRFslower_proxy = quantile(bootIRFs, alpha/2, 3)-IRFsmed_proxy+IRFs_proxy;

IRFsupper2_proxy = quantile(bootIRFs, 1-alpha2/2, 3)-IRFsmed_proxy+IRFs_proxy;  
IRFslower2_proxy = quantile(bootIRFs, alpha2/2, 3)-IRFsmed_proxy+IRFs_proxy;

% Compute percentile intervals. To make sure that bands lie within [0,1], we
% use a logistic transformation

if compFEVDs
    FEVDslogit = Logit(FEVDs_proxy);
    bootFEVDslogit = Logit(bootFEVDs);
    FEVDsmedlogit = quantile(bootFEVDslogit, 0.5, 3);

    FEVDsupper_proxy = InverseLogit(quantile(bootFEVDslogit, 1-alpha/2, 3)-FEVDsmedlogit+FEVDslogit); 
    FEVDslower_proxy = InverseLogit(quantile(bootFEVDslogit, alpha/2, 3)-FEVDsmedlogit+FEVDslogit);  
end
 

%% Figure
time = (0:horizon)';    % time horizon for IRFs and FEVDs

figure('Position',[100 100 700 600],'PaperPositionMode','Auto','DefaultAxesFontSize',13);
signIRFs = 1;
for j=1:nvar % variable
    h(j) = subplot(2,2,j); 
    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper_proxy(1,j); IRFslower_proxy(1:end,j); flipud([IRFsupper_proxy(1:end,j); IRFslower_proxy(end,j)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.2);
    set(hh,'edgecolor','none'); 

    hold on;
    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper2_proxy(1,j); IRFslower2_proxy(1:end,j); flipud([IRFsupper2_proxy(1:end,j); IRFslower2_proxy(end,j)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.4);
    set(hh,'edgecolor','none');

    plot(time, signIRFs*IRFs_proxy(:,j),'k', 'Linewidth', 1.5); hold on;
    if ~ismember(0,get(gca,'ylim'))
        line(get(gca,'xlim'),[0 0],'Color','k')
    end
    grid on ;hold off;
    title(varNames_paper{j}) 
    xlabel('Months');
    ylabel('\%');
    if j==3
        ylabel('MMBBL')
    end
    if j==4
        ylabel(' ')
    end
    xlim([0,horizon]);
end
pause(0.001)
h=axes('Position',[0,0,1,.5],'Xlim',[0 1],'Ylim',[0 1]);
set(h,'Visible','off');
if k==1
    string_1stage = ['First stage regression: F: ',num2str(olsEst.F,' %2.2f'),', robust F: ',num2str(olsEst.Frobust,' %2.2f'),', $R^2$: ',num2str(olsEst.R2*100,' %1.2f'),'\%, Adjusted $R^2$: ',num2str(olsEst.R2adj*100,' %1.2f'),'\%'];
    text('Position',[0.025 0.0],'string',string_1stage,'FontSize',14);
end
tightfig;
if saveFigs
    print('-dpdf', gcf, '../../results/appendix/figurea25');  
end
