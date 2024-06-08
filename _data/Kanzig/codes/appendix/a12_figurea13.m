%% Replication files for "The macroeconomic effects of oil supply news"
% This file generates figure A.13

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

% have to switch up order
data = data(:,[2 1 3:end]);

% number of variables in VAR
nvar = size(data,2);  

% names for paper
varNames_paper = {'World oil production','Real oil price','World oil inventories','World industrial production','U.S. industrial production','U.S. CPI'};
varNames_paperVD = {'Oil production','Real oil price','Oil inventories','World IP','U.S. IP','U.S. CPI'};

% select sample
smplStartInd = find(strcmp(sampleDates,smplStart));
smplEndInd   = find(strcmp(sampleDates,smplEnd));

data = data(smplStartInd:smplEndInd,:);
dataExo = dataExo(smplStartInd:smplEndInd,:);
sampleDates = sampleDates(smplStartInd:smplEndInd,:);
sampleDatesNum = sampleDatesNum(smplStartInd:smplEndInd,:);
    
%% Read in instruments

% get oil supply surprise series
smplStartProxyVARInd = find(strcmp(sampleDates,smplStartProxy));
smplEndProxyVARInd   = find(strcmp(sampleDates,smplEndProxy));

load('../../instrument/OilSurprisesMLog')

proxyRaw1 = [oilProxiesWTIM(:,ncontract)]; 

smplStartProxy1Ind = find(strcmp(sampleDatesProxy,smplStartProxy));
smplEndProxy1Ind   = find(strcmp(sampleDatesProxy,smplEndProxy));

proxyRaw1s = proxyRaw1(smplStartProxy1Ind:smplEndProxy1Ind,:);
proxy1s = proxyRaw1s;

% get oil supply shock proxy from Kilian
load('../../instrument/appendix/KilianInstruments')

proxyRaw2 = [oilProxyKilian(:,2)]; 

smplStartProxy2Ind = find(strcmp(sampleDatesProxyKilian,smplStartProxy));
smplEndProxy2Ind   = find(strcmp(sampleDatesProxyKilian,smplEndProxy));

proxyRaw2s = proxyRaw2(smplStartProxy2Ind:smplEndProxy2Ind,:);
proxy2s    = proxyRaw2s ;

proxy = [proxy2s proxy1s];
[T,np] = size(proxy);
k = 2; % index of variable(s) to be instrumented


%% External instruments VAR with 2 instruments

% run reduced-form VAR 
varEst = varxest(data,dataExo,p);

% only use proxy sample for identification
nexo = size(dataExo,2);
U = varEst.U(smplStartProxyVARInd-p:smplEndProxyVARInd-p,:);    % loose first p observations in estimation  
Sigma = U'*U/(T-p*nvar-nexo);

uhat = [];
for ik = 1:k
    inam = strcat('r',num2str(ik));
    olsEst.(inam) = olsest(proxy,U(:,ik),true,true);
    
    if verbo
        fprintf('F-stat: %4.3f, p-value: %4.3f, F-stat (robust): %4.3f, p-value: %4.3f, R^2: %4.3f, R^2 (adj): %4.3f \n',olsEst.(inam).F,olsEst.(inam).Fpval,olsEst.(inam).Frobust,olsEst.(inam).Frobustpval,olsEst.(inam).R2,olsEst.(inam).R2adj)
    end
    uhat = [uhat olsEst.(inam).yhat];
end

% Get structural impact matrix
b21ib11_2SLS   =   [ones(length(proxy),1) uhat]\U(:,k+1:end);   
b21ib11 = b21ib11_2SLS(2:end,:)';      % 2 SLS coefficients   
Sig11   = Sigma(1:k,1:k);
Sig21   = Sigma(k+1:nvar,1:k);
Sig22   = Sigma(k+1:nvar,k+1:nvar);
ZZp     = b21ib11*Sig11*b21ib11'-(Sig21*b21ib11'+b21ib11*Sig21')+Sig22;
b12b12p = (Sig21- b21ib11*Sig11)'*(ZZp\(Sig21- b21ib11*Sig11));
b11b11p = Sig11-b12b12p;
b22b22p = Sig22+b21ib11*(b12b12p-Sig11)*b21ib11';
L1      = chol(b11b11p,'lower');

b1 = [L1;b21ib11*L1];

% compute IRFs to shock
IRFs_proxy(:,:,1) = varirfsingle(varEst.B(:,1+nexo:end),b1(:,1),p,horizon);
IRFs_proxy(:,:,2) = varirfsingle(varEst.B(:,1+nexo:end),b1(:,2),p,horizon);

impactSurprise = IRFs_proxy(1,1,1);
impactNews = IRFs_proxy(1,2,2);

% standardize
IRFs_proxy(:,:,1) = -IRFs_proxy(:,:,1)./impactSurprise;
IRFs_proxy(:,:,2) = IRFs_proxy(:,:,2)./impactNews*10;


% compute FEVDs to shock
if compFEVDs
    varEst.Sigma = Sigma;   % use Sigma based on truncated sample
    FEVDs_proxy(:,:,1) = varxfevdsingle(varEst,b1(:,1),horizon+1); 
    FEVDs_proxy(:,:,2) = varxfevdsingle(varEst,b1(:,2),horizon+1); 
end

% compute the confidence bands using bootstrapping
bootIRFs = nan(horizon+1,nvar,k,nsim);
bootFEVDs = nan(horizon+1,nvar,k,nsim);
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

    %center the bootstrapped proxy variables
    Proxycentering = zeros(BlockSize,size(proxyLong,2));
    for j = 1:BlockSize 
        subProxy = proxyLong(j:T_est-BlockSize+j,:);
        % account for non-zero mean instrument:
        Proxycentering(j,:) = [mean(subProxy((subProxy(:,1) ~= 0),1),1)-mean(proxyLong((proxyLong(:,1) ~= 0),1),1) ...
                               mean(subProxy((subProxy(:,2) ~= 0),2),1)-mean(proxyLong((proxyLong(:,2) ~= 0),2),1)];
    end
    Proxycentering = repmat(Proxycentering,[nBlock,1]);
    Proxycentering = Proxycentering(1:T_est,:);
  
end

j = 1;
while j <= nsim
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
        for kk = 1:nBlock
            bootProxy(1+BlockSize*(kk-1):BlockSize*kk,:) = ProxyBlocks(:,:,index(kk,1));
        end
        bootProxy = bootProxy(1:T_est,:);

        %center the bootstrapped residuals and proxies
        bootU = bootU - VARcentering;
        for kk = 1:size(proxy,2)
            bootProxy((bootProxy(:,kk)~=0),kk) =...
                bootProxy((bootProxy(:,kk)~=0),kk) - Proxycentering((bootProxy(:,kk)~=0),kk);
        end
        
        % adjust for identification sample
        bootProxy = bootProxy(smplStartProxyVARInd-p:smplEndProxyVARInd-p,:); 
        
        % count the number proxy variables not censored to zero
        ProxyCount(j,:) = sum(abs(bootProxy) > 0,1);
        
        if ~all(ProxyCount(j,:)>=5)
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
    bootPhib    = [ones(length(proxy),1) bootProxy]\bootU;   
    bootPhib    = bootPhib(2:end,:);
    bootPhib11  = bootPhib(1:k,1:k);
    bootPhib21  = bootPhib(1:k,k+1:nvar);
    bootb21ib11 = (bootPhib11\bootPhib21)';                     % 2SLS coefficient

    bootSig11   = bootSigma(1:k,1:k);
    bootSig21   = bootSigma(k+1:nvar,1:k);
    bootSig22   = bootSigma(k+1:nvar,k+1:nvar);
    bootZZp     = bootb21ib11*bootSig11*bootb21ib11'-(bootSig21*bootb21ib11'+bootb21ib11*bootSig21')+bootSig22;
    bootb12b12p = (bootSig21- bootb21ib11*bootSig11)'*(bootZZp\(bootSig21- bootb21ib11*bootSig11));
    bootb11b11p = bootSig11-bootb12b12p;
    bootb22b22p = bootSig22+bootb21ib11*(bootb12b12p-bootSig11)*bootb21ib11';
    bootL1      = chol(bootb11b11p,'lower');

    bootb1 = [bootL1;bootb21ib11*bootL1];

    % compute IRFs
    bootIRFs(:,:,1,j) = varirfsingle(bootvarEst.B(:,1+nexo:end),bootb1(:,1),p,horizon);
    bootIRFs(:,:,2,j) = varirfsingle(bootvarEst.B(:,1+nexo:end),bootb1(:,2),p,horizon);

    % standardize
    bootIRFs(:,:,2,j) = bootIRFs(:,:,2,j)./bootIRFs(1,2,2,j)*10;
        
    if compFEVDs
        bootvarEst.Sigma = bootSigma;
        bootFEVDs(:,:,1,j) = varxfevdsingle(bootvarEst,bootb1(:,1),horizon+1); 
        bootFEVDs(:,:,2,j) = varxfevdsingle(bootvarEst,bootb1(:,2),horizon+1); 
    end
    
    j = j+1;
end

% rescale bootstrapped IRFs (center around sample estimates) and get quantiles
IRFsmed_proxy = quantile(bootIRFs, 0.5, 4);

IRFsupper_proxy = quantile(bootIRFs, 1-alpha/2, 4)-IRFsmed_proxy+IRFs_proxy;  % rescaling does not affect the ordering, thus it is fine to take quantile first
IRFslower_proxy = quantile(bootIRFs, alpha/2, 4)-IRFsmed_proxy+IRFs_proxy;

IRFsupper2_proxy = quantile(bootIRFs, 1-alpha2/2, 4)-IRFsmed_proxy+IRFs_proxy;
IRFslower2_proxy = quantile(bootIRFs, alpha2/2, 4)-IRFsmed_proxy+IRFs_proxy;

% rescale bootstrapped FEVDs 
if compFEVDs
    FEVDslogit = Logit(FEVDs_proxy);
    bootFEVDslogit = Logit(bootFEVDs);
    FEVDsmedlogit = quantile(bootFEVDslogit, 0.5, 4);

    FEVDsupper_proxy = InverseLogit(quantile(bootFEVDslogit, 1-alpha/2, 4)-FEVDsmedlogit+FEVDslogit); 
    FEVDslower_proxy = InverseLogit(quantile(bootFEVDslogit, alpha/2, 4)-FEVDsmedlogit+FEVDslogit);  
end

%% Figure
time = (0:horizon)';
shockNames = {'Oil supply surprise shock', 'Oil supply news shock'};
figure('Position',[100 100 1050 550],'PaperPositionMode','Auto','DefaultAxesFontSize',13);

for ik = 1:k
    j2 = 1;
    for j = [2 1 3 4] 
        subplot(k,nvar-2,(ik-1)*(nvar-2)+j2); 
        hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper_proxy(1,j,ik); IRFslower_proxy(1:end,j,ik); flipud([IRFsupper_proxy(1:end,j,ik); IRFslower_proxy(end,j,ik)])],[0.1, 0.4470, 0.7410]); 
        set(hh,'facealpha',.2);
        set(hh,'edgecolor','none'); 

        hold on;
        hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper2_proxy(1,j,ik); IRFslower2_proxy(1:end,j,ik); flipud([IRFsupper2_proxy(1:end,j,ik); IRFslower2_proxy(end,j,ik)])],[0.1, 0.4470, 0.7410]); 
        set(hh,'facealpha',.4);
        set(hh,'edgecolor','none');

        if ik==2 && includeBase
            load('../../results/IRFsbench')
            p2=plot(time, IRFs_base(:,j2), 'Linewidth', 2,'Color',[0.8500, 0.3250, 0.0980],'LineStyle',':');
        end
        p1=plot(time, IRFs_proxy(:,j,ik),'k', 'Linewidth', 1.5); hold on;

        if ~ismember(0,get(gca,'ylim'))
            line(get(gca,'xlim'),[0 0],'Color','k')
        end
        grid on ;hold off;
        title(varNames_paper{j}) 
        xlabel('Horizon');
        if j2==1
            ylabel(shockNames{ik});
        end
        xlim([0,horizon]);
        xticks([0:10:horizon]);
        if j2==1 && ik==2 && includeBase
            ylim([-5 16])
            legend([p1 p2],{'2-shock model','Baseline'},'Location','Southwest')
        end
        j2= j2+1;
    end
end
xlabel({'Horizon',' ',' '});
pause(0.001)
h=axes('Position',[0.25,0,.675,.5],'Xlim',[0 1],'Ylim',[0 1]);
set(h,'Visible','off');
tightfig;
if k==2
    string_1stage1 = ['First stage regression: F: ',num2str(olsEst.r1.F,' %2.2f'),', robust F: ',num2str(olsEst.r1.Frobust,' %2.2f'),', $R^2$: ',num2str(olsEst.r1.R2*100,' %1.2f'),'\%, Adjusted $R^2$: ',num2str(olsEst.r1.R2adj*100,' %1.2f'),'\%'];
    string_1stage2 = ['First stage regression: F: ',num2str(olsEst.r2.F,' %2.2f'),', robust F: ',num2str(olsEst.r2.Frobust,' %2.2f'),', $R^2$: ',num2str(olsEst.r2.R2*100,' %1.2f'),'\%, Adjusted $R^2$: ',num2str(olsEst.r2.R2adj*100,' %1.2f'),'\%'];
    text('Position',[-0.05 1.025],'string',string_1stage1,'FontSize',14);
    text('Position',[-0.05 0.05],'string',string_1stage2,'FontSize',14);
end
if saveFigs
    print('-dpdf', gcf, '../../results/appendix/figurea13');  
end

