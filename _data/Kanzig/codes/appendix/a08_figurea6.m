%% Replication files for "The macroeconomic effects of oil supply news"
% This figure generates Figure A.6

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


%% Compute IRFs based on selection of different models
%% VARs
% External instruments VAR

% load the proxy
loadProxy;

% lag specifications
pSet = [12 18 24];

% Proxy VAR
IRFs_proxys = zeros(horizon+1,nvar,length(pSet));
kk = 1;
for p = pSet
    % run reduced-form VAR 
    varEst = varxest(data,dataExo,p);

    % identification using the covariance structure between proxy and
    % reduced-form residuals a la Mertens and Ravn (2013)

    % only use proxy sample for identification (potentially a subset of the
    % estimation sample)
    nexo = size(dataExo,2);
    U = varEst.U;    % loose first p observations in estimation  
    Sigma = varEst.Sigma;

    % run first stage (k variables to be instrumented with size(proxy,2)
    % proxies. Need that size(proxy,2)>=k)
    uhat = [];
    for ik = 1:k
        inam = strcat('r',num2str(ik));
        olsEst.(inam) = olsest(proxyRaw(p+1:end),U(:,ik),true,true);

        uhat = [uhat olsEst.(inam).yhat];
    end

    % second stage
    b21ib11_2SLS    =   [ones(length(proxyRaw(p+1:end)),1) uhat]\U(:,k+1:end);  
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
    IRFs_proxys(:,:,kk) = varirfsingle(varEst.B(:,1+nexo:end),b1,p,horizon);
    if strcmp(shockType,'custom')
        IRFs_proxys(:,:,kk) = IRFs_proxys(:,:,kk)./IRFs_proxys(1,1,kk)*shockSize;
    end
    
    if kk == 1
        oilSupplyNewsShock = (b1unit'*inv2(Sigma)*U')'*inv2(b1unit'*inv2(Sigma)*b1unit);
    end
    
    kk = kk+1;
end

% Rigobon VAR
IRFs_rigos = zeros(horizon+1,nvar,length(pSet));
kk = 1;
for p = pSet
    proxyRig = proxyRaw(p+1:end);
    proxyRig(abs(proxyRig)<0.5)=0; 

    iregime = proxyRig~=0;

    % run reduced-form VAR 
    varEst = varxest(data,dataExo,p);

    T_est = varEst.T; % length of estimation sample
    nexo = size(dataExo,2); % Number of exog variables
    
    % Partition residuals
    U1 = varEst.U(~iregime,:);
    T1 = length(U1);
    U2 = varEst.U(iregime,:);
    T2 = length(U2);

    % Compute variances
    Sigma1 = U1'*U1/(T1);
    Sigma2 = U2'*U2/(T2);

    % Compute impact matrix by simultaneously diagonalizing Sigma1 and Sigma 2
    A0 = simdiag(Sigma1, Sigma2);
    Lam = A0*Sigma2*A0';
    invA0 = inv(A0);

    % Select shock as shock with maximum impact on price (or alternatively has
    % the maximum lambda?) max(abs((diag(Lam)-1)*100))
    [~,ishock] = max(invA0(1,:));

    % Compute impact vector of MP shock as in Wright
    D = dupmatrix(nvar);
    Dinv = inv2(D'*D)*D';

    VARvechSigma1 = 2*Dinv*(kron(Sigma1,Sigma1))*Dinv'/(T1); %-p*nvar-nexo);
    VARvechSigma2 = 2*Dinv*(kron(Sigma2,Sigma2))*Dinv'/(T2); %-p*nvar-nexo);

    objfun = @(b) (vech(Sigma2-Sigma1) - vech(b*b'))'*inv2(VARvechSigma1 + VARvechSigma2)*(vech(Sigma2-Sigma1) - vech(b*b'));

    b0 = fminunc(objfun,invA0(:,ishock)); 
    [b1,minb] = fminunc(objfun,b0);
    xi = vech(Sigma2-Sigma1) - vech(b1*b1');

    VarStat = vech(Sigma2-Sigma1)'*inv2(VARvechSigma1 + VARvechSigma2)*vech(Sigma2-Sigma1);

    IRFs_rigos(:,:,kk) = varirfsingle(varEst.B(:,1+nexo:end),b1,p,horizon);
    sd = IRFs_rigos(1,1,kk);
    if strcmp(shockType,'custom')
        IRFs_rigos(:,:,kk) = IRFs_rigos(:,:,kk)./sd*shockSize;
    end
    
    kk = kk+1;
end

% Heteroskedsaticity-based VAR à la Nakamura and Steinsson
% load the proxy
loadProxy;

% load placebo
load('../../instrument/OilSurprisesMLogControl')

proxyControlRaw = [oilPlaceboWTIM(:,ncontract)]; 

proxyControl = proxyControlRaw(smplStartProxyInd:smplEndProxyInd,:);   % we loose the first p values in the VAR

IRFs_rigosNaka = zeros(horizon+1,nvar,length(pSet));
kk = 1;
for p = pSet

    % run reduced-form VAR 
    varEst = varxest(data,dataExo,p);

    % only use proxy sample for identification (potentially a subset of the
    % estimation sample)
    nexo = size(dataExo,2);
    U = varEst.U(smplStartProxyVARInd-p:smplEndProxyVARInd-p,:);    % loose first p observations in estimation  
    Sigma = U'*U/(T-p*nvar-nexo);

    % Split data into treatment and control sample
    statementMindSel = statementMind(smplStartProxyInd:smplEndProxyInd,1);
    statementMindPlaceboSel = placeboMind(smplStartProxyInd:smplEndProxyInd,1);
    indsR1 = logical(statementMindSel);
    indsR2 = logical(statementMindPlaceboSel);

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
    IRFs_rigosNaka(:,:,kk) = varirfsingle(varEst.B(:,1+nexo:end),b1,p,horizon);
    if strcmp(shockType,'custom')
        IRFs_rigosNaka(:,:,kk) = IRFs_rigosNaka(:,:,kk)./IRFs_rigosNaka(1,1,kk)*shockSize;
    end
    
    if kk == 1
        Sigma1 = U(indsR1,:)'*U(indsR1,:)/(sum(indsR1)); 

        oilSupplyNewsShock1 = (b1'*inv2(Sigma1)*varEst.U')';
    end
    
    kk = kk+1;
end

%% Local projections

% settings
horizonLP = horizon;

colval = [0.8500, 0.3250, 0.0980];

% LPIVs
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
    
    %standardize
    IRFs_LPIVs(:,:,kk) = IRFs_LPIVs(:,:,kk)*shockSize;
    
    kk = kk+1;
end

% LPIVs flexible controls 
IRFs_LPIVsflex = zeros(horizonLP+1,nvar,length(pSet));
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
            IRFs_LPIVsflex(hh+1,ii,kk) = gmmLP.bhat(1);
        end
    end
    
    %standardize
    IRFs_LPIVsflex(:,:,kk) = IRFs_LPIVsflex(:,:,kk)*shockSize;
    
    kk = kk+1;
end

% Rigobon LP
IRFs_LPRigos = zeros(horizonLP+1,nvar,length(pSet));
kk = 1;
for pLP = pSet
    for hh = 0:horizonLP
        for ii = 1:nvar

            % first clean for lags on longer sample
            yi_pre = data(pLP+1+hh:end-horizonLP+hh,ii);

            Xr_pre = [];
            for jj = 1:pLP
                Xr_pre = [Xr_pre data(pLP+1-jj:end-horizonLP-jj,:)];
            end
            Xr_pre = [Xr_pre dataExo(pLP+1:end-horizonLP,1:end)];

            % construct shorter sample for LP on shock
            yi_all = [nan(pLP,1); yi_pre - Xr_pre*(Xr_pre\yi_pre); nan(horizonLP,1)];
            yi = yi_all(smplStartProxyVARInd:smplEndProxyVARInd-horizonLP);
            Xr = proxyRaw(smplStartProxyInd:smplEndProxyInd-horizonLP);                     
            XrC = proxyControlRaw(smplStartProxyInd:smplEndProxyInd-horizonLP); 

            indsR1 = logical(statementMind(smplStartProxyInd:smplEndProxyInd-horizonLP));
            indsR2 = logical(placeboMind(smplStartProxyInd:smplEndProxyInd-horizonLP));

            % compute using IV:
            T_OPEC = size(Xr(indsR1,:),1);
            T_Control = size(XrC(indsR2,:),1);
            
            XrIV = [(Xr(indsR1,:)-mean(Xr(indsR1,:)))/sqrt(T_OPEC); (XrC(indsR2,:)-mean(XrC(indsR2,:)))/sqrt(T_Control)];
            ZrIV = [(Xr(indsR1,:)-mean(Xr(indsR1,:)))/sqrt(T_OPEC); -(XrC(indsR2,:)-mean(XrC(indsR2,:)))/sqrt(T_Control)];
            yiIV = [(yi(indsR1,:)-mean(yi(indsR1,:)))/sqrt(T_OPEC); (yi(indsR2,:)-mean(yi(indsR2,:)))/sqrt(T_Control)];
            
            gmmLP = gmmest(XrIV,ZrIV,yiIV,false,2,hh+1);
            IRFs_LPRigos(hh+1,ii,kk) = gmmLP.bhat(1);

        end
    end
    
    % standardize
    estSize = IRFs_LPRigos(1,1,kk);
    IRFs_LPRigos(:,:,kk) = IRFs_LPRigos(:,:,kk)./estSize*shockSize;
    
    kk = kk+1;
end

% flexible controls
IRFs_LPRigosflex = zeros(horizonLP+1,nvar,length(pSet));
kk = 1;
for pLP = pSet
    for hh = 0:horizonLP
        for ii = 1:nvar
            
            % first clean for lags on longer sample
            yi_pre = data(pLP+1+hh:end-horizonLP+hh,ii);

            Xr_pre = [];
            for jj = 1:pLP
                Xr_pre = [Xr_pre data(pLP+1-jj:end-horizonLP-jj,[1:4])];
                if ii>4
                    Xr_pre = [Xr_pre data(pLP+1-jj:end-horizonLP-jj,ii)];
                end
            end
            Xr_pre = [Xr_pre dataExo(pLP+1:end-horizonLP,1:end)];

            % construct shorter sample for LP on shock
            yi_all = [nan(pLP,1); yi_pre - Xr_pre*(Xr_pre\yi_pre); nan(horizonLP,1)];
            yi = yi_all(smplStartProxyVARInd:smplEndProxyVARInd-horizonLP);
            Xr = proxyRaw(smplStartProxyInd:smplEndProxyInd-horizonLP);                      
            XrC = proxyControlRaw(smplStartProxyInd:smplEndProxyInd-horizonLP); 

            indsR1 = logical(statementMind(smplStartProxyInd:smplEndProxyInd-horizonLP));
            indsR2 = logical(placeboMind(smplStartProxyInd:smplEndProxyInd-horizonLP));

            % compute using IV:
            T_OPEC = size(Xr(indsR1,:),1);
            T_Control = size(XrC(indsR2,:),1);
            
            XrIV = [(Xr(indsR1,:)-mean(Xr(indsR1,:)))/sqrt(T_OPEC); (XrC(indsR2,:)-mean(XrC(indsR2,:)))/sqrt(T_Control)];
            ZrIV = [(Xr(indsR1,:)-mean(Xr(indsR1,:)))/sqrt(T_OPEC); -(XrC(indsR2,:)-mean(XrC(indsR2,:)))/sqrt(T_Control)];
            yiIV = [(yi(indsR1,:)-mean(yi(indsR1,:)))/sqrt(T_OPEC); (yi(indsR2,:)-mean(yi(indsR2,:)))/sqrt(T_Control)];
            
            gmmLP = gmmest(XrIV,ZrIV,yiIV,false,2,hh+1);
            IRFs_LPRigosflex(hh+1,ii,kk) = gmmLP.bhat(1);

        end
    end
    
    % standardize
    estSize = IRFs_LPRigosflex(1,1,kk);
    IRFs_LPRigosflex(:,:,kk) = IRFs_LPRigosflex(:,:,kk)./estSize*shockSize;
    
    kk = kk+1;
end

% LP on shock
IRFs_LPShocks = zeros(horizonLP+1,nvar,length(pSet));

kk = 1;
for pLP = pSet
    for hh = 0:horizonLP
        for ii = 1:nvar
 
            yi = data(pLP+1+hh:end-horizonLP+hh,ii);

            Xr = oilSupplyNewsShock(1+pLP-12:end-horizonLP);  
            for jj = 1:pLP
                Xr = [Xr data(pLP+1-jj:end-horizonLP-jj,:)];
            end
            Xr = [Xr dataExo(pLP+1:end-horizonLP,2:end)];

            olsLP = olsest(Xr,yi,true,2,hh+1);

            IRFs_LPShocks(hh+1,ii,kk) = olsLP.bhat(1);
         end
    end
    
    % standardize
    estSize = IRFs_LPShocks(1,1,kk);
    IRFs_LPShocks(:,:,kk) = IRFs_LPShocks(:,:,kk)./estSize*shockSize;
    
    kk = kk+1;
end

% flexible controls
IRFs_LPShocksflex = zeros(horizonLP+1,nvar,length(pSet));

kk = 1;
for pLP = pSet
    for hh = 0:horizonLP
        for ii = 1:nvar
 
            yi = data(pLP+1+hh:end-horizonLP+hh,ii);

            Xr = oilSupplyNewsShock(1+pLP-12:end-horizonLP);  
            for jj = 1:pLP
                Xr = [Xr data(pLP+1-jj:end-horizonLP-jj,[1:4])];
                if ii>4
                    Xr = [Xr data(pLP+1-jj:end-horizonLP-jj,ii)];
                end
            end
            Xr = [Xr dataExo(pLP+1:end-horizonLP,2:end)];

            olsLP = olsest(Xr,yi,true,2,hh+1);

            IRFs_LPShocksflex(hh+1,ii,kk) = olsLP.bhat(1);
         end
    end
    
    % standardize
    estSize = IRFs_LPShocksflex(1,1,kk);
    IRFs_LPShocksflex(:,:,kk) = IRFs_LPShocksflex(:,:,kk)./estSize*shockSize;
    
    kk = kk+1;
end


% shock from heteroskedasticity-based VAR
IRFs_LPShockshetero = zeros(horizonLP+1,nvar,length(pSet));

kk = 1;
for pLP = pSet
    for hh = 0:horizonLP
        for ii = 1:nvar

            yi = data(pLP+1+hh:end-horizonLP+hh,ii);

            Xr = oilSupplyNewsShock1(1+pLP-12:end-horizonLP);  
            for jj = 1:pLP
                Xr = [Xr data(pLP+1-jj:end-horizonLP-jj,:)];
            end
            Xr = [Xr dataExo(pLP+1:end-horizonLP,2:end)];

            olsLP = olsest(Xr,yi,true,2,hh+1);

            IRFs_LPShockshetero(hh+1,ii,kk) = olsLP.bhat(1);
         end
    end
    
    % standardize
    estSize = IRFs_LPShockshetero(1,1,kk);
    IRFs_LPShockshetero(:,:,kk) = IRFs_LPShockshetero(:,:,kk)./estSize*shockSize;
    
    kk = kk+1;
end

% flexible controls
IRFs_LPShocksheteroflex = zeros(horizonLP+1,nvar,length(pSet));

kk = 1;
for pLP = pSet
    for hh = 0:horizonLP
        for ii = 1:nvar
                
            yi = data(pLP+1+hh:end-horizonLP+hh,ii);

            Xr = oilSupplyNewsShock1(1+pLP-12:end-horizonLP);  
            for jj = 1:pLP
                Xr = [Xr data(pLP+1-jj:end-horizonLP-jj,[1:4])];
                if ii>4
                    Xr = [Xr data(pLP+1-jj:end-horizonLP-jj,ii)];
                end
            end
            Xr = [Xr dataExo(pLP+1:end-horizonLP,2:end)];
           
            olsLP = olsest(Xr,yi,true,2,hh+1);

            IRFs_LPShocksheteroflex(hh+1,ii,kk) = olsLP.bhat(1);
         end
    end
    
    % standardize
    estSize = IRFs_LPShocksheteroflex(1,1,kk);
    IRFs_LPShocksheteroflex(:,:,kk) = IRFs_LPShocksheteroflex(:,:,kk)./estSize*shockSize;
    
    kk = kk+1;
end


%% Figure
IRFs_All = IRFs_LPIVsflex;
IRFs_All = cat(3,IRFs_All,IRFs_LPIVs);
IRFs_All = cat(3,IRFs_All,IRFs_LPRigos);
IRFs_All = cat(3,IRFs_All,IRFs_LPRigosflex);
IRFs_All = cat(3,IRFs_All,IRFs_proxys);
IRFs_All = cat(3,IRFs_All,IRFs_LPShocks);
IRFs_All = cat(3,IRFs_All,IRFs_LPShockshetero);
IRFs_All = cat(3,IRFs_All,IRFs_LPShocksflex);
IRFs_All = cat(3,IRFs_All,IRFs_LPShocksheteroflex);
IRFs_All = cat(3,IRFs_All,IRFs_rigos);
IRFs_All = cat(3,IRFs_All,IRFs_rigosNaka);

time = (0:horizonLP)';
IRFsupper_LP = max(IRFs_All,[],3);
IRFslower_LP = min(IRFs_All,[],3);

figure('Position',[100 100 1000 600],'PaperPositionMode','Auto','DefaultAxesFontSize',13)
for j=1:nvar % variable
    subplot(2,ceil(nvar/2),j) 

    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper_LP(1,j); IRFslower_LP(1:end,j); flipud([IRFsupper_LP(1:end,j); IRFslower_LP(end,j)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.6);
    set(hh,'edgecolor','none');
    hold on;
    load('../../results/IRFsbench')
    for kk = 1:size(IRFs_All,3)
        p = plot(time, IRFs_All(1:horizonLP+1,j,kk), 'Linewidth', 0.25,'Color',[0.2 0.2 0.2]);
        p.Color(4) = 0.2;
    end
    plot(time, IRFs_base(1:horizonLP+1,j), 'Linewidth', 1.5,'Color','k');
    grid on ;hold off;
    if j==1
        ylim([-20 60])
    elseif j==4 || j==5
        ylim([-4 4])
    elseif j==2 
        ylim([-3 3])
    elseif j==6
        ylim([-.5 2])
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
    print('-dpdf', gcf, '../../results/appendix/figurea6');  
end
    