%% Replication files for "The macroeconomic effects of oil supply news"
% Step 2: This file creates table 1 in the paper

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
p          = 12;    % Lag order   

% proxy
ncontract = 14;     % Principal component of contracts spanning first year of term structure

% switches
saveFigs   = true;  % Save figures to disk


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

%% First stage of external instruments VAR

% baseline proxy
cbench = ncontract;

contracts = [2 3 4 7 10 13 14];
for ncontract = contracts
    
    % load the proxy
    loadProxy;

    % run reduced-form VAR 
    varEst = varxest(data,dataExo,p);

    % identification using the covariance structure between proxy and
    % reduced-form residuals a la Mertens and Ravn (2013) & Stock and Watson (2012)

    % only use proxy sample for identification (potentially a subset of the
    % estimation sample)
    nexo = size(dataExo,2);
    U = varEst.U(smplStartProxyVARInd-p:smplEndProxyVARInd-p,:);    % loose first p observations in estimation  
    Sigma = U'*U/(T-p*nvar-nexo);

    % run first stage
    cnam = strcat('c',num2str(ncontract));
    olsEst.(cnam) = olsest(proxy,U(:,k),true,true);

    fprintf('F-stat: %4.3f, p-value: %4.3f, F-stat (robust): %4.3f, p-value: %4.3f, R^2: %4.3f, R^2 (adj): %4.3f \n',olsEst.(cnam).F,olsEst.(cnam).Fpval,olsEst.(cnam).Frobust,olsEst.(cnam).Frobustpval,olsEst.(cnam).R2,olsEst.(cnam).R2adj)
    uhat = olsEst.(cnam).yhat;

end


%% Table
if saveFigs    
    FID = fopen('../results/table1.tex','w');
    fprintf(FID, strcat('\\begin{tabular}{l',repmat('c',1,7),'}\\toprule\\midrule  \n'));  
    fprintf(FID,' & 1M & 2M & 3M & 6M & 9M & 12M & COMP   \\\\ \\midrule  \n');
    tabVals = [ ];
    for ncontract = contracts
        cnam = strcat('c',num2str(ncontract));
        tabVals = [tabVals olsEst.(cnam).bhat(1)];
    end
    fprintf(FID,' Coefficient & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f   \\\\ \n',tabVals);
    tabVals = [ ];
    for ncontract = contracts
        cnam = strcat('c',num2str(ncontract));
        tabVals = [tabVals olsEst.(cnam).F];
    end
    fprintf(FID,' F-stat & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f   \\\\ \n',tabVals);
    tabVals = [ ];
    for ncontract = contracts
        cnam = strcat('c',num2str(ncontract));
        tabVals = [tabVals olsEst.(cnam).Frobust];
    end
    fprintf(FID,' F-stat (robust)& %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f   \\\\ \n',tabVals);
    tabVals = [ ];
    for ncontract = contracts
        cnam = strcat('c',num2str(ncontract));
        tabVals = [tabVals olsEst.(cnam).R2];
    end
    fprintf(FID,' $R^2$ & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f   \\\\ \n',tabVals*100);
    tabVals = [ ];
    for ncontract = contracts
        cnam = strcat('c',num2str(ncontract));
        tabVals = [tabVals olsEst.(cnam).R2adj];
    end
    fprintf(FID,' $R^2$ (adjusted) & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f & %8.2f   \\\\ \n',tabVals*100);
    tabVals = [ ];
    for ncontract = contracts
        cnam = strcat('c',num2str(ncontract));
        tabVals = [tabVals olsEst.(cnam).n];
    end
    fprintf(FID,' Observations & %8.0f & %8.0f & %8.0f & %8.0f & %8.0f & %8.0f & %8.0f   \\\\ \n',tabVals);
    fprintf(FID, '\\midrule\\bottomrule \n');
    fprintf(FID, '\\end{tabular}\n');
    fclose(FID);
end
