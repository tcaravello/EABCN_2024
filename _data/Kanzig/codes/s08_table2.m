%% Replication files for "The macroeconomic effects of oil supply news"
% Step x: This file creates table 2 in the paper

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
compFEVDs  = true;       % Compute variance decomposition
includeBase = false;     % Include baseline response in plots

verbo = false;           % Control verbosity
saveFigs   = true;       % Save figures to disk


%% Read in data
load("../data/dataQuantM")
% data: transformed endogenous variables 
% dataExo: exogenous variables (e.g. constant, trend)
% sampleDates: sample dates (string format)
% sampleDatesNum: sample dates (numeric format, e.g. 2000 = 2000M1)
% varNames: labels of variables

% number of variables in VAR
nvar = size(data,2);  

% names for paper
varNames_paper = varNames;
varNames_paperVD = {'Oil price','Oil production','Oil inventories','World IP','NEER','IP','CPI','FFR','VXO','TOT'};

% select sample
smplStartInd = find(strcmp(sampleDates,smplStart));
smplEndInd   = find(strcmp(sampleDates,smplEnd));

data = data(smplStartInd:smplEndInd,:);
dataExo = dataExo(smplStartInd:smplEndInd,:);
sampleDates = sampleDates(smplStartInd:smplEndInd,:);
sampleDatesNum = sampleDatesNum(smplStartInd:smplEndInd,:);


%% Large external isntruments VAR for quantitative importance

% load the proxy
loadProxy;

% run proxy VAR
runProxyVAR; 

%% FEVD table
if saveFigs    
    nvarbase = 5;
    FID = fopen('../results/table2.tex','w');
    fprintf(FID, strcat('\\begin{tabular}{l',repmat('c',1,6),'}\\toprule\\midrule  \n'));  
    fprintf(FID,'\\multicolumn{6}{l}{\\textit{Global variables and exchange rates:}} \\\\ \n');  
    nameString = [ ];
    for ii = 1:nvarbase-1
        nameString = [nameString ' %s &'];
    end
    nameString = [nameString ' %s '];
    fprintf(FID,strcat(' &  ',nameString,'  \\\\ \\midrule  \n'),varNames_paperVD{1:nvarbase});
    numString = ' %1.0f & ';
    for ii = 1:nvarbase-1
        numString = [numString '%8.2f & '];
    end
    numString = [numString '%8.2f'];
    numStringBands = [];
    for ii = 1:nvarbase-1
        numStringBands = [numStringBands ' [%1.2f,%8.2f] &'];
    end
    numStringBands = [numStringBands ' [%1.2f,%8.2f]'];
    for kk=[0 12 24 48] %48
        kk2=kk+1;
        if kk<36
            fprintf(FID,strcat(numString,'  \\\\  \n'),kk,FEVDs_proxy(kk2,1:nvarbase));
            fprintf(FID, strcat(' & ',numStringBands,' \\\\  \n'),vec([FEVDslower_proxy(kk2,1:nvarbase)' FEVDsupper_proxy(kk2,1:nvarbase)']'));
        else
            fprintf(FID,strcat(numString,'  \\\\  \n'),kk,FEVDs_proxy(kk2,1:nvarbase));
            fprintf(FID, strcat(' & ',numStringBands,' \\\\[2ex] \\midrule  \n'),vec([FEVDslower_proxy(kk2,1:nvarbase)' FEVDsupper_proxy(kk2,1:nvarbase)']'));
        end
    end
    fprintf(FID,'\\multicolumn{6}{l}{\\textit{U.S. variables:}} \\\\ \n'); 
    nameString = [ ];
    for ii = 1:nvar-nvarbase-1
        nameString = [nameString ' %s &'];
    end
    nameString = [nameString ' %s '];
    fprintf(FID,strcat(' &  ',nameString,'  \\\\ \\midrule  \n'),varNames_paperVD{nvarbase+1:end});
    numString = ' %1.0f & ';
    for ii = 1:nvar-nvarbase-1
        numString = [numString '%8.2f & '];
    end
    numString = [numString '%8.2f'];
    numStringBands = [];
    for ii = 1:nvar-nvarbase-1
        numStringBands = [numStringBands ' [%1.2f,%8.2f] &'];
    end
    numStringBands = [numStringBands ' [%1.2f,%8.2f]'];
    for kk=[0 12 24 48] 
        kk2=kk+1;
        fprintf(FID,strcat(numString,'  \\\\  \n'),kk,FEVDs_proxy(kk2,nvarbase+1:end));
        fprintf(FID, strcat(' & ',numStringBands,' \\\\  \n'),vec([FEVDslower_proxy(kk2,nvarbase+1:end)' FEVDsupper_proxy(kk2,nvarbase+1:end)']'));

    end
    fprintf(FID, '\\midrule\\bottomrule \n');
    fprintf(FID, '\\end{tabular}\n');
    fclose(FID);
end

