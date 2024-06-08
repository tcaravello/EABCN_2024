%% Replication files for "The macroeconomic effects of oil supply news"
% Preliminaries: This file creates the monthly data .mat files from the raw data
% stored in rawData.xlsx

% Diego R. Känzig
% LBS, October 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%% Read in the raw data
[dataRaw,labels] = xlsread('rawDataM.xlsx','Data');    % load the data set

header           = labels(1,2:end);
datesStringRaw   = labels(2:end,1);
datesNumRaw      = (str2double(datesStringRaw{1}(1:4))+(str2double(datesStringRaw{1}(end-1:end))-1)*1/12: ...
                    1/12:str2double(datesStringRaw{end}(1:4))+(str2double(datesStringRaw{end}(end-1:end))-1)*1/12)';

Tall = size(dataRaw,1);

% generate variables for series contained in data
for ii = 1:length(header)
    eval([header{ii} '= dataRaw(:,ii);'])
end
    
%% Define some relevant variables 
% base
POIL = WTISPLC; 
OILPROD = EIA1955;
OILSTOCKS = OECDSTOCKSKilianSA; 
WORLDIP = OECD6MNEIP;
IP   = INDPRO;
CPI  = CPIAUCSL;
UNEMP = UNRATE;

% extended
CoreCPI = CPILFESL;
CPIEnergy = CPIENGSL;
CPIND = CUSR0000SAN;		
CPID = CUSR0000SAD;		
CPIS = CUSR0000SAS;
PRATE = FF;
PSTOCK = SPCOMP;
POILEXP = BKEXP12MEXT;
INFLEXP = MICH;
CONS = PCE;
CONSD = PCEDG;
CONSND = PCEND;
CONSS = PCES;
CONSE = DNRGRC1M027SBEA;
CONSDEF = PCEPI;
CONSDDEF = DDURRG3M086SBEA;
CONSNDDEF = DNDGRG3M086SBEA;
CONSSDEF = DSERRG3M086SBEA;
CONSEDEF = DNRGRG3M086SBEA;
NEERB = TWEXBMTH;
NEERM = TWEXMMTH; 
VIX = VXOCLSEXT;
TOT = USTOTPRCF;

% create deterministic variables
const = ones(Tall,1);

%% Transform and save data for benchmark model
data = [log(POIL)*100-log(CPI/100)*100 log(OILPROD)*100 log(OILSTOCKS)*100 log(WORLDIP)*100 log(IP)*100 log(CPI)*100 UNEMP FF log(CoreCPI)*100];  

varNames = {'POIL','OILPROD','OILSTOCKS','WORLDIP','IP','CPI','UNEMP','FF'}; 

% select the sample
smplStart = '1974M01'; 
smplEnd   = '2017M12';

smplStartInd = find(strcmp(datesStringRaw,smplStart));
smplEndInd   = find(strcmp(datesStringRaw,smplEnd));

data = data(smplStartInd:smplEndInd,:);
dataExo = const(smplStartInd:smplEndInd,:);
sampleDates = datesStringRaw(smplStartInd:smplEndInd,:);
sampleDatesNum = datesNumRaw(smplStartInd:smplEndInd,:);

save dataBaseM data dataExo sampleDates sampleDatesNum varNames 

%% Transform and save data for wider effects
% main text
price_data = [log(CoreCPI)*100 log(CPIEnergy)*100 log(CPIND)*100 log(CPID)*100 log(CPIS)*100];
varNames_price = {'Core CPI','CPI energy','CPI non-durables','CPI durables','CPI services'}; 

activity_data = [UNRATE log(CONS./CONSDEF)*100];
varNames_activity = {'Unemployment rate', 'PCE'}; 

expenditure_data = [log(CONSE./CONSEDEF)*100 log(CONSND./CONSNDDEF)*100  ... 
                 log(CONSD./CONSDDEF)*100 log(CONSS./CONSSDEF)*100];
varNames_expenditure = {'PCE energy',' PCE non-durables' ...
                     ,'PCE durables','PCE services'}; 
                 
expectations_data = [log(POILEXP)*100 INFLEXP log(VIX)*100 log(GPR)*100];
varNames_expectations = {'Oil price expectations','Inflation expectations (Michigan)','VXO','Geopolitical risk'}; 

mp_data = [PRATE EBPG log(PSTOCK)*100];
varNames_mp = {'Fed funds rate','Excess bond premium','S\&P 500'}; 

FX_data = [log(NEERB)*100 log(NEERM)*100];
varNames_FX = {'NEER Broad','NEER Major'}; 

bilFX_data = [USCAN USDEN USEUR USINDO USJAP USMEX USNIG USNOR USRUS USSWE USSWI USUK];
bilFX_data = log(bilFX_data)*100;
varNames_bilFX = {'Canada','Denmark','Euro area','Indonesia','Japan','Mexico','Nigeria','Norway','Russia','Sweden','Switzerland','UK'};

data_extended = [price_data activity_data expenditure_data expectations_data mp_data FX_data bilFX_data log(TOT)*100]; 
varNames_paper_extended = [varNames_price varNames_activity varNames_expenditure varNames_expectations varNames_mp varNames_FX varNames_bilFX {'Terms of trade'}];

% select the sample
data_extended = data_extended(smplStartInd:smplEndInd,:);

save dataExtM data_extended sampleDates sampleDatesNum varNames_paper_extended


% appendix

data_extended_stocks = log([OILGSUS ELECTUS MNINGUS AUTOSUS RTAILUS TRLESUS])*100;

varNames_paper_extended_stocks = {'Oil and gas','Electricity','Precious metals','Automobiles and parts','Retail','Travel and leisure'}; 

% select the sample
data_extended_stocks = data_extended_stocks(smplStartInd:smplEndInd,:);

save appendix/dataExtStocksM data_extended_stocks sampleDates sampleDatesNum varNames_paper_extended_stocks


data_extended = [log(USCOCOIMA)*100-log(CPI)*100 GLOBACT];  

varNames_extended = {'RAC','GLOBACT'}; 

% select the sample
data_extended = data_extended(smplStartInd:smplEndInd,:);

% save to disk
save appendix/dataAppendixM data_extended sampleDates sampleDatesNum varNames_extended

%% Transform and save data for quantitative importance
data = [log(POIL)*100-log(CPI)*100 log(OILPROD)*100 log(OILSTOCKS)*100 log(WORLDIP)*100 log(NEERB)*100 log(IP)*100 log(CPI)*100 PRATE log(VIX)*100 log(TOT)*100];  

varNames = {'POIL','OILPROD','OILSTOCKS','WORLDIP','NEER','IP','CPI','FFR','VIX','TOT'}; 

% select the sample
data = data(smplStartInd:smplEndInd,:);

% save to disk
save dataQuantM data dataExo sampleDates sampleDatesNum varNames 

%% Transform and save data for robustness checks

% Granger causality tests
data_level = [log(POIL)*100 log(OILPROD)*100 log(OILSTOCKS)*100 log(WORLDIP)*100 log(IP)*100 log(CPI)*100 PRATE log(PSTOCK)*100 log(NEERB)*100 log(GPRH)*100];  % log(TOT)*100 , seems to be forecastable by it

varNames = {'POIL','OILPROD','OILSTOCKS','WORLDIP','IP','CPI','PRATE','PSTOCK','NEER','GPR'}; 

% number of variables 
nvar = size(data_level,2);    

% set integration order (only relevant if estimated in differences)
diffInd_diff  = [1 1 1 1 1 1 0 1 1 0];

% Transform data if requested
data_diff = nan(size(data_level,1)-1,nvar);
for ii = 1:nvar
    if diffInd_diff(ii)==0
        data_diff(:,ii) = data_level(2:end,ii);
    elseif diffInd_diff(ii)==1
        data_diff(:,ii) = diff(data_level(:,ii)); 
    end
end

% save to disk
save appendix/dataGrangerM data_diff sampleDates sampleDatesNum varNames 


% Stationary Model
data_level = [log(POIL)*100-log(CPI)*100 log(OILPROD)*100 OILSTOCKS log(WORLDIP)*100 log(IP)*100 log(CPI)*100 ];  

varNames = {'POIL','OILPROD','OILSTOCKS','WORLDIP','IP','CPI'}; 

% number of variables 
nvar = size(data_level,2);    

% set integration order (only relevant if estimated in differences)
diffInd_diff     = [0 1 1 1 1 1];

% Transform data if requested
data_diff = nan(size(data_level,1)-1,nvar);
for ii = 1:nvar
    if diffInd_diff(ii)==0
        data_diff(:,ii) = data_level(2:end,ii);
    elseif diffInd_diff(ii)==1
        data_diff(:,ii) = diff(data_level(:,ii)); 
    end
end

data = data_diff;

% save to disk
save appendix/dataStationaryM data dataExo sampleDates sampleDatesNum varNames 


% Model with Brent
data = [log(BRENTEXT)*100-log(CPI)*100 log(OILPROD)*100 log(OILSTOCKS)*100 log(WORLDIP)*100 log(IP)*100 log(CPI)*100];  

varNames = {'POIL','OILPROD','OILSTOCKS','WORLDIP','IP','CPI'}; 

% select the sample
data = data(smplStartInd:smplEndInd,:);

save appendix/dataBrentM data dataExo sampleDates sampleDatesNum varNames 


% BH Model
% replicate their change in stocks as a share of production series
dstock = [nan; OECDSTOCKSKilianSA(2:end)-OECDSTOCKSKilianSA(1:end-1)];
prodm = [nan; OILPROD(1:end-1)/1000*30];
BHSTOCKS = dstock./prodm*100;

data_level = [log(USCOCOIMA)*100-log(CPI/10)*100 log(OILPROD)*100 BHSTOCKS log(WORLDIP)*100];   

varNames = {'POIL','OILPROD','OILSTOCKS','WORLDIP'}; 

% number of variables in VAR
nvar = size(data_level,2);    

% set integration order (only relevant if estimated in differences)
diffInd_diff     = [1 1 0 1];

% Transform data if requested
data_diff = nan(size(data_level,1)-1,nvar);
for ii = 1:nvar
    if diffInd_diff(ii)==0
        data_diff(:,ii) = data_level(2:end,ii);
    elseif diffInd_diff(ii)==1
        data_diff(:,ii) = diff(data_level(:,ii)); 
    end
end

% select the sample
smplStart = '1975M02'; 
smplEnd   = '2016M12'; 

smplStartInd = find(strcmp(datesStringRaw,smplStart));
smplEndInd   = find(strcmp(datesStringRaw,smplEnd));

data = data_diff(smplStartInd-1:smplEndInd-1,:);
dataExo = const(smplStartInd:smplEndInd,:);
sampleDates = datesStringRaw(smplStartInd:smplEndInd,:);
sampleDatesNum = datesNumRaw(smplStartInd:smplEndInd,:);

save appendix/dataBHM data dataExo sampleDates sampleDatesNum varNames 

% Kilian Model
data_level = [log(USCOCOIMA)*100-log(CPI/10)*100 log(OILPROD)*100 OILSTOCKS GLOBACT];  

varNames = {'POIL','OILPROD','OILSTOCKS','GLOBACT'}; 

% number of variables in VAR
nvar = size(data_level,2);    

% set integration order (only relevant if estimated in differences)
diffInd_diff     = [0 1 1 0];

% Transform data if requested
data_diff = nan(size(data_level,1)-1,nvar);
for ii = 1:nvar
    if diffInd_diff(ii)==0
        data_diff(:,ii) = data_level(2:end,ii);
    elseif diffInd_diff(ii)==1
        data_diff(:,ii) = diff(data_level(:,ii)); 
    end
end

% select the sample
smplStart = '1974M02'; 
smplEnd   = '2014M12';

smplStartInd = find(strcmp(datesStringRaw,smplStart));
smplEndInd   = find(strcmp(datesStringRaw,smplEnd));

data = data_diff(smplStartInd-1:smplEndInd-1,:);
dataExo = const(smplStartInd:smplEndInd,:);
sampleDates = datesStringRaw(smplStartInd:smplEndInd,:);
sampleDatesNum = datesNumRaw(smplStartInd:smplEndInd,:);

save appendix/dataKilianM data dataExo sampleDates sampleDatesNum varNames 

