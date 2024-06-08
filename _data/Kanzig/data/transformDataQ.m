%% Replication files for "The macroeconomic effects of oil supply news"
% Preliminaries: This file creates the quarterly data .mat files from the raw data
% stored in rawData.xlsx

% Diego R. Känzig
% LBS, October 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%% Read in the raw data
[dataRaw,labels] = xlsread('rawDataQ.xlsx','Data');    % load the data set

header           = labels(1,2:end);
datesStringRaw   = labels(2:end,1);
datesNumRaw      = (str2double(datesStringRaw{1}(1:4))+(str2double(datesStringRaw{1}(6))-1)*0.25: ...
                        0.25:str2double(datesStringRaw{end}(1:4))+(str2double(datesStringRaw{end}(6))-1)*0.25)';

Tall = size(dataRaw,1);

% generate variables for series contained in data
for ii = 1:length(header)
    eval([header{ii} '= dataRaw(:,ii);'])
end

%% Define some relevant variables 

% base
POIL = WTISPLC; 
OILPROD = EIA1955;
OILSTOCKS = OECDSTOCKSSAAVG;
WORLDIP = OECD6MNEIP;
RGDP   = GDPC1;
CPI  = CPIAUCSL; 

% extended
SPF = CPI6;
CONS = PCECC96;				
INV = GPDIC1;

% create deterministic variables
const = ones(Tall,1);

%% Transform and save data for benchmark model
data = [log(POIL)*100-log(CPI)*100 log(OILPROD)*100 log(OILSTOCKS)*100 log(WORLDIP)*100 log(RGDP)*100 log(CPI)*100];

varNames = {'POIL','OILPROD','OILSTOCKS','WORLDIP','IP','CPI'}; 

% select the sample
smplStart = '1974Q1'; 
smplEnd   = '2017Q4'; 

smplStartInd = find(strcmp(datesStringRaw,smplStart));
smplEndInd   = find(strcmp(datesStringRaw,smplEnd));

data = data(smplStartInd:smplEndInd,:);
dataExo = const(smplStartInd:smplEndInd,:);
sampleDates = datesStringRaw(smplStartInd:smplEndInd,:);
sampleDatesNum = datesNumRaw(smplStartInd:smplEndInd,:);

save dataBaseQ data dataExo sampleDates sampleDatesNum varNames


%% Transform and save data for wider effects
data_extended = [SPF MICH log(CONS)*100 log(INV)*100 ((USBALGDSB/1000)./GDP)*100];

varNames_paper_extended = {'Inflation expectations (SPF)','Inflation expectations (Michigan survey)','Consumption','Investment','Trade balance'};

% select the sample
data_extended = data_extended(smplStartInd:smplEndInd,:);

save dataExtQ data_extended sampleDates sampleDatesNum varNames_paper_extended
