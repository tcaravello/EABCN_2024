%% Create instruments and control series from high frequency data on oil price futures
% Diego R. Kaenzig
% LBS, September 2020

clear all
close all
clc

addpath(genpath('../codes/auxfiles'))

%% Define sample
startYear = 1983;
endYear   = 2017;
sampleDatesProxy = cellstr(strcat(num2str(repelem((startYear:endYear)',12)),'M',num2str(repmat((1:12)',endYear-startYear+1,1))));
sampleDatesProxy = strrep(sampleDatesProxy,' ','0');
sampleDatesNum   = (str2double(sampleDatesProxy{1}(1:4))+(str2double(sampleDatesProxy{1}(end-1:end))-1)*1/12: ...
                    1/12:str2double(sampleDatesProxy{end}(1:4))+(str2double(sampleDatesProxy{end}(end-1:end))-1)*1/12)';
                
sampleDatesProxyQ = cellstr(strcat(num2str(repelem((startYear:endYear)',4)),'Q',num2str(repmat((1:4)',endYear-startYear+1,1))));
                
sampleDatesProxy = sampleDatesProxy(4:end); % futures only start end of march/ beginning of april
sampleDatesNum = sampleDatesNum(4:end);

sampleDatesProxyQ = sampleDatesProxyQ(2:end);

%% Get all announcement dates 
% OPEC statements
[OPECannouncementinfo,OPECannouncements] = xlsread('../data/OPECannouncements.xlsx','Daily');

statementDates = OPECannouncements(3:end,1);
statementDatesChar = char(statementDates);
statementMonths = strcat(cellstr(statementDatesChar(:,end-3:end)),'M',statementDatesChar(:,4:5));
statementQuarters = strcat(cellstr(statementDatesChar(:,end-3:end)),'Q',num2str(fix((str2num(statementDatesChar(:,4:5))-1)/3)+1));

statementDates_tradingDays = OPECannouncements(3:end,2);

statementInfo = OPECannouncementinfo(:,[1 2]);

% placebo/control series
[~,AllDates] = xlsread('../data/OPECplacebos.xlsx','Control'); 
placeboDates = AllDates(2:end,3);
placeboMonths = AllDates(2:end,1);

% only keep relevant dates
emptyCells = cellfun('isempty',placeboDates);
placeboDates(emptyCells) = [];
placeboMonths(emptyCells) = [];

%% Get futures data
% Oil futures (WTI contracts)
% Have to download data from Datastream (Mnemnonics: NCLC.01-NCLC.13, datatypes: PS & VM) 
% and save data for settlement prices in sheet POIL and data on volumes in sheet VOL

% price
[POilRaw,textData] = xlsread('../data/Oilfutures.xlsx','POIL');
PWTI = log(POilRaw(:,1:24))*100;
futuresDates = textData(3:end,1);
% PWTI should be a 9068x24 matrix containing daily log oil futures prices from
% 30/03/1983 to 29/12/2017

% volume
VOilRaw = xlsread('../data/Oilfutures.xlsx','VOL');
VWTI = VOilRaw(:,1:24);
% VWTI should be a 9068x24 matrix containing daily volumes for oil futures 
% from 30/03/1983 to 29/12/2017

%% Construct oil price proxies

% treatment:
% set maximum maturity for contract
contractMax = 13;

% windows
wf = 0; % days after the announcement
wb = 1; % days before the announcement

% Compute it for events
% Preallocate
oilProxiesWTIFutures = zeros(length(statementDates),contractMax);
proxyVolumesWTI = zeros(length(statementDates),contractMax,2);

% Compute surprises: futures (do it on selected trading days)
contracts = 1:contractMax;
for jj = 1:length(statementDates)
    dateInd = find(strcmp(futuresDates,statementDates_tradingDays(jj)));
    oilProxiesWTIFutures(jj,contracts) = PWTI(dateInd+wf,contracts)-PWTI(dateInd-wb,contracts);
    proxyVolumesWTI(jj,contracts,1) = VWTI(dateInd+wf,contracts);
    proxyVolumesWTI(jj,contracts,2) = VWTI(dateInd-wb,contracts);
end
oilProxiesWTIFutures(isnan(oilProxiesWTIFutures)) = 0;

% principal component(s)
proxySel = [oilProxiesWTIFutures(:,2:13)]; 
[coefs,score,latent,~,explained] = pca(zscore(proxySel));
oilProxiesPC = (score(:,1));
oilProxiesPC =oilProxiesPC./std(oilProxiesPC)*mean(std(proxySel));

% collect
oilProxiesWTI = [oilProxiesWTIFutures oilProxiesPC];

% aggregate to the monthly frequency
% check for multiplicity
[UniqueSM,~,iSM] = unique(statementMonths);
nmax = max(histc(iSM,1:numel(UniqueSM)));

% Preallocate 
oilProxiesWTIM = zeros(length(sampleDatesProxy),size(oilProxiesWTI,2));
statementInfoM = struct;
statementInfoM.decision = nan(length(sampleDatesProxy),nmax);
statementInfoM.type = nan(length(sampleDatesProxy),nmax);
statementMind = zeros(length(sampleDatesProxy),1);

jj = 1;
for ii = 1:length(sampleDatesProxy)                             % go through all sample dates/months
    if sum(strcmp(statementMonths,sampleDatesProxy(ii)))==1     % if there is only one announcement within this month

        oilProxiesWTIM(ii,:) = oilProxiesWTI(jj,:);
        
        statementInfoM.decision(ii,1) = statementInfo(jj,1);
        statementInfoM.type(ii,1) = statementInfo(jj,2);
        
        statementMind(ii) = 1;
        
        jj = jj + 1;
    elseif sum(strcmp(statementMonths,sampleDatesProxy(ii)))>1     % if there are multiple announcements in the month, sum the shocks
        tempProxy = zeros(1,size(oilProxiesWTI,2));
        for kk = 1:sum(strcmp(statementMonths,sampleDatesProxy(ii)))

            tempProxy = tempProxy + oilProxiesWTI(jj,:);
            
            statementInfoM.decision(ii,kk) = statementInfo(jj,1);
            statementInfoM.type(ii,kk)     = statementInfo(jj,2);
            jj = jj + 1;
        end
        oilProxiesWTIM(ii,:) = tempProxy;
        
        statementMind(ii) = 1;
        
    else
        oilProxiesWTIM(ii,:) = zeros(1,size(oilProxiesWTI,2));
        
        statementInfoM.decision(ii,1) = nan;
        statementInfoM.type(ii,1) = nan;
    end
end


% control:
usePCcoefs = false;

% Compute it for events
% Preallocate
oilPlaceboWTIFutures = zeros(length(placeboDates),contractMax);
oilPlaceboWTISpot = zeros(length(placeboDates),1);
placeboVolumesWTI = zeros(length(placeboDates),contractMax,2);

% Compute surprises: futures (do it on selected trading days)
contracts = 1:contractMax;
for jj = 1:length(placeboDates)
    dateInd = find(strcmp(futuresDates,placeboDates(jj)));
    oilPlaceboWTIFutures(jj,contracts) = PWTI(dateInd+wf,contracts)-PWTI(dateInd-wb,contracts);
    placeboVolumesWTI(jj,contracts,1) = VWTI(dateInd+wf,contracts);
    placeboVolumesWTI(jj,contracts,2) = VWTI(dateInd-wb,contracts);
end
oilPlaceboWTIFutures(isnan(oilPlaceboWTIFutures)) = 0;

% principal component(s)
if usePCcoefs
    placeboSel = [oilPlaceboWTIFutures(:,2:13)];
    oilPlaceboPC = (zscore(placeboSel)*coefs(:,1));
    oilPlaceboPC =oilPlaceboPC./std(oilPlaceboPC)*mean(std(placeboSel));
else
    placeboSel = [oilPlaceboWTIFutures(:,2:13)];
    [~,score,latent,~,explained] = pca(zscore(placeboSel));
    oilPlaceboPC = (score(:,1));
    oilPlaceboPC =oilPlaceboPC./std(oilPlaceboPC)*mean(std(placeboSel));
end

% collect
oilPlaceboWTI = [oilPlaceboWTIFutures oilPlaceboPC];

% aggregate to the monthly frequency
% Preallocate 
oilPlaceboWTIM = zeros(length(sampleDatesProxy),size(oilPlaceboWTI,2));
placeboMind = zeros(length(sampleDatesProxy),1);

jj = 1;
for ii = 1:length(sampleDatesProxy)                             % go through all sample dates/months
    if sum(strcmp(placeboMonths,sampleDatesProxy(ii)))==1     % if there is only one announcement within this month

        oilPlaceboWTIM(ii,:) = oilPlaceboWTI(jj,:);
        
        placeboMind(ii) = 1;
        
        jj = jj + 1;
    elseif sum(strcmp(placeboMonths,sampleDatesProxy(ii)))>1     % if there are multiple announcements in the month, sum the shocks
        tempPlacebo = zeros(1,size(oilPlaceboWTI,2));
        for kk = 1:sum(strcmp(placeboMonths,sampleDatesProxy(ii)))

            tempPlacebo = tempPlacebo + oilPlaceboWTI(jj,:);
            
            jj = jj + 1;
        end
        oilPlaceboWTIM(ii,:) = tempPlacebo;
        
        placeboMind(ii) = 1;
        
    else
        oilPlaceboWTIM(ii,:) = zeros(1,size(oilPlaceboWTI,2));
        
    end
end

% extend the proxy to span the entire estimation sample (back to 1974M1)
startYearInit = 1974;
endYearInit   = 1983;
sampleDatesProxyInit = cellstr(strcat(num2str(repelem((startYearInit:endYearInit)',12)),'M',num2str(repmat((1:12)',endYearInit-startYearInit+1,1))));
sampleDatesProxyInit = strrep(sampleDatesProxyInit,' ','0');
sampleDatesProxyInit = sampleDatesProxyInit(1:end-9); % get rid of last 3 quarters in 1983
sampleDatesNumInit   = (str2double(sampleDatesProxyInit{1}(1:4))+(str2double(sampleDatesProxyInit{1}(end-1:end))-1)*1/12: ...
                    1/12:str2double(sampleDatesProxyInit{end}(1:4))+(str2double(sampleDatesProxyInit{end}(end-1:end))-1)*1/12)';

sampleDatesProxy = [sampleDatesProxyInit; sampleDatesProxy];
sampleDatesNumExt = [sampleDatesNumInit; sampleDatesNum];
oilProxiesWTIM = [zeros(size(sampleDatesProxyInit,1),size(oilProxiesWTIM,2)); oilProxiesWTIM];

statementInfoM.decision = [nan(size(sampleDatesProxyInit,1),nmax); statementInfoM.decision];
statementInfoM.type = [nan(size(sampleDatesProxyInit,1),nmax); statementInfoM.type];

statementMind = [zeros(size(sampleDatesProxyInit,1),1); statementMind];

oilPlaceboWTIM = [zeros(size(sampleDatesProxyInit,1),size(oilPlaceboWTIM,2)); oilPlaceboWTIM];

placeboMind = [zeros(size(sampleDatesProxyInit,1),1); placeboMind];

save('OilSurprisesMLogControl','sampleDatesProxy','oilProxiesWTI','oilProxiesWTIM','statementMind','oilPlaceboWTI','oilPlaceboWTIM','placeboMind')



