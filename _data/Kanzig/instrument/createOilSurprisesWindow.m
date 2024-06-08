%% Create instruments from high frequency data on oil price futures
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

% Information: (1) decision type, (2) meeting type
statementInfo = OPECannouncementinfo(:,[1 2]);

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
[~,score,latent,~,explained] = pca(zscore(proxySel));
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

save('OilSurprisesMLog','sampleDatesProxy','oilProxiesWTI','oilProxiesWTIM','statementInfoM','statementMind')


% aggregate to quarterly
sampleDatesProxyQInit = cellstr(strcat(num2str(repelem((startYearInit:endYearInit)',4)),'Q',num2str(repmat((1:4)',endYearInit-startYearInit+1,1))));
sampleDatesProxyQInit = sampleDatesProxyQInit(1:end-3);
sampleDatesProxyQ = [sampleDatesProxyQInit; sampleDatesProxyQ];

oilProxiesWTIQ = nan(size(sampleDatesProxyQ,1),size(oilProxiesWTIM,2));
for ncontract = 1:size(oilProxiesWTIM,2)
    oilProxiesWTIQ(:,ncontract) = sum(reshape(oilProxiesWTIM(:,ncontract),3,[]))';
end

save('OilSurprisesQLog','sampleDatesProxyQ','oilProxiesWTIQ')


%% Construct informationally robust instrument for appendix

% Read in OPEC FCs
actCorrStartInd = find(strcmp(sampleDatesProxy,'2001M01'));
actCorrEndInd = find(strcmp(sampleDatesProxy,'2017M12'));
oilProxiesWTIMShort = oilProxiesWTIM(actCorrStartInd:actCorrEndInd,:);
sampleDatesProxyShort = sampleDatesProxy(actCorrStartInd:actCorrEndInd,:);
statementDatesShort = statementDates(find(strcmp(statementDates,'17/01/2001')):end);
statementDatesShortChar = char(statementDatesShort);
statementMonthsShortTemp = strcat(cellstr(statementDatesShortChar(:,end-3:end)),'M',statementDatesShortChar(:,4:5));
statementMonthsShort = unique(statementMonthsShortTemp); 

[oilDemandRaw,oilDemandLabelsRaw] = xlsread('../data/appendix/OPECreports.xlsx','FCs');
oilDemandLabels = oilDemandLabelsRaw(2,2:end);
oilDemandCurrent = oilDemandRaw(:,1:5);
oilDemandGrowthCurrent = log(oilDemandCurrent(:,2:end))-log(oilDemandCurrent(:,1:end-1));
oilDemandPrevious = oilDemandRaw(:,6:end);
oilDemandGrowthPrevious = log(oilDemandPrevious(:,2:end))-log(oilDemandPrevious(:,1:end-1));
oilDemandRevision = oilDemandCurrent - oilDemandPrevious;
oilDemandGrowthRevision = oilDemandGrowthCurrent - oilDemandGrowthPrevious;

% construct informationally robust instrument by regressing instruments on
% FCs and FC revisions
xcorr = [oilDemandGrowthCurrent oilDemandGrowthRevision]; 

oilProxiesWTIMrefined = zeros(size(oilProxiesWTIM,1),size(oilProxiesWTIM,2));
announcmentsInd = ismember(sampleDatesProxyShort,statementMonthsShort);
for ncontract = 1:size(oilProxiesWTIM,2) 
    proxyShort = oilProxiesWTIMShort(:,ncontract);
    olsProxy = olsest(xcorr(announcmentsInd,:),proxyShort(announcmentsInd,1),true);
    proxyShortRefined = zeros(length(proxyShort),1);
    proxyShortRefined(announcmentsInd,1) = olsProxy.resid;
    oilProxiesWTIMrefined(actCorrStartInd:actCorrEndInd,ncontract) = proxyShortRefined;
end

oilProxiesWTIMcensored = oilProxiesWTIM;
oilProxiesWTIMcensored(1:actCorrStartInd-1,:) = 0;

save('appendix/OilSurprisesMLogRefined','sampleDatesProxy','oilProxiesWTIMcensored','oilProxiesWTIMrefined')

