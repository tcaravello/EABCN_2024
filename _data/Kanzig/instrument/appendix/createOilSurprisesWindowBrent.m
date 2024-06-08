%% Create instruments from high frequency data on oil price futures (Brent)
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
[OPECannouncementinfo,OPECannouncements] = xlsread('../../data/appendix/OPECannouncementsEU.xlsx','Daily'); 

statementDates = OPECannouncements(3:end,1);
statementDatesChar = char(statementDates);
statementMonths = strcat(cellstr(statementDatesChar(:,end-3:end)),'M',statementDatesChar(:,4:5));
statementQuarters = strcat(cellstr(statementDatesChar(:,end-3:end)),'Q',num2str(fix((str2num(statementDatesChar(:,4:5))-1)/3)+1));

statementDates_tradingDays = OPECannouncements(3:end,2);

% Information: (1) decision type, (2) meeting type
statementInfo = OPECannouncementinfo(:,[1 2]);

%% Get futures data
% Oil futures (Brent contracts)
% Have to download data from Datastream (Mnemnonics: LLCC.01-LLCC.13, datatypes: PS & VM) 
% and save data for settlement prices in sheet POIL and data on volumes in sheet VOL

% price
[POilRaw,textData] = xlsread('../../data/Oilfutures.xlsx','POIL');
PBrent = log(POilRaw(:,25:end))*100;
futuresDates = textData(3:end,1);
% PBrent should be a 9068x24 matrix containing daily log oil futures prices from
% 30/03/1983 to 29/12/2017

% volume
VOilRaw = xlsread('../../data/Oilfutures.xlsx','VOL');
VBrent = VOilRaw(:,25:end);
% VBrent should be a 9068x24 matrix containing daily volumes for oil futures 
% from 30/03/1983 to 29/12/2017

%% Construct oil price proxies

% set maximum maturity for contract
contractMax = 13;

% windows
wf = 0; % days after the announcement
wb = 1; % days before the announcement

% Compute it for events
% Preallocate
oilProxiesBrentFutures = zeros(length(statementDates),contractMax);
proxyVolumesBrent = zeros(length(statementDates),contractMax,2);

% Compute surprises: futures (do it on selected trading days)
contracts = 1:contractMax;
for jj = 1:length(statementDates)
    dateInd = find(strcmp(futuresDates,statementDates_tradingDays(jj)));
    oilProxiesBrentFutures(jj,contracts) = PBrent(dateInd+wf,contracts)-PBrent(dateInd-wb,contracts);
    proxyVolumesBrent(jj,contracts,1) = VBrent(dateInd+wf,contracts);
    proxyVolumesBrent(jj,contracts,2) = VBrent(dateInd-wb,contracts);
end
oilProxiesBrentFutures(isnan(oilProxiesBrentFutures)) = 0;
oilProxiesBrentFutures(abs(oilProxiesBrentFutures)>15) = 0; % correct for outlier

% principal component(s)
proxySel = [oilProxiesBrentFutures(18:end,2:13)]; 
[~,score,latent,~,explained] = pca(zscore(proxySel));
oilProxiesPC = (score(:,1));
oilProxiesPC = [zeros(17,1); oilProxiesPC./std(oilProxiesPC)*mean(std(proxySel))];

% collect
oilProxiesBrent = [oilProxiesBrentFutures oilProxiesPC];

% aggregate to the monthly frequency
% check for multiplicity
[UniqueSM,~,iSM] = unique(statementMonths);
nmax = max(histc(iSM,1:numel(UniqueSM)));

% Preallocate 
oilProxiesBrentM = zeros(length(sampleDatesProxy),size(oilProxiesBrent,2));
statementInfoM = struct;
statementInfoM.decision = nan(length(sampleDatesProxy),nmax);
statementInfoM.type = nan(length(sampleDatesProxy),nmax);
statementMind = zeros(length(sampleDatesProxy),1);

jj = 1;
for ii = 1:length(sampleDatesProxy)                             % go through all sample dates/months
    if sum(strcmp(statementMonths,sampleDatesProxy(ii)))==1     % if there is only one announcement within this month

        oilProxiesBrentM(ii,:) = oilProxiesBrent(jj,:);
        
        statementInfoM.decision(ii,1) = statementInfo(jj,1);
        statementInfoM.type(ii,1) = statementInfo(jj,2);
        
        statementMind(ii) = 1;
        
        jj = jj + 1;
    elseif sum(strcmp(statementMonths,sampleDatesProxy(ii)))>1     % if there are multiple announcements in the month, sum the shocks
        tempProxy = zeros(1,size(oilProxiesBrent,2));
        for kk = 1:sum(strcmp(statementMonths,sampleDatesProxy(ii)))

            tempProxy = tempProxy + oilProxiesBrent(jj,:);
            
            statementInfoM.decision(ii,kk) = statementInfo(jj,1);
            statementInfoM.type(ii,kk)     = statementInfo(jj,2);
            jj = jj + 1;
        end
        oilProxiesBrentM(ii,:) = tempProxy;
        
        statementMind(ii) = 1;
        
    else
        oilProxiesBrentM(ii,:) = zeros(1,size(oilProxiesBrent,2));
        
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
oilProxiesBrentM = [zeros(size(sampleDatesProxyInit,1),size(oilProxiesBrentM,2)); oilProxiesBrentM];

statementInfoM.decision = [nan(size(sampleDatesProxyInit,1),nmax); statementInfoM.decision];
statementInfoM.type = [nan(size(sampleDatesProxyInit,1),nmax); statementInfoM.type];

statementMind = [zeros(size(sampleDatesProxyInit,1),1); statementMind];

save('OilSurprisesMLogBrent','sampleDatesProxy','oilProxiesBrent','oilProxiesBrentM', ...
     'statementInfoM','statementMind')

 