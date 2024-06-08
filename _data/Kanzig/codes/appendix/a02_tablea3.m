%% Replication files for "The macroeconomic effects of oil supply news"
% This file creates Table A.3. in the appendix
% The file calls 'shockseries.xlsx',
% which I cannot share publicly because it was shared with me by another researcher

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
set(0,'defaulttextinterpreter','latex')

saveFigs = true;

%% Monthly
% load shock
load('../../instrument/OilSurprisesMLog')
ncontract = 14;

% load data
[~,~,infoRaw] = xlsread('../../data/appendix/shockseries.xlsx','M');
header = infoRaw(2,2:end);
datesStringRaw = infoRaw(3:end,1);
datesNumRaw      = (str2double(datesStringRaw{1}(1:4))+(str2double(datesStringRaw{1}(end-1:end))-1)*1/12: ...
                    1/12:str2double(datesStringRaw{end}(1:4))+(str2double(datesStringRaw{end}(end-1:end))-1)*1/12)';
dataRaw = infoRaw(3:end,2:end);
dataRaw(strcmp(dataRaw,'ActiveX VT_ERROR: ')) = {nan};
dataRaw = cell2mat(dataRaw);

% select sample
smplStartProxy = '1974M01'; 
smplEndProxy = '2017M12';

smplStartInd = find(strcmp(datesStringRaw,smplStartProxy));
smplEndInd = find(strcmp(datesStringRaw,smplEndProxy));
data = dataRaw(smplStartInd:smplEndInd,:);

smplStartProxyInd = find(strcmp(sampleDatesProxy,smplStartProxy));
smplEndProxyInd   = find(strcmp(sampleDatesProxy,smplEndProxy));
proxy = oilProxiesWTIM(smplStartProxyInd:smplEndProxyInd,ncontract);
sampleDatesSel = sampleDatesProxy(smplStartProxyInd:smplEndProxyInd);

disp('Monthly')
for ii = [1 3:length(header)] % 1:length(header); does not work for second column as series cannot be shared
    [RHO,PVAL] =corr(proxy(~isnan(data(:,ii))),data(~isnan(data(:,ii)),ii));
    fprintf('Coeff: %4.2f, p-val: %4.2f, nobs: %4.0f, sample: %s to %s \n',RHO,PVAL,length(data(~isnan(data(:,ii)),ii)),sampleDatesSel{find(~isnan(data(:,ii)), 1)},sampleDatesSel{find(~isnan(data(:,ii)), 1, 'last')})
end
    
%% Quarterly
load('../../instrument/OilSurprisesQLog')

% load data
[~,~,infoRawQ] = xlsread('../../data/appendix/shockseries.xlsx','Q');
headerQ = infoRawQ(2,2:end);
datesStringRawQ = infoRawQ(3:end,1);
dataRawQ = infoRawQ(3:end,2:end);
dataRawQ(strcmp(dataRawQ,'ActiveX VT_ERROR: ')) = {nan};
dataRawQ = cell2mat(dataRawQ);

% select sample
smplStartProxyQ = mtoqdate(smplStartProxy); 
smplEndProxyQ = mtoqdate(smplEndProxy);

smplStartIndQ = find(strcmp(datesStringRawQ,smplStartProxyQ));
smplEndIndQ = find(strcmp(datesStringRawQ,smplEndProxyQ));
dataQ = dataRawQ(smplStartIndQ:smplEndIndQ,:);

smplStartProxyIndQ = find(strcmp(sampleDatesProxyQ,smplStartProxyQ));
smplEndProxyIndQ   = find(strcmp(sampleDatesProxyQ,smplEndProxyQ));
proxyQ = oilProxiesWTIQ(smplStartProxyIndQ:smplEndProxyIndQ,ncontract);
sampleDatesSelQ = sampleDatesProxyQ(smplStartProxyIndQ:smplEndProxyIndQ);

disp(' ')
disp('Quarterly')
for ii = 1:length(headerQ)
    [RHO,PVAL] =corr(proxyQ(~isnan(dataQ(:,ii))),dataQ(~isnan(dataQ(:,ii)),ii));
    fprintf('Coeff: %4.2f, p-val: %4.2f, nobs: %4.0f, sample: %s to %s \n',RHO,PVAL,length(dataQ(~isnan(dataQ(:,ii)),ii)),sampleDatesSelQ{find(~isnan(dataQ(:,ii)), 1)},sampleDatesSelQ{find(~isnan(dataQ(:,ii)), 1, 'last')})
end


%% Create table

if saveFigs    
    FID = fopen('../../results/appendix/tablea3.tex','w');
    fprintf(FID, strcat('\\begin{tabular}{llcccl}\\toprule\\midrule  \n'));  
    fprintf(FID,'Shock & \\multicolumn{1}{c}{Source} & $\\rho$ & p-value & $n$ & \\multicolumn{1}{c}{Sample}    \\\\ \\midrule  \n');
    fprintf(FID,' \\multicolumn{6}{l}{\\textit{Panel A: Oil shocks}} \\\\  \n');
    % Oil
    ii = 1;    
    [RHO,PVAL] =corr(proxy(~isnan(data(:,ii))),data(~isnan(data(:,ii)),ii));
    fprintf(FID,' Oil price & \\cite{hamilton2003oil} & %4.2f & %4.2f &  %4.0f & %s-%s \\\\  \n',RHO,PVAL,length(data(~isnan(data(:,ii)),ii)),sampleDatesSel{find(~isnan(data(:,ii)), 1)},sampleDatesSel{find(~isnan(data(:,ii)), 1, 'last')});  
    ii = 2;    
    % [RHO,PVAL] =corr(proxy(~isnan(data(:,ii))),data(~isnan(data(:,ii)),ii));
    % fprintf(FID,' Oil supply & \\cite{kilian2008exogenous} & %4.2f & %4.2f &  %4.0f & %s-%s  \\\\  \n',RHO,PVAL,length(data(~isnan(data(:,ii)),ii)),sampleDatesSel{find(~isnan(data(:,ii)), 1)},sampleDatesSel{find(~isnan(data(:,ii)), 1, 'last')});
    fprintf(FID,' Oil supply & \\cite{kilian2008exogenous} &  &  &   &   \\\\  \n');
    ii = 3;    
    [RHO,PVAL] =corr(proxy(~isnan(data(:,ii))),data(~isnan(data(:,ii)),ii));
    fprintf(FID,'  & \\cite{caldara2018oil} & %4.2f & %4.2f &  %4.0f & %s-%s  \\\\  \n',RHO,PVAL,length(data(~isnan(data(:,ii)),ii)),sampleDatesSel{find(~isnan(data(:,ii)), 1)},sampleDatesSel{find(~isnan(data(:,ii)), 1, 'last')});
    ii = 4;   
    [RHO,PVAL] =corr(proxy(~isnan(data(:,ii))),data(~isnan(data(:,ii)),ii));
    fprintf(FID,'  & \\cite{baumeister2017structural} & %4.2f & %4.2f &  %4.0f & %s-%s  \\\\  \n',RHO,PVAL,length(data(~isnan(data(:,ii)),ii)),sampleDatesSel{find(~isnan(data(:,ii)), 1)},sampleDatesSel{find(~isnan(data(:,ii)), 1, 'last')});
    ii = 5;    
    [RHO,PVAL] =corr(proxy(~isnan(data(:,ii))),data(~isnan(data(:,ii)),ii));
    fprintf(FID,'  & \\cite{kilian2009not} & %4.2f & %4.2f &  %4.0f & %s-%s  \\\\  \n',RHO,PVAL,length(data(~isnan(data(:,ii)),ii)),sampleDatesSel{find(~isnan(data(:,ii)), 1)},sampleDatesSel{find(~isnan(data(:,ii)), 1, 'last')});
    ii = 6;    
    [RHO,PVAL] =corr(proxy(~isnan(data(:,ii))),data(~isnan(data(:,ii)),ii));
    fprintf(FID,' Global demand & \\cite{kilian2009not} & %4.2f & %4.2f &  %4.0f & %s-%s  \\\\  \n',RHO,PVAL,length(data(~isnan(data(:,ii)),ii)),sampleDatesSel{find(~isnan(data(:,ii)), 1)},sampleDatesSel{find(~isnan(data(:,ii)), 1, 'last')});
    ii = 7;    
    [RHO,PVAL] =corr(proxy(~isnan(data(:,ii))),data(~isnan(data(:,ii)),ii));
    fprintf(FID,'Oil-specific demand & \\cite{kilian2009not} & %4.2f & %4.2f &  %4.0f & %s-%s  \\\\  \n',RHO,PVAL,length(data(~isnan(data(:,ii)),ii)),sampleDatesSel{find(~isnan(data(:,ii)), 1)},sampleDatesSel{find(~isnan(data(:,ii)), 1, 'last')});

    fprintf(FID,' \\multicolumn{6}{l}{} \\\\  \n'); 
    fprintf(FID,' \\multicolumn{6}{l}{\\textit{Panel B: Other shocks}} \\\\  \n'); 
    % Productivity
    ii = 8;    
    [RHO,PVAL] =corr(proxyQ(~isnan(dataQ(:,ii))),dataQ(~isnan(dataQ(:,ii)),ii));
    fprintf(FID,'Productivity & \\cite{basu2006technology}  & %4.2f & %4.2f &  %4.0f & %s-%s  \\\\  \n',RHO,PVAL,length(dataQ(~isnan(dataQ(:,ii)),ii)),sampleDatesSelQ{find(~isnan(dataQ(:,ii)), 1)},sampleDatesSelQ{find(~isnan(dataQ(:,ii)), 1, 'last')});
    ii = 7;    
    [RHO,PVAL] =corr(proxyQ(~isnan(dataQ(:,ii))),dataQ(~isnan(dataQ(:,ii)),ii));
    fprintf(FID,' & \\cite{smets2007shocks} & %4.2f & %4.2f &  %4.0f & %s-%s  \\\\  \n',RHO,PVAL,length(dataQ(~isnan(dataQ(:,ii)),ii)),sampleDatesSelQ{find(~isnan(dataQ(:,ii)), 1)},sampleDatesSelQ{find(~isnan(dataQ(:,ii)), 1, 'last')});
     % News
    ii = 10;    
    [RHO,PVAL] =corr(proxyQ(~isnan(dataQ(:,ii))),dataQ(~isnan(dataQ(:,ii)),ii));
    fprintf(FID,'News & \\cite{barsky2011news}  & %4.2f & %4.2f &  %4.0f & %s-%s  \\\\  \n',RHO,PVAL,length(dataQ(~isnan(dataQ(:,ii)),ii)),sampleDatesSelQ{find(~isnan(dataQ(:,ii)), 1)},sampleDatesSelQ{find(~isnan(dataQ(:,ii)), 1, 'last')});
    ii = 11;    
    [RHO,PVAL] =corr(proxyQ(~isnan(dataQ(:,ii))),dataQ(~isnan(dataQ(:,ii)),ii));
    fprintf(FID,' & \\cite{kurmann2013news} & %4.2f & %4.2f &  %4.0f & %s-%s  \\\\  \n',RHO,PVAL,length(dataQ(~isnan(dataQ(:,ii)),ii)),sampleDatesSelQ{find(~isnan(dataQ(:,ii)), 1)},sampleDatesSelQ{find(~isnan(dataQ(:,ii)), 1, 'last')});
    ii = 12;    
    [RHO,PVAL] =corr(proxyQ(~isnan(dataQ(:,ii))),dataQ(~isnan(dataQ(:,ii)),ii));
    fprintf(FID,'& \\cite{beaudry2014news}  & %4.2f & %4.2f &  %4.0f & %s-%s  \\\\  \n',RHO,PVAL,length(dataQ(~isnan(dataQ(:,ii)),ii)),sampleDatesSelQ{find(~isnan(dataQ(:,ii)), 1)},sampleDatesSelQ{find(~isnan(dataQ(:,ii)), 1, 'last')});
    % MP
    ii = 8;    
    [RHO,PVAL] =corr(proxy(~isnan(data(:,ii))),data(~isnan(data(:,ii)),ii));
    fprintf(FID,'Monetary policy & \\cite{gertler2015monetary} & %4.2f & %4.2f &  %4.0f & %s-%s  \\\\  \n',RHO,PVAL,length(data(~isnan(data(:,ii)),ii)),sampleDatesSel{find(~isnan(data(:,ii)), 1)},sampleDatesSel{find(~isnan(data(:,ii)), 1, 'last')});
    ii = 9;    
    [RHO,PVAL] =corr(proxy(~isnan(data(:,ii))),data(~isnan(data(:,ii)),ii));
    fprintf(FID,' & \\cite{romer2004new} & %4.2f & %4.2f &  %4.0f & %s-%s  \\\\  \n',RHO,PVAL,length(data(~isnan(data(:,ii)),ii)),sampleDatesSel{find(~isnan(data(:,ii)), 1)},sampleDatesSel{find(~isnan(data(:,ii)), 1, 'last')});
    ii = 1;    
    [RHO,PVAL] =corr(proxyQ(~isnan(dataQ(:,ii))),dataQ(~isnan(dataQ(:,ii)),ii));
    fprintf(FID,' & \\cite{smets2007shocks} & %4.2f & %4.2f &  %4.0f & %s-%s  \\\\  \n',RHO,PVAL,length(dataQ(~isnan(dataQ(:,ii)),ii)),sampleDatesSelQ{find(~isnan(dataQ(:,ii)), 1)},sampleDatesSelQ{find(~isnan(dataQ(:,ii)), 1, 'last')});
    % Uncertainty
    ii = 10;    
    [RHO,PVAL] =corr(proxy(~isnan(data(:,ii))),data(~isnan(data(:,ii)),ii));
    fprintf(FID,' Uncertainty & \\cite{bloom2009impact} & %4.2f & %4.2f &  %4.0f & %s-%s \\\\  \n',RHO,PVAL,length(data(~isnan(data(:,ii)),ii)),sampleDatesSel{find(~isnan(data(:,ii)), 1)},sampleDatesSel{find(~isnan(data(:,ii)), 1, 'last')});  
    ii = 11;    
    [RHO,PVAL] =corr(proxy(~isnan(data(:,ii))),data(~isnan(data(:,ii)),ii));
    fprintf(FID,'  & \\cite{baker2016measuring} & %4.2f & %4.2f &  %4.0f & %s-%s  \\\\  \n',RHO,PVAL,length(data(~isnan(data(:,ii)),ii)),sampleDatesSel{find(~isnan(data(:,ii)), 1)},sampleDatesSel{find(~isnan(data(:,ii)), 1, 'last')});
    % Financial
    ii = 12;    
    [RHO,PVAL] =corr(proxy(~isnan(data(:,ii))),data(~isnan(data(:,ii)),ii));
    fprintf(FID,' Financial  & \\cite{gilchrist2012credit} & %4.2f & %4.2f &  %4.0f & %s-%s  \\\\  \n',RHO,PVAL,length(data(~isnan(data(:,ii)),ii)),sampleDatesSel{find(~isnan(data(:,ii)), 1)},sampleDatesSel{find(~isnan(data(:,ii)), 1, 'last')});
    ii = 9;    
    [RHO,PVAL] =corr(proxyQ(~isnan(dataQ(:,ii))),dataQ(~isnan(dataQ(:,ii)),ii));
    fprintf(FID,' & \\cite{bassett2014changes}  & %4.2f & %4.2f &  %4.0f & %s-%s  \\\\  \n',RHO,PVAL,length(dataQ(~isnan(dataQ(:,ii)),ii)),sampleDatesSelQ{find(~isnan(dataQ(:,ii)), 1)},sampleDatesSelQ{find(~isnan(dataQ(:,ii)), 1, 'last')});
    % Fiscal policy
    ii = 2;    
    [RHO,PVAL] =corr(proxyQ(~isnan(dataQ(:,ii))),dataQ(~isnan(dataQ(:,ii)),ii));
    fprintf(FID,'Fiscal policy & \\cite{romer2010macroeconomic}  & %4.2f & %4.2f &  %4.0f & %s-%s  \\\\  \n',RHO,PVAL,length(dataQ(~isnan(dataQ(:,ii)),ii)),sampleDatesSelQ{find(~isnan(dataQ(:,ii)), 1)},sampleDatesSelQ{find(~isnan(dataQ(:,ii)), 1, 'last')});
    ii = 3;    
    [RHO,PVAL] =corr(proxyQ(~isnan(dataQ(:,ii))),dataQ(~isnan(dataQ(:,ii)),ii));
    fprintf(FID,' & \\cite{ramey2011identifying} & %4.2f & %4.2f &  %4.0f & %s-%s  \\\\  \n',RHO,PVAL,length(dataQ(~isnan(dataQ(:,ii)),ii)),sampleDatesSelQ{find(~isnan(dataQ(:,ii)), 1)},sampleDatesSelQ{find(~isnan(dataQ(:,ii)), 1, 'last')});
    ii = 4;    
    [RHO,PVAL] =corr(proxyQ(~isnan(dataQ(:,ii))),dataQ(~isnan(dataQ(:,ii)),ii));
    fprintf(FID,'& \\cite{fisher2010using} & %4.2f & %4.2f &  %4.0f & %s-%s  \\\\  \n',RHO,PVAL,length(dataQ(~isnan(dataQ(:,ii)),ii)),sampleDatesSelQ{find(~isnan(dataQ(:,ii)), 1)},sampleDatesSelQ{find(~isnan(dataQ(:,ii)), 1, 'last')});

    fprintf(FID, '\\midrule\\bottomrule \n');
    fprintf(FID, '\\end{tabular}\n');
    fclose(FID);
end
