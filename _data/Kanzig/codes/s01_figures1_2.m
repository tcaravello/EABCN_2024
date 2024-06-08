%% Replication files for "The macroeconomic effects of oil supply news"
% Step 1: This file creates figures 1 and 2 in the paper

% Diego R. Känzig
% LBS, September 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

% add tools directory
addpath(genpath('auxfiles'))

% initialize random number generator (bad practice but done to be able to
% replicate the results)
rng default

% Set text interpreter to latex
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');


%% Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data frequency
dataFrequency = 'M'; % specify data frequency ('M' or 'Q')

% Instrument sample
smplStartProxy = '1983M04'; 
smplEndProxy   = '2017M12'; 

% select proxy
ncontract = 14;     % Principal component of contracts spanning first year of term structure

% switches
saveFigs = true;    % Save figures to disk


%% Figure 1: The oil supply surprise series

% get sample dates
load('../data/dataBaseM','sampleDates','sampleDatesNum')

% load the proxy
loadProxy;

% plot proxy
figure('DefaultAxesFontSize',13); 
hold on
sampleDatesNumSel = sampleDatesNum(smplStartProxyVARInd:smplEndProxyVARInd);
plot(sampleDatesNumSel,proxy)
plot(sampleDatesNumSel(abs(proxy)>7 & sampleDatesNumSel<2015),proxy(abs(proxy)>7 & sampleDatesNumSel<2015),'ro')
xlim([sampleDatesNum(smplStartProxyVARInd) sampleDatesNum(smplEndProxyVARInd)])
ylim([-15 15])
text(1987,10.5,'5 Aug 1986','fontsize',11)
text(1996,-8.6,'14 Nov 2001','fontsize',11)
text(2009,-11.25,'27 Nov 2014','fontsize',11)
ylabel('Revision in oil price expectations [\%]')
line(get(gca,'xlim'),[0 0],'Color','k')
grid on
box on
if saveFigs
    saveas(gcf,'../results/figure1','epsc2')
end

%% Figure 2: The oil supply surprise series and the control series

% load placebo
load('../instrument/OilSurprisesMLogControl')

proxyControlRaw = [oilPlaceboWTIM(:,ncontract)]; 

proxyControl = proxyControlRaw(smplStartProxyInd:smplEndProxyInd,:);   % we loose the first p values in the VAR


% variance ratio
vratio = var(oilProxiesWTI(:,ncontract))/var(oilPlaceboWTI(:,ncontract));
disp(vratio)


% Brown-Forsythe test
temp = [[oilProxiesWTI(:,ncontract); nan(length(oilPlaceboWTI(:,ncontract))-length(oilProxiesWTI(:,ncontract)),1)] oilPlaceboWTI(:,ncontract)];
[pvalue,Fstat] = vartestn(temp,'TestType','BrownForsythe','Display','off'); 


% plot proxy and placebo
smplStartProxyIndPlot = find(strcmp(sampleDatesProxy,'1984M04'));
smplEndProxyIndPlot   = find(strcmp(sampleDatesProxy,smplEndProxy));
smplStartProxyVARIndPlot = find(strcmp(sampleDates,'1984M04'));
smplEndProxyVARIndPlot   = find(strcmp(sampleDates,smplEndProxy));
   

% compute pdfs
[fproxy,xproxy,bw] = ksdensity(oilProxiesWTI(:,ncontract));
x = xproxy;
pd_kernel = fitdist(oilProxiesWTI(:,ncontract),'kernel','Kernel','epanechnikov'); % 
feproxy = pdf(pd_kernel,x);

[fcontrol,xcontrol,bw] = ksdensity(oilPlaceboWTI(:,ncontract));
pd_kernel = fitdist(oilPlaceboWTI(:,ncontract),'kernel','Kernel','epanechnikov','Bandwidth',0.45); % 
fecontrol = pdf(pd_kernel,xcontrol);


% figure
figure('Position',[100 100 900 350],'PaperPositionMode','Auto','DefaultAxesFontSize',13);
subplot(1,2,1)
hold on
sampleDatesNumSel = sampleDatesNum(smplStartProxyVARIndPlot:smplEndProxyVARIndPlot);
p2=plot(sampleDatesNumSel,proxyControlRaw(smplStartProxyIndPlot:smplEndProxyIndPlot),'color',[0.8500, 0.3250, 0.0980]);
p1=plot(sampleDatesNumSel,proxyRaw(smplStartProxyIndPlot:smplEndProxyIndPlot),'color',[0, 0.4470, 0.7410]);
xlim([sampleDatesNum(smplStartProxyVARIndPlot) sampleDatesNum(smplEndProxyVARIndPlot)])
ylim([-15 15])
xticks([1985:5:2017])
ylabel('\%')
line(get(gca,'xlim'),[0 0],'Color','k')
title('Level')
grid on
box on
str = {strcat('$V_{ann.}/V_{control}: ',sprintf('%2.2f',vratio),'$'), ...
       'H$_0$: $V_{ann.} \leq V_{control}$', ...
       strcat('F-stat: $',sprintf('%2.2f',Fstat.fstat),'$')};
annotation('textbox',[0.14, 0.23, 0.1, 0.1],'String',str,'FitBoxToText','on','EdgeColor','none','Interpreter','latex','fontsize',11);
subplot(1,2,2)
hold on
h1=plot(xproxy,feproxy,'LineWidth',1.5);
h2=plot(xcontrol,fecontrol,'LineWidth',1.5);
grid on 
box on
title('PDF')
xlabel('$x$')
ylabel('$f(x)$')
legend([h1 h2],'Announcement','Control','Autoupdate','off','Interpreter','latex')
uistack(h1,'top')
ylim([0 0.45])
xlim([-15 15])
if saveFigs
    saveas(gcf,'../results/figure2','epsc2')
end

%% Appendix figures

% check autocorrelation
figure('DefaultAxesFontSize',13)
autocorr(proxy)
if saveFigs
    saveas(gcf,'../results/appendix/figurea1','epsc2')
end

% Granger causality
load('../data/appendix/dataGrangerM')
nvar = size(data_diff,2);
varNames_paper = {'Oil price','World oil production','World oil inventories','World industrial production','US industrial production','US CPI','Fed funds rate','S\&P 500','NEER','Geopolitical risk'};

pardl = 12;
x_ardl = data_diff(smplStartProxyVARInd:smplEndProxyVARInd,:);
nx = size(x_ardl,2);
pvals = zeros(nx+2,1);

ardlEst = ardlest(proxy,x_ardl,pardl,pardl,true);

% own lags
R = zeros(pardl,length(ardlEst.bhat));
R(:,1:pardl) = eye(pardl);
q0 = zeros(pardl,1);

% compute Wald statistic
W = (R*ardlEst.bhat-q0)'*inv(R*ardlEst.varbhat*R')*(R*ardlEst.bhat-q0);

Chicrit = chi2inv(0.95,pardl);

pval = 1-chi2cdf(W,pardl);
pvals(1,1) = pval;

% other variables
kk = 1; 
for j = 1:nx

    R = zeros(pardl,length(ardlEst.bhat));
    R(:,pardl*(kk)+1:pardl*(kk)+pardl) = eye(pardl);
    q0 = zeros(pardl,1);

    % compute Wald statistic
    W = (R*ardlEst.bhat-q0)'*inv(R*ardlEst.varbhat*R')*(R*ardlEst.bhat-q0);

    Chicrit = chi2inv(0.95,pardl);

    pval = 1-chi2cdf(W,pardl);
    pvals(1+j,1) = pval;
    kk = kk+1;
end

% all other variables
R = zeros(pardl*nx,length(ardlEst.bhat));
R(:,pardl+1:pardl*nx+pardl) = eye(pardl*nx);
q0 = zeros(pardl*nx,1);

% compute Wald statistic
W = (R*ardlEst.bhat-q0)'*inv(R*ardlEst.varbhat*R')*(R*ardlEst.bhat-q0);

Chicrit = chi2inv(0.95,pardl*nx);

pval = 1-chi2cdf(W,pardl*nx);
pvals(1+nx+1,1) = pval;

% create table
if saveFigs
    fId=fopen('../results/appendix/tablea2.tex','w');
    fprintf(fId,'\\begin{tabular}{lc}\\toprule\\midrule \n');
    fprintf(fId,' Variable & p-value  \\\\  \\midrule \n');
    fprintf(fId,'%s & %4.4f  \\\\ \n','Instrument',pvals(1));
    for ii =1:nvar
        fprintf(fId,'%s & %4.4f  \\\\ \n',varNames_paper{ii},pvals(1+ii));
    end
    fprintf(fId,'%s & %4.4f  \\\\ \n','Joint',pvals(1+nvar+1));
    fprintf(fId,'\\midrule\\bottomrule \n');
    fprintf(fId,'\\end{tabular} \n');
    fclose(fId);
end

