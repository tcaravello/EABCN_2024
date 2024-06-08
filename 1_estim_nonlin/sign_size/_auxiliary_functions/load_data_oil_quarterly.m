Tall = size(DATES,1);
dataFrequency = 'Q'; 

load oil_data_new.mat


% Data manipulations

smplStart = '1983M01'; 
smplEnd   = '2019M12'; 
% select data
data_aux = [log(POIL)*100-log(CPI/100)*100 log(IP)*100 log(CPI)*100 FF UNEMP 100*log(WAGES)];
%data_aux = [log(POIL)*100-log(CPI/100)*100 log(IP)*100 log(CPI)*100 FF];
data = data_aux(find(strcmp(DATES,smplStart)):find(strcmp(DATES,smplEnd)),:);
dates_label_aux = (1974:1/12:2023)';
dates_labels_use = dates_label_aux(find(strcmp(DATES,smplStart)):find(strcmp(DATES,smplEnd)),:);
z = proxy(find(strcmp(DATES,smplStart)):find(strcmp(DATES,smplEnd)));
% set series names
%varNames_level = {'POIL','IP','CPI','FF','UNEMP','RPCE','RWAGES'}; 
%varNames_level = {'POIL','IP','CPI','FF'};
varNames_level = {'POIL','IP','CPI','FF','UNEMP','WAGES'};
varNames_diff  = varNames_level;
%varNames_paper = {'Real oil price','U.S. industrial production','U.S. CPI','Fed Funds Rate','Unemployment','Real Personal Consumption Expenditures','Real Wages'};
varNames_paper = {'Real oil price','U.S. industrial production','U.S. CPI','Fed Funds Rate','Unemployment','Nominal Wages'};

%% Plot data

% close all
% % plot proxy
% figure('DefaultAxesFontSize',13); 
% hold on
% sampleDatesNumSel = dates_labels_use;
% plot(dates_labels_use,z)
% hold on
% plot(dates_labels_use(abs(z)>7), z(abs(z)>7),'ro')
% hold on
% %xlim([dates_labels_use sampleDatesNum(smplEndProxyVARInd)])
% ylim([-15 15])
% xlim([1983 2020])
% text(1987,10.5,'5 Aug 1986','fontsize',11)
% text(1996,-8.6,'14 Nov 2001','fontsize',11)
% text(2009,-11.25,'27 Nov 2014','fontsize',11)
% text(2009,7,'30 Nov 2016','fontsize',11)
% ylabel('Revision in oil price expectations [\%]')
% line(get(gca,'xlim'),[0 0],'Color','k')
% grid on
% box on
% hold off
% subplot(2,1,1)
% plot(dates_labels_use(13:end) ,data(13:end,6)-data(1:end-12,6))
% xlim([1984 2020])
% subplot(2,1,2)
% plot(dates_labels_use(13:end) ,z(13:end))
% xlim([1984 2020])
