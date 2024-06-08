%% Replication files for "The macroeconomic effects of oil supply news"
% This file creates Table A.1. in the appendix
% The file calls 'FED announcements.xlsx' and 'macronewskey_updated.xls',
% which I cannot share publicly because they were shared with me by another researcher

% Diego R. Känzig
% LBS, September 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define sample
startYear = 1983;
endYear   = 2017;
sampleDatesProxy = cellstr(strcat(num2str(repelem((startYear:endYear)',12)),'M',num2str(repmat((1:12)',endYear-startYear+1,1))));
sampleDatesProxy = strrep(sampleDatesProxy,' ','0');
sampleDatesNum   = (str2double(sampleDatesProxy{1}(1:4))+(str2double(sampleDatesProxy{1}(end-1:end))-1)*1/12: ...
                    1/12:str2double(sampleDatesProxy{end}(1:4))+(str2double(sampleDatesProxy{end}(end-1:end))-1)*1/12)';
                            
sampleDatesProxy = sampleDatesProxy(4:end); % futures only start end of march/ beginning of april
sampleDatesNum = sampleDatesNum(4:end);


% load announcements
[~,OPECannouncements] = xlsread('../../data/OPECannouncements.xlsx','Daily'); 
 
statementDatesOPEC = OPECannouncements(3:end,1); % 1
statementDatesChar = char(statementDatesOPEC);
statementMonths = strcat(cellstr(statementDatesChar(:,end-3:end)),'M',statementDatesChar(:,4:5));

%% FOMC

[~,FEDannouncements] = xlsread('../../data/appendix/FED announcements.xlsx','Final');

statementDatesFed = FEDannouncements(2:end,1);

flagFED = zeros(length(statementDatesOPEC),1);

for ii=1:length(statementDatesOPEC)
    flagFED(ii) = sum(strcmp(statementDatesFed,statementDatesOPEC(ii)));
end

sum(flagFED)

%% Other news
% GDP Advanced
[~,temp] = xlsread('../../data/appendix/macronewskey_updated.xls','GDP Advanced');

datesGDP = temp(2:end,1);

flagGDP = zeros(length(statementDatesOPEC),1);

for ii=1:length(statementDatesOPEC)
    flagGDP(ii) = sum(strcmp(datesGDP,statementDatesOPEC(ii)));
end

sum(flagGDP)

% Unemployment Rate
[~,temp] = xlsread('../../data/appendix/macronewskey_updated.xls','Unemployment Rate');

datesUR = temp(2:end,1);

flagUR = zeros(length(statementDatesOPEC),1);

for ii=1:length(statementDatesOPEC)
    flagUR(ii) = sum(strcmp(datesUR,statementDatesOPEC(ii)));
end

sum(flagUR)

% Nonfarm payroll
[~,temp] = xlsread('../../data/appendix/macronewskey_updated.xls','Nonfarm Payroll');

datesNFP = temp(2:end,1);

flagNFP = zeros(length(statementDatesOPEC),1);

for ii=1:length(statementDatesOPEC)
    flagNFP(ii) = sum(strcmp(datesNFP,statementDatesOPEC(ii)));
end

sum(flagNFP)

% Retail Sales
[~,temp] = xlsread('../../data/appendix/macronewskey_updated.xls','Retail Sales');

datesRS = temp(2:end,1);

flagRS = zeros(length(statementDatesOPEC),1);

for ii=1:length(statementDatesOPEC)
    flagRS(ii) = sum(strcmp(datesRS,statementDatesOPEC(ii)));
end

sum(flagRS)

% Industrial Production
[~,temp] = xlsread('../../data/appendix/macronewskey_updated.xls','Industrial Production');

datesIP = temp(2:end,1);

flagIP = zeros(length(statementDatesOPEC),1);

for ii=1:length(statementDatesOPEC)
    flagIP(ii) = sum(strcmp(datesIP,statementDatesOPEC(ii)));
end

sum(flagIP)

% Durable Goods Orders
[~,temp] = xlsread('../../data/appendix/macronewskey_updated.xls','Durable Goods Orders');

datesDGO = temp(2:end,1);

flagDGO = zeros(length(statementDatesOPEC),1);

for ii=1:length(statementDatesOPEC)
    flagDGO(ii) = sum(strcmp(datesDGO,statementDatesOPEC(ii)));
end

sum(flagDGO)

% Trade Balance
[~,temp] = xlsread('../../data/appendix/macronewskey_updated.xls','Trade Balance');

datesTB = temp(2:end,1);

flagTB = zeros(length(statementDatesOPEC),1);

for ii=1:length(statementDatesOPEC)
    flagTB(ii) = sum(strcmp(datesTB,statementDatesOPEC(ii)));
end

sum(flagTB)

% CPI
[~,temp] = xlsread('../../data/appendix/macronewskey_updated.xls','CPI');

datesCPI = temp(2:end,1);

flagCPI = zeros(length(statementDatesOPEC),1);

for ii=1:length(statementDatesOPEC)
    flagCPI(ii) = sum(strcmp(datesCPI,statementDatesOPEC(ii)));
end

sum(flagCPI)

% PPI
[~,temp] = xlsread('../../data/appendix/macronewskey_updated.xls','PPI');

datesPPI = temp(2:end,1);

flagPPI = zeros(length(statementDatesOPEC),1);

for ii=1:length(statementDatesOPEC)
    flagPPI(ii) = sum(strcmp(datesPPI,statementDatesOPEC(ii)));
end

sum(flagPPI)


%% Export table

FID = fopen('../../results/appendix/tablea1.tex','w');
fprintf(FID, strcat('\\begin{tabular}{lrlllr}\\toprule\\midrule  \n'));  
fprintf(FID,' Announcement &  \\multicolumn{1}{c}{Observations} & Source & \\multicolumn{1}{c}{Dates} & Frequency & \\multicolumn{1}{c}{Overlaps}  \\\\ \\midrule  \n');
fprintf(FID,' GDP &  %4.0f \\ \\ \\ \\ \\ & BEA & 4/1987-12/2017 & quarterly & %2.0f \\ \\ \\ \\ \\ \\ \\\\  \n',length(datesGDP),sum(flagGDP));
fprintf(FID,' Unemployment rate &  %4.0f \\ \\ \\ \\ \\ & BLS & 1/1983-12/2017 & monthly & %2.0f \\ \\ \\ \\ \\ \\ \\\\   \n',length(datesUR),sum(flagUR));
fprintf(FID,' Nonfarm payrolls &  %4.0f \\ \\ \\ \\ \\ & BLS & 2/1985-12/2017 & monthly & %2.0f \\ \\ \\ \\ \\ \\ \\\\   \n',length(datesNFP),sum(flagNFP));
fprintf(FID,' Retail sales &  %4.0f \\ \\ \\ \\ \\ & BC & 12/1986-12/2017 & monthly & %2.0f \\ \\ \\ \\ \\ \\ \\\\   \n',length(datesRS),sum(flagRS));
fprintf(FID,' Industrial production &  %4.0f \\ \\ \\ \\ \\ & FRB & 12/1986-12/2017 & monthly & %2.0f \\ \\ \\ \\ \\ \\ \\\\   \n',length(datesIP),sum(flagIP));
fprintf(FID,' Durable goods orders &  %4.0f \\ \\ \\ \\ \\ & BC & 4/1983-12/2017 & monthly & %2.0f \\ \\ \\ \\ \\ \\ \\\\  \n',length(datesDGO),sum(flagDGO));
fprintf(FID,' Trade balance &  %4.0f \\ \\ \\ \\ \\ & BEA & 12/1986-12/2017 & monthly & %2.0f \\ \\ \\ \\ \\ \\ \\\\   \n',length(datesTB),sum(flagTB));
fprintf(FID,' CPI &  %4.0f \\ \\ \\ \\ \\ & BLS & 1/1983-12/2017 & monthly & %2.0f \\ \\ \\ \\ \\ \\ \\\\   \n',length(datesCPI),sum(flagCPI));
fprintf(FID,' PPI &  %4.0f \\ \\ \\ \\ \\ & BLS & 12/1986-12/2017 & monthly & %2.0f \\ \\ \\ \\ \\ \\ \\\\   \n',length(datesPPI),sum(flagPPI));
fprintf(FID,' FOMC &  %4.0f \\ \\ \\ \\ \\ & FED & 3/1983-12/2017 & six-week & %2.0f \\ \\ \\ \\ \\ \\ \\\\  \n',length(statementDatesFed),sum(flagFED));
fprintf(FID, '\\midrule\\bottomrule \n');
fprintf(FID, '\\end{tabular}\n');
fclose(FID);

