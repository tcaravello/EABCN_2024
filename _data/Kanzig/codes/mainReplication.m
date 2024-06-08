%% Replication files for "The macroeconomic effects of oil supply news"
% This is the main shell that can be used to reproduce all results in the
% paper

% All code was written and tested in Matlab R2019b
% The results are saved in ../results of the current directory
% Auxiliary functions are located in /auxfiles

% Diego R. Känzig
% LBS, September 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

% Produce figures 1 and 2: 
s01_figures1_2;

% Produce table 1: 
s02_table1;

% Produce figures 3 and 5:
s03_figures3_5;

% Produce figures 4a and 4b:
s04_figure4a;
s05_figure4b;

% Produce figures 6, 8, 9a, 10, and 11:
s06_figures6_8_9a_10_11;

% Produce figures 7, 9b, and 12:
s07_figures7_9b_12;

% Produce table 2:
s08_table2;

toc