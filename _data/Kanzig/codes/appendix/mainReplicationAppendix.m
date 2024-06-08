%% Replication files for the online appendix of "The macroeconomic effects of oil supply news"
% This is the main shell that can be used to reproduce all results in the
% online appendix

% All code was written and tested in Matlab R2019b
% The results are saved in ../results of the current directory
% Auxiliary functions are located in /auxfiles

% Note that not all data could be publicly shared. The files that do not
% run without this data have been commented out.

% Diego R. Känzig
% LBS, September 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

% Produce table A.1
a01_tablea1;

% Produce table A.3
a02_tablea3;

% Produce figure A.2 (Panels A & B)
a03_figurea2a;
a04_figurea2b;

% Produce figure A.4 (Panels A & B)
a05_figurea4a;
a06_figurea4b;

% Produce figure A.5
a07_figurea5;

% Produce figure A.6
a08_figurea6;

% Produce figure A.9
a09_figurea9;

% Produce figure A.11
a10_figurea11;

% Produce figure A.12
a11_figurea12;

% Produce figure A.13 
% a12_figurea13;

% Produce figure A.14
a13_figurea14;

% Produce figure A.16
a14_figurea16;

% Produce figure A.17
a15_figurea17;

% Produce figure A.18
a16_figurea18;

% Produce figure A.19
a17_figurea19;

% Produce figures A.20-A.23
a18_figuresa20_21_22_23;

% Produce figure A.24
a19_figurea24;

% Produce figure A.25
a20_figurea25;

% Produce figure A.26
a21_figurea26;

% Produce figures A.27-A.30
a22_figuresa27_28_29_30;

close all;

toc
