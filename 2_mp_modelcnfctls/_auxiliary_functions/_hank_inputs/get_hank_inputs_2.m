%% GET HANK INPUTS FOR POLICY COUNTERFACTUAL ILLUSTRATION
% Tomas Caravello & Alisdair McKay & Christian Wolf
% this version: 07/10/2023

%% HOUSEKEEPING

%clc
%clear all
%close all

warning('off','MATLAB:dispatcher:nameConflict')
path_2 = [path session '/_auxiliary_functions/_hank_inputs']; 
addpath(genpath([path_2 '/_aux']))
addpath([path_2 '/_income_process'])
addpath([path_2 '/_steady_state'])
addpath([path_2 '/_jacobians'])

cd(path_2);

%% Calibration

%global beta alpha gamma delta varphi tau_l_rate A_SS r_SS L_SS Y_SS I_SS W_SS D_SS Tau_SS G_SS C_SS B_Y_SS Tau_Y

% beta        = 0.99;
% alpha       = 1/3;
% eis	        = 1;
% delta	    = 0.03;
% varphi	    = 1;
% 
% tau_l_rate  = 0.3;
% Tau_Y       = -0.05;
% B_Y_SS      = 1.04; 
% 
% r_SS	= 1/beta - 1;
% L_SS	= 1;
% K_SS	= (alpha/(r_SS + delta))^(1/(1-alpha)) * L_SS;
% Y_SS	= K_SS^alpha * L_SS^(1-alpha);
% I_SS	= delta * K_SS;
% W_SS	= (1-alpha) * K_SS^alpha * L_SS^(-alpha);
% D_SS	= Y_SS - W_SS * L_SS - I_SS;
% B_SS    = B_Y_SS* Y_SS;
% A_SS	= B_SS;
% Tau_SS = Tau_Y * Y_SS;
% G_SS = tau_l_rate*(1-alpha)*Y_SS + Tau_SS - r_SS * A_SS;
% C_SS	= Y_SS - I_SS - G_SS;


%% COMPUTE STEADY STATE

disp('I am solving for the steady state of the model.')

get_ss

disp('Done!')

%% COMPUTE JACOBIANS

disp('I am computing the PE Jacobians.')

get_jacobians

disp('Done!')
%%
load inputs_hank.mat
