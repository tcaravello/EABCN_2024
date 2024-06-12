%% GENERATE SPENDER-OLG-HYBRID DERIVATIVE MATRICES

% Tomas Caravello and Christian Wolf

% this version: June 8, 2024

%% HOUSEKEEPING

clc
clear all
close all


path = '/Users/tomyc/Dropbox (MIT)/EACBN_2024/code';
session = '/3_micro_macro';
experiment = '/olg';

save_fig = 0;
only_final_figure = 1; % 0 if you want figures for slides, 1 if you only want final figure.


addpath([path session experiment '/_auxiliary_functions']);
cd([path session experiment]);

%% PARAMETERS

%----------------------------------------------------------------
% Preferences
%----------------------------------------------------------------

global beta omega sigma MPC_target mu

beta  = 0.99^(0.25);
omega = 0.88;
sigma = 1;

MPC_target = 0.22;
mu         = (MPC_target - (1 - beta*omega))/(beta*omega);

%----------------------------------------------------------------
% Steady-State Objects
%----------------------------------------------------------------

global tau_y Y_SS C_SS r_SS D_SS D_OLG_SS

tau_y    = 1/3;
Y_SS     = 1;
C_SS     = Y_SS;
r_SS     = 1/beta - 1;
D_SS     = 1.04;
D_OLG_SS = 1/(1-mu) * D_SS;

%% SETTINGS

global T

T    = 500;
step = 1;

exo.y_hat     = zeros(T,1);
exo.i_hat     = zeros(T,1);
exo.pi_hat    = zeros(T,1);
exo.zeta_hat  = zeros(T,1);

r_NPV = zeros(T,1);
for t = 1:T
    r_NPV(t) = (1/(1+r_SS))^(t-1);
end

%% DERIVATIVE MATRICES

%----------------------------------------------------------------
% Auxiliary Matrix
%----------------------------------------------------------------

A = zeros(2*T,2*T); % order (c_t,a_{t+1}) & (BC, EE)

% Budget Constraint

for t = 1:T
    A(t,t) = 1; % c_{t}
    A(t,T+t) = beta; % beta a_{t+1} 
    if t > 1
        A(t,T+t-1) = -1; % - a_{t}
    end
end

% Euler Equation

for t = T+1:2*T
    % in the first equation, the previous periods assets are given.
    if t == T+1
        A(t,t-T) = 1 - omega * (1-beta*omega); % (1 - omega (1 - beta omega)) c_t
        A(t,t-T+1) = - beta * omega; % - beta omega c_{t+1}
    elseif t < 2*T
        A(t,t-T) = 1 - omega * (1-beta*omega); % (1 - omega (1 - beta omega)) c_t
        A(t,t-T+1) = - beta * omega; % - beta omega c_{t+1}
        A(t,t-1) = - (1-beta*omega) * (1-omega); % - (1-beta omega)(1-omega) a_{t} 
    elseif t == 2*T
        A(t,t-T) = 1 - omega * (1-beta*omega);
        A(t,t-1) = 1; %impose zero assets at T
    end
end

A_inv = A^(-1);

%----------------------------------------------------------------
% Baseline
%----------------------------------------------------------------

[c_base,d_base] = c_hybrid_fn(exo,A_inv);

%----------------------------------------------------------------
% Income
%----------------------------------------------------------------

C_y = NaN(T,T);
D_y = NaN(T,T);

for t = 1:T
    exo.y_hat(t) = step;
    [c_shock,d_shock] = c_hybrid_fn(exo,A_inv);
    C_y(:,t) = (c_shock-c_base)/step;
    D_y(:,t) = (d_shock-d_base)/step;
    exo.y_hat(t) = 0;
end

% add the spender fraction
C_y = mu * eye(T) + (1-mu) * C_y;
D_y = mu * zeros(T) + (1-mu) * D_y;
%% BEHAVIORAL ADJUSTMENTS

%----------------------------------------------------------------
% Select Friction
%----------------------------------------------------------------

hh_R      = 0; % RE
hh_myopic = 0; % myopic
hh_CD     = 1; % cognitive discounting
hh_SI     = 0; % sticky information

%----------------------------------------------------------------
% Construct Expectations Matrices
%----------------------------------------------------------------

% rational expectations

E_R = ones(T,T);

% full myopia

E_M = zeros(T,T);

for i = 1:T
    for j = 1:i
        E_M(i,j) = 1;
    end
end

% cognitive discounting

theta_CD = 0.8;

E_CD = zeros(T,T);

for i = 1:T
    for j = 1:T
        if j < i
            E_CD(i,j) = 1;
        else
            E_CD(i,j) = theta_CD^(j-i);
        end
    end
end

% sticky information

theta_SI = 0.95;

E_SI = zeros(T,T);

for i = 1:T
    for j = 1:T
        if j <= i
            E_SI(i,j) = 1;
        else
            E_SI(i,j) = 1 - theta_SI^i;
        end
    end
end

% select

if hh_R == 1

    E = E_R;
    behav_string = '';
    
elseif hh_myopic == 1
    
    E = E_M;
    behav_string = '_myopic';
elseif hh_CD == 1
    
    E = E_CD;
    behav_string = '_cog_disc';

elseif hh_SI == 1
    
    E = E_SI;
    behav_string = '_sticky_info';
end

%----------------------------------------------------------------
% Adjust Derivative Matrices
%----------------------------------------------------------------

% income

C_y_be = ME_fn(C_y,E);
D_y_be = ME_fn(D_y,E);
%% REPORT STATISTICS FOR MODEL CALIBRATION

%----------------------------------------------------------------
% iMPCs
%----------------------------------------------------------------

disp(['The impact MPC is ' num2str(C_y(1,1))])
disp(['The iMPC decay rate is ' num2str(C_y(3,1)/C_y(2,1))])

%----------------------------------------------------------------
% Anticipation Effect Decay
%----------------------------------------------------------------

disp(['The iMPC anticipation decay rates are ' num2str(C_y(1,2)/C_y(1,1)) ' and ' num2str(C_y(1,3)/C_y(1,2))])

%----------------------------------------------------------------
% Spending Share
%----------------------------------------------------------------

cons_share = cumsum(r_NPV .* C_y(:,1));
disp(['The fraction of income spent within 10 years is ' num2str(cons_share(40))])

cons_NPV = r_NPV .* C_y(:,1);
cons_share_annual = zeros(T/4,1);

for t = 1:T/4
    cons_share_annual(t) = sum(cons_NPV(1+(t-1)*4:t*4,1));
end
cons_share_annual = cumsum(cons_share_annual);

% empirical evidence on iMPCs from Fagergan-Holm-Natvik
load fhn_results

%% PLOT iMPCs

close all

% colors

settings.colors.black  = [0 0 0];
settings.colors.grey   = [230/255 230/255 230/255];
settings.colors.orange = [204/255 102/255 0/255];
settings.colors.green = [37/255 152/255 14/255];
settings.colors.navyblue = [0/255 0/255 50/255];
settings.colors.lnavyblue = 0.25 * [0/255 0/255 50/255] + 0.75 * [1 1 1];
settings.colors.purple = [160/255 32/255 240/255];
settings.colors.lpurple = 0.25 * [160/255 32/255 240/255] + 0.75 * [1 1 1];
settings.colors.beige  = [116/255 158/255 178/255];
settings.colors.lbeige = 0.25 * settings.colors.beige + 0.75 * [1 1 1];

weights = [1 0.8 0.5 0.2];
settings.colors.orange_all = zeros(length(weights),3);
settings.colors.blue_all = zeros(length(weights),3);
for i = 1:4
    settings.colors.orange_all(i,:) = weights(i) * settings.colors.orange + (1-weights(i)) * [1 1 1];
    settings.colors.blue_all(i,:) = weights(i) * settings.colors.black + (1-weights(i)) * [1 1 1];
end

% plot size

plotwidth = 0.385;
gapsize = 0.055;
gapsize_edges = (1-2*plotwidth-1*gapsize)/2;
left_pos = [gapsize_edges, gapsize_edges + gapsize + plotwidth];

% plot iMPC matrix and cumulative MPCs

select_hor_aux = [1 6 11 16];

max_hor = 21;

cd([path session experiment '/_figures']);

if only_final_figure == 1
    figs_indicator = [2];
else
    figs_indicator = [1;2];
end

% Full Information

for fig_n = figs_indicator

if fig_n == 1
select_hor = select_hor_aux(1);
else
select_hor = select_hor_aux;
end

figure(fig_n)

subplot(1,2,1)
pos = get(gca, 'Position');
pos(1) = left_pos(1);
pos(3) = plotwidth;
set(gca,'Position', pos)
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex')
hold on
for i = 1:length(select_hor)
    plot(0:1:max_hor-1,C_y(1:max_hor,select_hor(i)),'linewidth',5,'linestyle','-','color',settings.colors.blue_all(i,:))
    hold on
end
set(gcf,'color','w')
title('Intertemporal MPCs','interpreter','latex','fontsize',21)
xlabel('Quarter','interpreter','latex','FontSize',18)
ylabel('MPC','interpreter','latex','FontSize',18)%,'Rotation',0)
legend({'$t=0$','$t=5$','$t=10$','$t=15$'},'NumColumns',2,'Location','Southeast','fontsize',18,'interpreter','latex')
ylim([-0.04 0.25])
grid on
hold off
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) pos(3) pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');

subplot(1,2,2)
pos = get(gca, 'Position');
pos(1) = left_pos(2);
pos(3) = plotwidth;
set(gca,'Position', pos)
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex')
hold on
jbfill(-4:1:5,lowCI,uppCI,...
    0.9 * [1 1 1],0.9 * [1 1 1],0,1);
hold on
% plot(-4:1:5,cumC1,'Color',settings.colors.beige,'LineWidth',4)
% hold on
plot(-4:1:5,[zeros(4,1);cons_share_annual(1:6)],'linewidth',5,'linestyle','-','color',settings.colors.blue_all(2,:))
hold off
set(gcf,'color','w')
title('Cumulative MPCs','interpreter','latex','fontsize',21)
xlabel('Year','interpreter','latex','FontSize',18)
% ylabel('MPC','interpreter','latex','FontSize',18)%,'Rotation',0)
legend({'Data','Model'},'Location','Southeast','fontsize',18,'interpreter','latex')
ylim([-0.19 1])
xlim([-4 5])
grid on
hold off

pos = get(gcf, 'Position');
set(gcf, 'Position', [0.1*pos(1) 0.5*pos(2) 2.2*pos(3) 1.5*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
if save_fig == 1
if fig_n == 1
    print('iMPCs_hybrid_1','-dpng');
else
    print('iMPCs_hybrid_2','-dpng');
end
end
end
%%
% Behavioral

if only_final_figure == 1
    figs_indicator = [2];
else
    figs_indicator = [1;2];
end

for fig_n = figs_indicator

if fig_n == 1
select_hor = select_hor_aux(1);
else
select_hor = select_hor_aux;
end

figure(fig_n+2)

subplot(1,2,1)
pos = get(gca, 'Position');
pos(1) = left_pos(1);
pos(3) = plotwidth;
set(gca,'Position', pos)
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex')
hold on
for i = 1:length(select_hor)
    plot(0:1:max_hor-1,C_y_be(1:max_hor,select_hor(i)),'linewidth',5,'linestyle','-','color',settings.colors.blue_all(i,:))
    hold on
end
set(gcf,'color','w')
title('Intertemporal MPCs','interpreter','latex','fontsize',21)
xlabel('Quarter','interpreter','latex','FontSize',18)
ylabel('MPC','interpreter','latex','FontSize',18)%,'Rotation',0)
legend({'$t=0$','$t=5$','$t=10$','$t=15$'},'NumColumns',2,'Location','Southeast','fontsize',18,'interpreter','latex')
ylim([-0.04 0.25])
grid on
hold off
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) pos(3) pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');

subplot(1,2,2)
pos = get(gca, 'Position');
pos(1) = left_pos(2);
pos(3) = plotwidth;
set(gca,'Position', pos)
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex')
hold on
jbfill(-4:1:5,lowCI,uppCI,...
    0.9 * [1 1 1],0.9 * [1 1 1],0,1);
hold on
% plot(-4:1:5,cumC1,'Color',settings.colors.beige,'LineWidth',4)
% hold on
plot(-4:1:5,[zeros(4,1);cons_share_annual(1:6)],'linewidth',5,'linestyle','-','color',settings.colors.blue_all(2,:))
hold off
set(gcf,'color','w')
title('Cumulative MPCs','interpreter','latex','fontsize',21)
xlabel('Year','interpreter','latex','FontSize',18)
% ylabel('MPC','interpreter','latex','FontSize',18)%,'Rotation',0)
legend({'Data','Model'},'Location','Southeast','fontsize',18,'interpreter','latex')
ylim([-0.19 1])
xlim([-4 5])
grid on
hold off
pos = get(gcf, 'Position');
set(gcf, 'Position', [0.1*pos(1) 0.5*pos(2) 2.2*pos(3) 1.5*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
if save_fig == 1
if fig_n == 1
    print(['iMPCs_hybrid_1' behav_string],'-dpng');
else
    print(['iMPCs_hybrid_2' behav_string],'-dpng');
end
end

end


cd([path session experiment]);