function  [irf_coefs,irf_se,irf_t_stats,irf_r_sq,residualized_z] = run_local_projections_predictability(y,x,lags_controls,horizon, level, opt)

%%% INPUTS
% y: T x n_y, output variables to explain at different horizons
% z: T x n_f, matrix of observations for the shock and non-linear transformations 
% x: T x n_c, matrix of variables to use as controls. DO NOT INCLUDE A
% CONSTANT, WILL BE ADDED LATER.
% lags_controls: 1 x 1, integer. Number of lags to use for the control
% variables.
% horizon: 1 x 1, integer. Number of months ahead to compute the IRF
% level: 1 x 1 logical, =1 if you want LHS variable to be the level y_{t+h},
%                       =0 if you want LHS variable to be y_{t+h}-y_{t-1},
%
% opt: 1 x 1, logical, = 1
%%% OUTPUTS
% irf_coefs = coefficents of the predictability regressions;
% irf_se = se of the predictability regressions
% irf_t_stats = t-stats of the predictability regressions
% irf_r_sq = R^2 of the predictability regerssions;
% residualized_z = residualized z and f(z);

T_x = size(x,1);
[T_y,n_y] = size(y);
if (T_x ~= T_y); error('y and x have different number of observations'); end;

%y_leads = flip(lagmatrix(y, -horizon:0),2); %order: first column is y_t, second y_t+1 and so on.
x_mat = [ones(T_x,1), lagmatrix(x,1:lags_controls)];
x_use = [x_mat(lags_controls+1:end,:)];
y_use = y(lags_controls+1:end,:);

irf_coefs = [];
irf_se = [];
irf_t_stats =[];
irf_r_sq = [];
residualized_z = [];

for j=1:n_y %for each output variable
    for i=1:horizon % for each horizon
         if level==1
            yy=y_use(i:end,j); %use the data in levels.
         else
            yy=(y_use(i+1:end,j)-y_use(1:end-i,j)); %use the long difference
         end
         results=nwest(yy, x_use(1:end-i+1,:),i);
         betas(:,i)=results.beta;
         rsquared(:,i) = results.rsqr;
         resid(:,i) = results.resid;

         if opt==0
            se(:,i)=results.se';
         else
            [EstCov, hacse, coeff]=hac_alt(x_use(1:end-i+1,:), yy, 'intercept', false, 'smallT', false, 'display', 'off'); 
            se(:,i)=hacse';
         end
         t_stats(:,i) = betas(:,i)./se(:,i);
    end
% store results for each variable
irf_coefs = [irf_coefs;betas];
irf_se = [irf_se; se];
irf_t_stats = [irf_t_stats;t_stats];
irf_r_sq = [irf_r_sq,rsquared];
residualized_z = [residualized_z,resid];
end

end



