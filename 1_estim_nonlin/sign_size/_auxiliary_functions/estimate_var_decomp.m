n_x_mom  = 1000000;

% generate a simulated sample to compute moments

obs_shock = rand(1000000,1);
obs_shock_2 = obs_shock>prop_zeros;

if strcmp(distribution_simul,'normal')
    x_moments_aux = sigma_oil * randn(n_x_mom,1);
    x_moments = obs_shock_2.*x_moments_aux;
elseif strcmp(distribution_simul,'exponential')
    x_moments_aux = sigma_oil * exprnd(ones(n_x_mom,1))-1;
    x_moments = obs_shock_2.*x_moments_aux;
end


if strcmp(type,'sign')
if order_sim == 1
    f_shock = coef_linear * x_moments + coef_non_linear * abs(x_moments);
else
    f_shock = coef_linear * x_moments + coef_non_linear*(abs(x_moments).^(order_sim));
end
elseif strcmp(type,'size')
    if order_sim == 1
        bound = cut_off*sigma_oil;
        fun_size = @(x) ((x>bound).*(x-bound)+(x<-bound).*(x+bound));
        f_shock = coef_linear * x_moments + coef_non_linear * fun_size(x_moments);
    else
        f_shock = coef_linear * x_moments + x_moments.*(abs(x_moments).^(order_sim-1));
    end
else
    f_shock = x_moments;
end

var_x_mom = var(f_shock);

[x_solve, f_solve] = fminsearch(@(x) var_decom_error(x),[0.999; 0.5]);
%%

rho_confound = x_solve(1); %rho of the confounder
sigma_c = x_solve(2); %unconditional standard deviation of the confounder shock
coefs_confound_aux = rho_confound.^(hor_coefs);

[error, var_shares] = var_decom_error(x_solve);