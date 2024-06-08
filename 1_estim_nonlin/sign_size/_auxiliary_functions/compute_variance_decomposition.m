n_x_mom  = 1000000;

% generate a simulated sample to compute moments
if strcmp(distribution_simul,'normal')
    x_moments = randn(n_x_mom,1);
elseif strcmp(distribution_simul,'exponential')
    x_moments = exprnd(ones(n_x_mom,1))-1;
end





terms_oil = prop_big *(sigma_oil * ma_coefs_vec_aux).^2 + (1-prop_big) *(small_factor*sigma_oil * ma_coefs_vec_aux).^2;
terms_confund =(sigma_c * coefs_confound_aux).^2;

for i_hor = 1:length(hor_var)
    hor_var_aux = hor_var(i_hor);
    oil_term_aux = sum(terms_oil(1:hor_var_aux));
    confound_term_aux = sum(terms_confund(1:hor_var_aux));
    var_share(i_hor) = oil_term_aux/(oil_term_aux + confound_term_aux);
end