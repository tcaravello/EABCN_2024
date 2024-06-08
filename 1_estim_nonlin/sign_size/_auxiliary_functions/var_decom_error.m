function [error, varargout] = var_decom_error(x)

global var_decomp_target var_decomp_weights hor_coefs ma_coefs_vec_aux var_x_mom hor_var

rho_confound = x(1); %rho of the confounder
sigma_c = x(2); %unconditional standard deviation of the confounder shock
 

coefs_confound_aux = rho_confound.^(hor_coefs);

terms_oil = var_x_mom*(ma_coefs_vec_aux).^2;
terms_confund =(sigma_c * coefs_confound_aux).^2;

var_share = zeros(length(hor_var),1);

for i_hor = 1:length(hor_var)
    hor_var_aux = hor_var(i_hor);
    oil_term_aux = sum(terms_oil(1:hor_var_aux));
    confound_term_aux = sum(terms_confund(1:hor_var_aux));
    var_share(i_hor) = oil_term_aux/(oil_term_aux + confound_term_aux);
end

error = sum(var_decomp_weights.*(var_share-var_decomp_target).^2);

if nargout>1
    varargout{1} = var_share;
end