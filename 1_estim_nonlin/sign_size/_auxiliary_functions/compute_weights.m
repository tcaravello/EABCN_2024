function  [weights, weights_se] = compute_weights(a_grid,z,f_z,var_z)
%%% INPUTS
% a_grid n_grid x 1, grid of values of a to compute the weights
% z: T x 1, values of the shock or proxy
% f_z: values of the function of the shock to be used. In general,
% something like z-beta f(z) or f(z)-beta z where beta is the multiple
% regression cofficient of one variable on the other + controls.
% var_z: residual variance of the variable once residualized wrt all controls. 

%%% OUTPUTS
% weights: n_grid x 1, vector of regression weights.
% written by Tomas Caravello 1/23/2023


size_grid = size(a_grid,1);
demeaned_f_z = f_z-mean(f_z);
raw_weights = zeros(size_grid,1);
raw_weights_var = zeros(size_grid,1);
for i = 1:size_grid;
    % define vector of indicators
    ind_vector_z = z>=a_grid(i);
    % compute the weights
    raw_weights(i) = mean(ind_vector_z.*demeaned_f_z);
    raw_weights_var(i) = var((ind_vector_z-mean(ind_vector_z)).*demeaned_f_z);
end
weights = raw_weights/var_z;
weights_se = sqrt(raw_weights_var)/var_z;