function [results] = bootstrap_weights(z,f_z,nrep,a_grid,qe,linear_comb_mat)
% this function computes confidence intervals for the weights derived from
% z via wild bootstrap. Bootstrap validity is established by noting that,
% since the weights are a sample average, under standard boundedness
% assumptions on moments the central limit theorem holds, so the
% distribution is asymptotically normal.
% INPUTS
% z: T x 1, vector of values of the shock
% f_z: T x n_func, vector of values of the function of shock. Must be in the
% same order as z (i.e, each element in f_z is in the same position as it's
% corresponding element in z)
% nrep: 1 x 1, number of bootstrap repetitions
% a_grid: size_grid x 1: values of a to compute the weights
% qe: column vector of significance levels to use. Usually [.1,.05]
% linear_comb_mat: m x x n_func: matrix of lienar combinations of coefficeints to obtain weights and bootstrapped CI 


T_tot = size(z,1);
n_funcs = size(f_z,2); %number of different functions to use. 
size_grid = size(a_grid,1);
R_array = unidrnd(T_tot,T_tot,nrep); %bootstrap with resampling
weights_mat = zeros(size_grid,nrep,n_funcs);

for fn = 1:n_funcs
for n = 1:nrep
[weights_mat(:,n,fn), ~]  = compute_weights(a_grid,z(R_array(:,n),1),f_z(R_array(:,n),fn),var(f_z(R_array(:,n),fn)));
end
end
weights_mat_even = (weights_mat+flip(weights_mat,1))/2;
weights_mat_odd = (weights_mat-flip(weights_mat,1))/2;

n_lin = size(linear_comb_mat,1); % number of linear combinations to use
weights_transformed = zeros(size_grid,nrep,n_lin); %weights on the linear transformations
for n_l =1:n_lin
    for fn=1:n_funcs
        weights_transformed(:,:,n_l) = weights_transformed(:,:,n_l)+ weights_mat(:,:,fn)*linear_comb_mat(n_l,fn);
    end
end
weights_transformed_even = (weights_transformed+flip(weights_transformed,1))/2;
weights_transformed_odd = (weights_transformed-flip(weights_transformed,1))/2;


% compute quantiles
n_quantiles = size(qe,1);
quantiles_to_use = zeros(2*n_quantiles,1);
for jj = 1:n_quantiles
quantiles_to_use(2*jj-1:2*jj,1) = [qe(jj);1-qe(jj)] ;
end
n_quant = size(quantiles_to_use,1);
q_total = zeros(size_grid,n_quant,n_funcs);
q_even = zeros(size_grid,n_quant,n_funcs);
q_odd = zeros(size_grid,n_quant,n_funcs);
q_total_trans = zeros(size_grid,n_quant,n_lin);
q_even_trans = zeros(size_grid,n_quant,n_lin);
q_odd_trans = zeros(size_grid,n_quant,n_lin);

for fn = 1:n_funcs
q_total(:,:,fn) = (quantile(weights_mat(:,:,fn)',quantiles_to_use))';
q_even(:,:,fn)  = (quantile(weights_mat_even(:,:,fn)',quantiles_to_use))';
q_odd(:,:,fn)  = (quantile(weights_mat_odd(:,:,fn)',quantiles_to_use))';
end
for n_l = 1:n_lin
q_total_trans(:,:,n_l) = (quantile(weights_transformed(:,:,n_l)',quantiles_to_use))';
q_even_trans(:,:,n_l) = (quantile(weights_transformed_even(:,:,n_l)',quantiles_to_use))';
q_odd_trans(:,:,n_l) = (quantile(weights_transformed_odd(:,:,n_l)',quantiles_to_use))';
end

% compute weights in the actual sample
weights_actual = zeros(size_grid,n_funcs);
for fn = 1:n_funcs
[weights_actual(:,fn), ~] = compute_weights(a_grid,z,f_z(:,fn),var(f_z(:,fn)));
end
weights_actual_even = (weights_actual+flip(weights_actual,1))/2;
weights_actual_odd = (weights_actual-flip(weights_actual,1))/2;

weights_actual_trans = zeros(size_grid,n_lin);
for n_l = 1:n_lin
    for fn = 1:n_funcs
        weights_actual_trans(:,n_l) = weights_actual_trans(:,n_l) + linear_comb_mat(n_l,fn)*weights_actual(:,fn);
    end
end
weights_actual_even_trans = (weights_actual_trans+flip(weights_actual_trans,1))/2;
weights_actual_odd_trans = (weights_actual_trans-flip(weights_actual_trans,1))/2;


results = struct;
results.weights = weights_actual;
results.weights_even = weights_actual_even;
results.weights_odd = weights_actual_odd;
results.weights_trans = weights_actual_trans;
results.weights_even_trans = weights_actual_even_trans;
results.weights_odd_trans = weights_actual_odd_trans;
results.ci_total = q_total;
results.ci_even = q_even;
results.ci_odd = q_odd;
results.ci_total_trans = q_total_trans;
results.ci_even_trans = q_even_trans;
results.ci_odd_trans = q_odd_trans;


end