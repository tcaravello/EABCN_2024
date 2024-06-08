function inv_hessian = compute_hessian_at_opt(x,full_inertia,target,target_Sigma_inv,model)


end

function error_fit = handle_full_intertia(x,full_inertia,target,target_Sigma_inv,model)
if full_inertia = 0
    error_fit = solve_model(x,target,target_Sigma_inv,model);
else
    param_aux = zeros(7,1);
    param_aux([1 2 4 6 7]) = x;
    param_aux([3 5]) = 0.99999;
    error_fit = solve_model(param_aux,target,target_Sigma_inv,model);
end