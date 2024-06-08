% this code was in the main file.
% get hessian. THIS IS NUMERICALLY UNSTABLE. WILL CHANGE TO ADAPTATIVE SAMPLING
% x0_hessian = param_sol([1 2 4 6 7]);
% options_min_hessian = optimset('Display','iter','MaxFunEvals',10000,'Algorithm','trust-region');
% [x_sol_2,fval,exitflag,output_2,grad,hessian] = compute_hessian_at_opt(x0_hessian,target,target_Sigma_inv,model,options_min_hessian);
% inv_hessian = inv(hessian);


function [x,fval,exitflag,output_2,grad,hessian] = compute_hessian_at_opt(param_sol,target,target_Sigma_inv,model,options_min)

[x,fval,exitflag,output_2,grad,hessian] = fminunc(@(x) handle_full_intertia(x,target,target_Sigma_inv,model),param_sol, options_min);

end

function error_fit = handle_full_intertia(x,target,target_Sigma_inv,model)
if length(x) == 7
    error_fit = solve_model(x,target,target_Sigma_inv,model);
else
    param_aux = zeros(7,1);
    param_aux([1 2 4 6 7]) = x;
    param_aux([3 5]) = 0.9999999;
    error_fit = solve_model(param_aux,target,target_Sigma_inv,model);
end
end