function calib_error = calib_fn_aux(param);

global beta MPC_target

% set parameters

omega_1 = param(1);
omega_2 = param(2);
chi_1   = param(3);
chi_2   = param(4);

% compute quarterly MPCs

max_hor = 4 * 6;

MPC_path = zeros(max_hor,1);
MPC_path(1) = chi_1 * (1 - beta * omega_1) + chi_2 * (1 - beta * omega_2) ...
    + (1 - chi_1 - chi_2);
for i_hor = 2:max_hor
    MPC_path(i_hor) = chi_1 * (1 - beta * omega_1) * omega_1^(i_hor-1) ... 
        + chi_2 * (1 - beta * omega_2) * omega_2^(i_hor-1);
end

beta_vec = zeros(max_hor,1);
for i_hor = 1:max_hor
    beta_vec(i_hor) = beta^(i_hor-1);
end

MPC_NPV = beta_vec .* MPC_path;
cons_share_annual = zeros(max_hor/4,1);
for t = 1:max_hor/4
    cons_share_annual(t) = sum(MPC_NPV(1+(t-1)*4:t*4,1));
end
cons_share_annual = cumsum(cons_share_annual);

% compute discrepancy

calib_error = sqrt(sum((cons_share_annual - MPC_target).^2));