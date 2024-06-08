function [c,d] = c_hybrid_fn(exo,A_inv);

global beta omega sigma tau_y Y_SS C_SS D_OLG_SS r_SS ...
    T

% get inputs

y_hat     = exo.y_hat;
i_hat     = exo.i_hat;
pi_hat    = exo.pi_hat;
zeta_hat  = exo.zeta_hat;

T = length(y_hat);

% solve equilibrium system

b = [y_hat + D_OLG_SS/Y_SS * ([0;i_hat(1:T-1)] - pi_hat); ...
    (1-beta*omega) * (1-omega) * y_hat ...
    + (1-beta*omega) * (1-omega) * D_OLG_SS/Y_SS * ([0;i_hat(1:T-1)] - pi_hat) ...
    - (sigma * beta * omega - (1-beta*omega) * (1-omega) * beta * D_OLG_SS/Y_SS) * (i_hat - [pi_hat(2:T);0]) ...
    - sigma * beta * omega * zeta_hat];

sol = A_inv * b;

c = sol(1:T);
d = sol(T+1:2*T);

end