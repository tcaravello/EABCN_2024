function M_tilde = ME_fn(M,E);

T = size(M,1);

M_tilde = 0 * M;

for t = 1:T
    for s = 1:T
        for tau = 1:min(t,s)
            if tau == 1
                M_tilde(t,s) = M_tilde(t,s) + (E(tau,s) - 0) * M(t-tau+1,s-tau+1);
            else
                M_tilde(t,s) = M_tilde(t,s) + (E(tau,s) - E(tau-1,s)) * M(t-tau+1,s-tau+1);
            end
        end
    end
end

end