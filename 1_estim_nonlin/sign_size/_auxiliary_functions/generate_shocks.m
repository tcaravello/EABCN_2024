function [base_shocks,f_shock] = generate_shocks(type,order, T_s)

global coef_linear coef_non_linear cut_off prop_zeros sigma_oil T_burn N_sim

% Generate Shocks

obs_shock = sigma_oil*rand(T_burn+T_s,N_sim);
obs_shock_2 = obs_shock>prop_zeros;
shocks = randn(T_burn+T_s,N_sim).*obs_shock_2;


if strcmp(type,'sign')
if order == 1
    f_shock = coef_linear * shocks + coef_non_linear * abs(shocks);
else
    f_shock = coef_linear * shocks + coef_non_linear*(abs(shocks).^(order));
end
elseif strcmp(type,'size')
    if order == 1
        bound = cut_off*sigma_oil;
        fun_size = @(x) ((x>bound).*(x-bound)+(x<-bound).*(x+bound));
        f_shock = coef_linear * shocks + coef_non_linear * fun_size(shocks);
    else
        f_shock = coef_linear * shocks + shocks.*((abs(shocks)).^(order-1));
    end
else
    f_shock = x_moments;
end


% if strcmp(type,'sign')
% if order == 1
%     shocks_plus =  shocks.*(shocks>0); 
%     shocks_minus = shocks.*(shocks<0); 
%     f_shock = sigma_oil*(shocks_plus + small_factor*shocks_minus);
% else
%     f_shock = shocks + (1-small_factor)*(sigma_oil*abs(shocks)).^(order);
% end
% elseif strcmp(type,'size')
%     if order == 1
%     shocks_small =  shocks.*(abs(shocks)<cut_off); 
%     shocks_big = shocks.*(abs(shocks)>=cut_off); 
%     shocks_big_use = shocks_big + cut_off * (shocks_big<=-cut_off) - cut_off * (shocks_big>=cut_off);
%     f_shock = sigma_oil*(shocks_big_use + small_factor*shocks_small);
%     else
%         f_shock = shocks + (1-small_factor)*shocks.*(sigma_oil*abs(shocks)).^(order-1);
%     end
% else
%     f_shock = sigma_oil*shocks;
% end
base_shocks = shocks(T_burn+1:end, :);


end