function f_shock = transform_shocks(shocks,type,order)

global small_factor cut_off prop_zeros sigma_oil

if strcmp(type,'sign')
if order == 1
    shocks_plus =  shocks.*(shocks>0); 
    shocks_minus = shocks.*(shocks<0); 
    f_shock = sigma_oil*(shocks_plus + small_factor*shocks_minus);
else
    f_shock = (sigma_oil*abs(shocks)).^(order);
end
elseif strcmp(type,'size')
    if order == 1
    shocks_small =  shocks.*(abs(shocks)<cut_off); 
    shocks_big = shocks.*(abs(shocks)>=cut_off); 
    shocks_big_use = shocks_big + cut_off * (shocks_big<=-cut_off) - cut_off * (shocks_big>=cut_off);
    f_shock = sigma_oil*(shocks_big_use + small_factor*shocks_small);
    else
        f_shock = shocks.*(sigma_oil*abs(shocks)).^(order-1);
    end
else
    f_shock = shocks;
end



end