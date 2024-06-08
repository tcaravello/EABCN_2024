% import VAR data
load var_data.mat

date         = data_var(:,1);
tfp         = data_var(:,2);
lpcom       = data_var(:,3);
rr_shock    = data_var(:,4);
gk_0        = data_var(:,5);
gk_1        = data_var(:,6);
gk_2        = data_var(:,7);
gk_3        = data_var(:,8);
gk_4        = data_var(:,9);
mar_raw     = data_var(:,10);
mar_svariv  = data_var(:,11);
ad_shock    = data_var(:,12);

unemp        = data_var(:,13);
gdp_per_cap  = 100*data_var(:,14);
dur_cons     = 100*data_var(:,15);
inv          = 100*data_var(:,16);
non_dur_cons = 100*data_var(:,17);
serv_cons    = 100*data_var(:,18);
cons         = 100*log(exp(data_var(:,15))+exp(data_var(:,17))+exp(data_var(:,18)));
cons_other   = 100*log(exp(data_var(:,17))+exp(data_var(:,18)));
lab          = data_var(:,19);
lab_prod    = data_var(:,20);
lab_share   = data_var(:,21); 
inf         = data_var(:,22);
ffr         = data_var(:,23);
tb3         = data_var(:,24);
cpi_inf     = data_var(:,25);
cpi_level   = cumsum(data_var(:,25)); 
rpce        = 100*data_var(:,26);
pce_inf     = data_var(:,27);
real_oil    = data_var(:,28); 

gdp_tot     = 100*data_var(:,29); %Real GDP
gdp_growth = 4*[0;(gdp_tot(2:end)-gdp_tot(1:end-1))]; %multiply times 4 to annualize.
gdp_tot_filtered = stat_transform(gdp_tot,1);

inv_gap = stat_transform(inv,1);
cons_gap = stat_transform(cons,1);
lab_gap       = 100 * stat_transform(lab,1);
lab_prod_gap  = 100 * stat_transform(lab_prod,1);
lab_share_gap = 100 * stat_transform(lab_share,1);
tfp_gap       = nan(length(tfp),1); 
tfp_gap(1:end-2)  = 100 * stat_transform(tfp(1:end-2),1);

[inf_detrended,inf_trend] = regcyc_2(pce_inf);
%output_detrended = stat_transform(gdp_tot,4); %usa quadratic determinisitc trend.
%[output_gap, nat_output] = regcyc_2(output_detrended); % separate the remainder in "natural" and "gap".

[output_gap, nat_output] = regcyc_2(gdp_tot); % separate in "natural" and "gap".

dur_cons_tot     = 100*data_var(:,30);
inv_tot          = 100*data_var(:,31);
non_dur_cons_tot = 100*data_var(:,32);
serv_cons_tot    = 100*data_var(:,33);
cons_tot         = 100*log(exp(data_var(:,30))+exp(data_var(:,32))+exp(data_var(:,33)));
cons_other_tot  = 100*log(exp(data_var(:,32))+exp(data_var(:,33)));
pce_core        = data_var(:,34);
cpi_shelter     = data_var(:,35); 
cpi_non_shelter = data_var(:,36);
cpi_energy = data_var(:,37);
pce_inf_shelter = data_var(:,38);
pce_inf_non_shelter = data_var(:,39);

cpi_food = data_var(:,40);
cpi_nondur = data_var(:,41);
cpi_dur = data_var(:,42);
cpi_serv = data_var(:,43);

mort_30yr = data_var(:,44);
bond_30yr = data_var(:,45);
bond_10yr = data_var(:,46);
auto_4yr  = data_var(:,47); 
bond_5yr  = data_var(:,48);
ebp = data_var(:,49);
fci = data_var(:,50);
fci_2 = data_var(:,51);
fci_goldman = data_var(:,52);
fci_2_1y = data_var(:,53);
pot_gdp = 100*data_var(:,54);
r_star_lm = data_var(:,55);
r_star_hlw = data_var(:,56);

output_gap_cbo = gdp_tot-pot_gdp;
fci_nan_index = find(isnan(fci_2),1,'last');
fci_3 = fci_2;
fci_3(1:fci_nan_index) = fci(1:fci_nan_index);

fci_3_1y = fci_2_1y;
fci_3_1y(1:fci_nan_index) = fci(1:fci_nan_index);

