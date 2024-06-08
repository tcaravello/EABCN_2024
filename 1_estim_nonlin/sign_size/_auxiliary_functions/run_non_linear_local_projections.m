function  [irf_coefs,irf_se,irf_t_stats,irf_var_mat,extra_stuff] = run_non_linear_local_projections(y,z,x,lags_controls,horizon, non_linearity ,varargin)
%%% INPUTS
% y: T x n_y, output variables to explain at different horizons
% z: T x 1, matrix of observations for the shock 
% x: T x n_c, matrix of variables to use as controls in all equations. DO NOT INCLUDE A
% CONSTANT, WILL BE ADDED LATER. Can be left empty.

% lags_controls: 1 x 1, integer. Number of lags to use for the control
% variables.

% horizon: 1 x 1, integer. Number of months ahead to compute the IRF

% non_linearity: string. Determines if you want sign or size.

% TO RUN THE LP WITH ALL LAGS, SPECIFY X = Y. IF YOU WANNA USE OWN LAGS
% ONLY, SET X = []


% options: struct, including: level, opt, har, trend, order

% level: 1 x 1 logical, =1 if you want LHS variable to be the level y_{t+h},
%                       =0 if you want LHS variable to be y_{t+h}-y_{t-1},
%
% opt: 1 x 1, logical, = 1 if you want automatic bandwidth selection for
% HAR
% har: 1 x 1, logical, = 1 if you want to use the horizon of the LP for the
% coviarance matrix, 0 if you wanna use only heteroskedasticity consistent
% SE as argued by Plagborg-Moeller and Montiel Olea (ECTA, 2022)
% trend: 1 x 1, logical, = 1 if you want to add a time trend to the
% regressions
% order: vector of integers, column number of the variables used as
% controls in all equations. 


% first, check that the number of observations is correct

[T_y,n_y] = size(y);
[T_z,n_z] = size(z);
if (T_z ~= T_y); error('y and z have different number of observations'); end;


if length(varargin) == 1
    try
        level   = varargin{1}.level;
    catch
        level = 0;
    end
    
    try
        opt     = varargin{1}.opt;
    catch
        opt = 0;
    end
    
    try
        har     = varargin{1}.har;
    catch
        har = 0;
    end
    
    try
       trend   = varargin{1}.trend;
    catch
       trend = 0;
    end
     
    try
        order = varargin{1}.order;
    catch

    if isempty(x)
        order = 0;
    else
        order = 1:(n_y);
    end

     end
    

    try
        transformation = varargin{1}.transformation;
    catch
        transformation =1;
    end
   

    try
        quantiles_size = varargin{1}.quantiles_size;
    catch
        quantiles_size = [0.2, 0.8];
    end
   
    try
       only_non_zeros = varargin{1}.only_non_zeros;
    catch
        only_non_zeros = 0;
    end

else %defaults
    level = 0;
    opt   = 0;
    har   = 0;
    trend = 0;
    only_non_zeros = 0;
    transformation = 1; %1 is leave the default (linear), 1 lets you select if you want order 2 or 3 polynomial.
    quantiles_size = [0.2, 0.8];
    if isempty(x)
        order = 0;
    else
        order = 1:(n_y);
    end
end
%% Build matrix of non-linear transformations:

if strcmp(non_linearity,'sign')
    if transformation == 1
        z_vec = [z,abs(z)];
    elseif transformation == 2
        z_vec = [z,z.^2];
    elseif transformation == 3
        z_vec = [z,abs(z).*z.^2];
    end

elseif strcmp(non_linearity,'size')
    if transformation == 1
         % get quantiles for the linear case
         z_vec_quantiles = quantile(z(z~=0),quantiles_size);
         bound = max(abs(z_vec_quantiles));
        fun_size = @(x) ((x>bound).*(x-bound)+(x<-bound).*(x+bound));
        z_vec = [z,fun_size(z)];
    elseif transformation == 2
        z_vec = [z,z.*abs(z)];
    elseif transformation == 3
        z_vec = [z,z.^3];
    end


elseif strcmp(non_linearity,'linear')
    z_vec = z;
else
 error('Select size, sign or linear');
end

[T_z,n_z] = size(z_vec);

%% Once built, create the matrices of data to run the local projections.

if order == 0
if trend ==1 %add a trend
    x_mat = [ones(T_y,1), cumsum(ones(T_y,1))/T_y,];
elseif trend ==2 
    quadratic_trend = (1:T_y).^2/(T_y^2);
    x_mat = [ones(T_y,1), cumsum(ones(T_y,1))/T_y, quadratic_trend'];
else  % do not add trend
    x_mat = [ones(T_y,1)];
end

else

if trend ==1 %add a linear trend
    x_mat = [ones(T_y,1), cumsum(ones(T_y,1))/T_y, lagmatrix(x(:,order),1:lags_controls)];
elseif trend ==2 
    quadratic_trend = (1:T_y).^2/(T_y^2);
    x_mat = [ones(T_y,1), cumsum(ones(T_y,1))/T_y, quadratic_trend', lagmatrix(x(:,order),1:lags_controls)];
elseif trend ==3 
    quadratic_trend = (1:T_y).^2/(T_y^2);
    cubic_trend = (1:T_y).^3/(T_y^3);
    x_mat = [ones(T_y,1), cumsum(ones(T_y,1))/T_y, quadratic_trend', cubic_trend', lagmatrix(x(:,order),1:lags_controls)];
elseif trend == 4 
    quadratic_trend = (1:T_y).^2/(T_y^2);
    cubic_trend = (1:T_y).^3/(T_y^3);
    quartic_trend = (1:T_y).^4/(T_y^4);
    x_mat = [ones(T_y,1), cumsum(ones(T_y,1))/T_y, quadratic_trend', cubic_trend',quartic_trend', lagmatrix(x(:,order),1:lags_controls)];
else  % do not add trend
    x_mat = [ones(T_y,1), lagmatrix(x(:,order),1:lags_controls)];
end
end

x_use = [z_vec(lags_controls+1:end,:),x_mat(lags_controls+1:end,:)];
y_use = y(lags_controls+1:end,:);


n_x_use = size(x_use,2);
irf_coefs = zeros(n_y*n_z,horizon);
irf_coefs_v2 = zeros(n_y*n_z,horizon);
irf_se = zeros(n_y*n_z,horizon);
irf_t_stats =zeros(n_y*n_z,horizon);
irf_var_mat = zeros(n_z,n_z,horizon,n_y);
coefs_all = zeros(size(x_use,2)+lags_controls,horizon,n_y);
se_all = zeros(size(x_use,2)+lags_controls,horizon,n_y);
t_stats_all = zeros(size(x_use,2)+lags_controls,horizon,n_y);


for j=1:n_y %for each output variable
    betas = zeros(n_x_use+lags_controls,horizon);
    se = zeros(n_x_use+lags_controls,horizon);
    t_stats = zeros(n_x_use+lags_controls,horizon);

    if sum(order==j)==0 && only_non_zeros ~=1 %if the variable is not used as control for all, add lags of the same variable.
        y_lags_to_use = lagmatrix(y(:,j),1:lags_controls);
        x_use_2 = [x_use y_lags_to_use(lags_controls+1:end,:)];
    else
        x_use_2 = x_use; %do not add lags for the same variable, because it is already included in x.
    end 

    for i=1:horizon % for each horizon
         if level==1
            yy=y_use(i:end,j); %use the data in levels.
         else
            yy=y(lags_controls+i:end,j)-y(lags_controls:end-i,j); %use the long difference
         end

    if only_non_zeros ~= 1
         if har ==1
            results=nwest(yy(1:end-horizon+i), x_use_2(1:end-horizon+1,:),i);
         else
             results=nwest(yy(1:end-horizon+i), x_use_2(1:end-horizon+1,:),0);
         end
    else
        x_use_3 = x_use_2(1:end-i+1,:);
        index_use = (x_use_3(:,1) ~=0);
        results=nwest(yy(index_use),x_use_3(index_use,:),0);

    end
if sum(order==j)==0 %if the variable is not used as control for all
         betas(:,i)=results.beta;
         if opt==0
            se(:,i)=results.se;
         else
            [~, hacse, ~]=hac_alt(x_use_2(1:end-i+1,:), yy, 'intercept', false, 'smallT', false, 'display', 'off'); 
            se(:,i)=hacse';
         end
         t_stats(:,i) = betas(:,i)./se(:,i);
else
    betas(1:end-lags_controls,i)=results.beta;
    se(1:end-lags_controls,i)=results.se;
    t_stats(1:end-lags_controls,i) = results.beta./(results.se);
end
results_2.resid = results.resid;
results_collect(i,j) = results_2;
         irf_var_mat(:,:,i,j) = results.V(1:n_z,1:n_z);
    end

% store results for each variable
irf_coefs((j-1)*n_z+1:n_z*j,:) = betas(1:n_z,:);
irf_se((j-1)*n_z+1:n_z*j,:) = se(1:n_z,:);
irf_t_stats((j-1)*n_z+1:n_z*j,:) = t_stats(1:n_z,:);
coefs_all(:,:,j) = betas;
se_all(:,:,j) = se;
t_stats_all(:,:,j) = t_stats;
end
extra_stuff.coefs_all = coefs_all;
extra_stuff.se_all = se_all;
extra_stuff.t_stats_all = t_stats_all;
extra_stuff.results_collect = results_collect;

end


