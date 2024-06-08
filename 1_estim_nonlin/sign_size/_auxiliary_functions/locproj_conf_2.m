function obj = locproj_conf_2(varargin)

    switch length(varargin)
        case 3
            obj     = varargin{1};
            H       = varargin{2};
            nlag    = varargin{3};
            lambda  = 0.0;
        case 4
            obj     = varargin{1};
            H       = varargin{2};
            nlag    = varargin{3};
            lambda  = varargin{4};
            
        otherwise
            error('wrong number of input arguments')
    end 

    % POINT ESTIMATE
    XXP   = ( obj.X'*obj.X + lambda * obj.P );
    theta = XXP \ ( obj.X'*obj.Y );

    % COMPUTE NW ESTIMATOR
    % BREAD
    bread = XXP^-1;
    n_x = obj.n_x;
    K = obj.K;
    HR = obj.HR;
    % MEAT
    %nlag    = H;
    T       = obj.T;
    npar    = length(obj.theta);
    weights = [ 0.5 (nlag+1-(1:nlag))/(nlag+1) ];    
    idx     = obj.idx;
    U       = obj.Y - obj.X * obj.theta;    
    X       = obj.X ;
    V       = zeros( npar , npar );
    X_by_t  = zeros(obj.HR,npar,T-obj.HR-1); 
    U_by_t  = zeros(obj.HR,T-obj.HR-1); 
    
    % reorder these matrices by t to avoid having to find indexes multiple
    % times
    for t = 1:(T-HR-1) 
            X_by_t(:,:,t) = X( idx(:,1)==t , : );
            U_by_t(:,t) = U( idx(:,1)==t );
    end
% paralelize for lags
% see if you can write the sum as a tensor product.
V_collector = zeros(npar,npar,nlag+1);

    for l = 0:nlag
        tic
        GplusGprime = zeros( npar , npar );
        for t = (l+1):(T-HR-1) 
            S1 = X_by_t(:,:,t)' * U_by_t(:,t);
            S2 = X_by_t(:,:,t-l)' * U_by_t(:,t-l);
            GplusGprime = GplusGprime + S1 * S2' + S2 * S1';
        end
        V_collector(:,:,l+1) = GplusGprime;
        toc
    end

for l = 0:nlag
V = V + weights(l+1) * V_collector(:,:,l+1);
end

    meat = V;

    
%     for l = 0:nlag
%         GplusGprime = zeros( npar , npar );
%         for t = (l+1):(T-obj.HR-1) 
%             S1 = X( idx(:,1)==t , : )' * U( idx(:,1)==t );
%             S2 = X( idx(:,1)==(t-l) , : )' * U( idx(:,1)==(t-l) );
%             GplusGprime = GplusGprime + S1 * S2' + S2 * S1';
%         end
%         V = V + weights(l+1) * GplusGprime;
%     end
%     meat = V;

    VC = bread * meat * bread;
    
    conf = nan( obj.H_max+1 , n_x*2);
    ster = nan( obj.H_max+1 , n_x);
    n_hor = obj.H_max+1-obj.H_min; %horizon.

    for i_x = 1:n_x
    if strcmp(obj.type,'reg')==1
        ster(1+obj.H_min:end,i_x) = sqrt( diag( VC( n_hor*(i_x-1)+1:i_x*n_hor , n_hor*(i_x-1)+1:i_x*n_hor ) ) );                
        conf(1+obj.H_min:end,(i_x-1)*n_x+1) = theta((i_x-1)*obj.K+1:i_x*obj.K)*obj.delta + ster(:,i_x)*obj.delta*norminv(0.05);
        conf(1+obj.H_min:end,(i_x-1)*n_x+2) = theta((i_x-1)*obj.K+1:i_x*obj.K)*obj.delta + ster(:,i_x)*obj.delta*norminv(0.95);
    else
        ster(1+obj.H_min:end,i_x) = sqrt( diag( obj.B*VC( (i_x-1)*obj.K+1:i_x*obj.K , (i_x-1)*obj.K+1:i_x*obj.K )*obj.B' ) );                
        conf(1+obj.H_min:end,(i_x-1)*n_x+1) = obj.B*theta((i_x-1)*obj.K+1:i_x*obj.K)*obj.delta + ster(:,i_x)*obj.delta*norminv(0.05);
        conf(1+obj.H_min:end,(i_x-1)*n_x+2) = obj.B*theta((i_x-1)*obj.K+1:i_x*obj.K)*obj.delta + ster(:,i_x)*obj.delta*norminv(0.95);
    end
    end
    obj.ster = ster;
    obj.conf = conf;
    obj.var_cov_mat = VC(1:n_x*n_hor,1:n_x*n_hor);
end
