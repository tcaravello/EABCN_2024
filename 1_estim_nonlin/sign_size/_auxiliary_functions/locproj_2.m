function obj = locproj_2(varargin)
% 
% locproj creates a struct object. Once you create the object, you can see it in the workspace. 
% You can see all the properties the object contains by clicking on it. 
% 
% obj = locproj(y,x,w,H_min,H_max,type)
% obj = locproj(y,x,w,H_min,H_max,type,r,lambda)
% 
% Description
% obj = locproj(y,x,w,H_min,H_max,type) returns a local projection struct with the 
% instrumental impulse response coefficients for a normal local projection.
% 
% obj = locproj(y,x,w,H_min,H_max,type,r,lambda) returns a local projection struct with
% instrumental impulse response coefficients for a smooth local projection.
% 
% Input Arguments
% y - Response (LHS) vector
% x - Endogenous impulse vector
% w - Optional matrix of all other RHS vectors, including the constant, instrumental variables, and lagged variables (if applicable)
% H_min - First period the shock hits (0 or 1)
% H_max - Number of periods the local projection spans
% type - tring value that dictates whether the local projection ran is regular or smooth ('reg' or 'smooth')
% r - Used in smooth local projection: order of the limit polynomial  (for example, r=2 implies the impulse response is shrunk towards a line)
% lambda - Used in smooth local projection: numeric value for the shrinking parameter lambda
% 
% Output Arguments
% obj.T - Number of observations
% obj.H_min - First period the shock hits (0 or 1)
% obj.H_max - Number of periods the local projection spans
% obj.HR - Number of periods of impulse responses (H_max + 1 - H_min)
% obj.K - Number of columns of the basis spline
% obj.B - Basis spline used in smooth local projection
% obj.P - Penalty matrix
% obj.lambda - lambda used above as an input
% obj.type - Type of local project, smooth or reg
% obj.delta - numeric value containing the residual standard error from the first stage OLS regression of the endogenous impulse variable on other endogenous and exogenous RHS variables
% obj.idx - Observation horizon matrix
% obj.Y - Shifted response (LHS) vector ran in local projection simulatenously over all horizons
% obj.X - Shifted RHS matrix ran in local projection simultaenously over all horizons (includes all RHS vectors including the constant)
% obj.theta - Theta matrix containing coefficient estimates for all variables, for each value of labmda
% obj.IR - matrix containing estimations of the instrumented impulse coefficient for each value of lambda 
% obj.n_x = n_x; %number of variables in the IRF 


    switch length(varargin)
        case 6
            y     = varargin{1};
            x     = varargin{2};
            w     = varargin{3};
            H_min = varargin{4};
            H_max = varargin{5};
            type  = varargin{6};

        case 8
            y       = varargin{1};
            x       = varargin{2};
            w       = varargin{3};
            H_min   = varargin{4};
            H_max   = varargin{5};
            type    = varargin{6};            
            r       = varargin{7};
            lambda  = varargin{8};
            
        otherwise
            error('wrong number of input arguments')
    end    

    obj = struct();
    
%     if isempty(w)
%         delta = std(x); %scale everything by the coefficient of this
%     else
%         delta = std( x-w*inv(w'*w)*w'*x ); %scale everything by the coefficient of this, residualized.
%     end

delta =1;
    
    % isreg == 1 if it is a normal regression, 0 otherwise.
    isreg = strcmp('reg',type);
    
    T  = length(y); %number of obs
    HR = H_max + 1 - H_min; %adjustment for building matrices
    n_x = size(x,2);
    % construct the B-spline basis functions. One knot per horizon
    if ~isreg
        B = bspline( (H_min:H_max)' , H_min , H_max+1 , H_max+1-H_min , 3 );
        K = size( B , 2 ); 
    else
        K = HR;
    end

    % building up the regression representation of the local projection
    idx = nan( (H_max+1)*T , 2 );
    Y   = nan( (H_max+1)*T , 1 );
    Xb  = zeros( (H_max+1)*T , n_x*K );
    Xc  = zeros( (H_max+1)*T , HR , size(w,2)+1 );
    
    w = [ ones(T,1) w ];
    ;
    for t = 1:T-H_min
        
        idx_beg = (t-1)*HR + 1;
        idx_end = t*HR;

        idx( idx_beg:idx_end , 1 ) = t;
        idx( idx_beg:idx_end , 2 ) = H_min:H_max;
        
        % y
        y_range = (t+H_min) : min((t+H_max),T)';
        Y( idx_beg:idx_end ) = [ y( y_range ) ; nan(HR-length(y_range),1) ];

        % regressors of interest. HERE YOU SHOULD MODIFY
        if isreg
            Xb( idx_beg:idx_end , : ) = kron(x(t,:),eye(HR));
        else
            Xb( idx_beg:idx_end , : ) =kron(x(t,:),B);
        end

        % controls. As we can see, those are not being shrinked. In these
        % loops, i is which control is begin added.
        for i = 1:size(w,2)
            Xc( idx_beg:idx_end , : , i ) = eye(HR)*w(t,i);
        end
        
    end
    
    X = Xb;
    % THIS PART IS BASICALLY STACKING ALL CONTROLS TOGHETHER. CONTROLS ARE
    % NOT BEING SHRINKED.
    for i = 1:size(w,2)
        X = [X Xc(:,:,i)];
    end
    
    select = isfinite(Y);  %this part creates an index that equals 1 if the variable has a value, zero otherwise
    idx = idx(select,:);
    Y   = Y(select); %this part creates an index that equals 1 if the variable has a value, zero otherwise
    X   = X(select,:);%this part creates an index that equals 1 if the variable has a value, zero otherwise
    X   = sparse(X);%this part creates an index that equals 1 if the variable has a value, zero otherwise

    % estimation
    IR  = zeros(H_max+1,n_x);
    
    if isreg

        theta     = ( X'*X )\( X'*Y );
        for n = 1:n_x
        IR((H_min+1):end,n) = theta((n-1)*K+1:n*K);
        end
    else
        % this part creates the penalization matrix. 
        P = zeros( size(X,2) );
        % build the D matrix
        D = eye(K);
        for k = 1:r 
            D = diff(D);
        end
        % penalty matrix. just put non-zeros in the IRF part. The controls
        % are not being shrinked.
        for n = 1:n_x;
            P((n-1)*K+1:n*K,(n-1)*K+1:n*K) = D' * D;
        end
        % estimated coefficients
        theta = ( X'*X + lambda*P )\( X'*Y );
        % IR is the rotation of the coefficients.
        %IR((1+H_min):end) = B * theta(1:n_x*K) * delta; %the delta is just rescaling.
        for n = 1:n_x
        IR((H_min+1):end,n) = B* theta((n-1)*K+1:n*K);
        end
    end
    
    % pack everything up
    obj.T     = T;
    obj.H_min = H_min;
    obj.H_max = H_max;
    obj.HR    = HR;
    obj.K     = K;
    
    if isreg
        obj.B      = 0;
        obj.P      = zeros( size(X,2) );
        obj.lambda = 0;
    else
        obj.B      = B;
        obj.P      = P;
        obj.lambda = lambda;
    end
    
    obj.type  = type;
    obj.delta = delta;
    obj.idx   = idx;
    obj.Y     = Y;
    obj.X     = X;
    obj.theta = theta;
    obj.IR    = IR;
    obj.delta = delta;
    obj.n_x = n_x; %number of variables in the IRF
    
    % debug stuff
    % Bs is for display / debugging pourposes only
    % Bs = bspline( (H_min:0.1:H)' , H_min , H+1 , H+1-H_min , 3 );
    % obj.Bs = Bs;
end

function B = bspline(x, xl, xr, ndx, bdeg)
    dx = (xr - xl) / ndx;
    t = xl + dx * [-bdeg:ndx-1];
    T = (0 * x + 1) * t;
    X = x * (0 * t + 1);
    P = (X - T) / dx;
    B = (T <= X) & (X < (T + dx));
    r = [2:length(t) 1];
    for k = 1:bdeg
        B = (P .* B + (k + 1 - P) .* B(:, r)) / k;
    end
end
