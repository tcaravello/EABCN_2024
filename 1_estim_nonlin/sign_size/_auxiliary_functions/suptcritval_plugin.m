function [suptplugin]=suptcritval_plugin(confidence,Sigma,I)
            %This function reports the sup-t critical value, based on
            %multivariate normal draws.
            %INPUT
            %------
            %a)confidence: confidence level                 (1 x 1)
            %b)     Sigma: cov matrix of theta              (k x k)
            %              (this matrix can be singular)
            %c)         I: number of draws to compute the 
            %              supt critical value              (1 x 1)   
            %------
            %OUTPUT
            %a)suptplugin: supt critical value              (1 x 1)
            %              (we use plug-in, as we the input
            %              Sigma in this function is a consistent
            %              estimator of the true covariance matrix)
            %------
            d           = diag(Sigma).^(.5);     % Extracts the diagonal 
                                                 % of Vartheta
                                               
            d_nonzero   = (d>eps);               % Elements with nonzero 
                                                 % variance
                                                 % This line implicitly 
                                                 % drops out all the zero
                                                 % variance components of
                                                 % Signa
                                               
            Corrmat     = bsxfun(@ldivide, ...
                                 d(d_nonzero), ...
                                 bsxfun(@rdivide,...
                                 Sigma(d_nonzero,d_nonzero),...
                                 d(d_nonzero)')); %Matrix of Correlations
                             
            % Compute square root using eigendecomposition
            [Corrmat_V,Corrmat_D] ...
                         = eig(Corrmat);
            
            Corrmat_sqrt = sqrt(Corrmat_D)*Corrmat_V';
            
            t            = abs(randn(I,size(Corrmat_sqrt,1))*Corrmat_sqrt); %limit of tstats
            
            suptplugin   = quantile(max(t,[],2),confidence,1);
            
        end