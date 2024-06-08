%% Run the internal instruments VAR

% run reduced-form VAR 
if cumProxy
    dataext = [cumsum(proxy) data];
else
    dataext = [proxy data];
end

varEst = varxest(dataext,dataExo,p);
nvar = nvar+1;

% identification using Cholesky decomposition of instrument-augmented VAR

nexo = size(dataExo,2);

A0 = chol(varEst.Sigma)';
b1      = A0(:,1);

% compute IRFs to shock
IRFs_proxy = varirfsingle(varEst.B(:,1+nexo:end),b1,p,horizon);
if strcmp(shockType,'custom')
    IRFs_proxy = IRFs_proxy./IRFs_proxy(1,2)*shockSize;
end

% compute FEVDs to shock
if compFEVDs
    FEVDs_proxy = varxfevdsingle(varEst,b1,horizon+1); 
end

% compute the confidence bands using bootstrapping
bootIRFs = nan(horizon+1,nvar,nsim);
if compFEVDs
    bootFEVDs = nan(horizon+1,nvar,nsim);
end
T_est = varEst.T; % length of estimation sample

if strcmp(bootType,'mbb1block')
    % if identification sample is shorter that estimation sample, censor
    % unobserved values to zero

    BlockSize = round(5.03*T_est^0.25);
    nBlock = ceil(T_est/BlockSize);
    VARBlocks = zeros(BlockSize,nvar,T_est-BlockSize+1);
    for j = 1:T_est-BlockSize+1
        VARBlocks(:,:,j) = varEst.U(j:BlockSize+j-1,:);
    end

    % center the bootstrapped VAR errors
    VARcentering = zeros(BlockSize,nvar);
    for j = 1:BlockSize
        VARcentering(j,:) = mean(varEst.U(j:T_est-BlockSize+j,:),1);
    end
    VARcentering = repmat(VARcentering,[nBlock,1]);
    VARcentering = VARcentering(1:T_est,:);
    
end

j = 1;
while j <= nsim
    % generate artificial data
    
    if strcmp(bootType,'mbb1block')
        % Moving block bootstrap (Lundsford and Jentsch) using one block
        % type
        
        % draw bootstrapped residuals and proxies
        index = ceil((T_est - BlockSize + 1)*rand(nBlock,1));
        bootU = zeros(nBlock*BlockSize,nvar);
        for kk = 1:nBlock
            bootU(1+BlockSize*(kk-1):BlockSize*kk,:) = VARBlocks(:,:,index(kk,1));
        end
        bootU = bootU(1:T_est,:);
        
        %center the bootstrapped residuals and proxies
        bootU = bootU - VARcentering;

        % simulate VAR starting from initial values
        Xexo = varEst.Xexo;
        
        bootData = zeros(T_est+p,nvar); 
        bootData(1:p,:) = [proxy(1:p) data(1:p,:)]; %initial values of y, same for all j
        for i = p+1:T_est+p
            bootData(i,:)= varEst.B*[Xexo(i-p,:)'; vec(fliplr(bootData(i-p:i-1,:)'))] ...
                             + bootU(i-p,:)'; % bootstrap
        end
        
        % re-estimate the VAR
        bootvarEst = varxest(bootData,dataExo,p);

    end
        
    % structural impact matrix
    % 2SLS of U1 on U2 using proxy as an instrument (include a constant for the case that proxy is not demeaned)

    % first stage
    bootA0 = chol(bootvarEst.Sigma)';
    bootb1     = bootA0(:,1);

    % compute IRFs
    bootIRFs(:,:,j)  = varirfsingle(bootvarEst.B(:,1+nexo:end),bootb1,p,horizon);
    if strcmp(shockType,'custom')
        bootIRFs(:,:,j) = bootIRFs(:,:,j)./bootIRFs(1,2,j)*shockSize;
    end 
    if compFEVDs
        bootFEVDs(:,:,j) = varxfevdsingle(bootvarEst,bootb1,horizon+1); 
    end
    
    j = j+1;
end

% rescale bootstrapped IRFs (center around sample estimates) and get quantiles
IRFsmed_proxy = quantile(bootIRFs, 0.5, 3);

IRFsupper_proxy = quantile(bootIRFs, 1-alpha/2, 3)-IRFsmed_proxy+IRFs_proxy;  % rescaling does not affect the ordering, thus it is fine to take quantile first
IRFslower_proxy = quantile(bootIRFs, alpha/2, 3)-IRFsmed_proxy+IRFs_proxy;

IRFsupper2_proxy = quantile(bootIRFs, 1-0.32/2, 3)-IRFsmed_proxy+IRFs_proxy;  % rescaling does not affect the ordering, thus it is fine to take quantile first
IRFslower2_proxy = quantile(bootIRFs, 0.32/2, 3)-IRFsmed_proxy+IRFs_proxy;
