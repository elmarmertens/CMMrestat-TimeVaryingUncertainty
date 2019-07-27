function [feSVdraws, feDraws, etaSVdraws] = mcmcsamplerConst(eta, etaNanny, windowSize, T, horizons, MCMCdraws, Nfedraws, doZeroSlopes, doDiffuseSlopes, rndStream)

eta         = eta(T - windowSize + 1 : T,:);
etaNanny    = etaNanny(T - windowSize + 1 : T,:);
Nsurvey     = size(eta,2);



doDiffusePrior = true;
%#ok<*UNRCH>

zeroNsurvey = zeros(Nsurvey); 

SigmaDof = Nsurvey + 2;
SigmaT   = eye(Nsurvey) * (SigmaDof - Nsurvey - 1); % if doDiffusePrior is true, this prior only needed to init missing values

SIGMA = NaN(Nsurvey,Nsurvey,MCMCdraws);

sigma = iwishdraw(SigmaT, SigmaDof, 1, [], rndStream);


if doZeroSlopes
    error('doZeroslopes not yet implemented');
end

if ~doDiffuseSlopes
    error('doDiffuseSlopes == false: not yet implemented');
end


for n = -MCMCdraws : MCMCdraws
    
    
    if any(etaNanny(:))
        
        %% draw missing values
        for t = find(any(etaNanny,2))'
            thisNdx = ~etaNanny(t,:);
            Nmissing = sum(~thisNdx);
            
            if Nmissing < Nsurvey
                etaData = eta(t,thisNdx);
                varData = sigma(thisNdx, thisNdx);
                covData = sigma(~thisNdx, thisNdx);
                beta    = covData / varData;
                
                varMissing  = sigma(~thisNdx, ~thisNdx);
                residVar    = varMissing - beta * varData * beta'; % bit clumys but maybe more stable numerically than beta * covData
                missingData = beta * etaData' + chol(residVar)' * randn(rndStream, Nmissing, 1);
                eta(t, ~thisNdx) = missingData;
            else
                missingData = randn(rndStream, 1, Nsurvey) * chol(sigma);
            end
            eta(t, ~thisNdx) = missingData;
        end
    end
    
    count = 0;
    OK = false;
    while ~OK && count < 1e2
        if doDiffusePrior
            sigma = bayesVCVgibbsDraw1(zeroNsurvey, 0, eta, rndStream, false);
        else
            sigma = bayesVCVgibbsDraw1(SigmaT, SigmaDof, eta, rndStream, false); 
        end
        OK = abs(det(sigma)) > 1e-6;
        count = count + 1;
    end
    
    %% store draws post burnin
    if n > 0
        SIGMA(:,:,n)     = sigma;
    end
    
    % progressbar((n + MCMCdraws + 1) / (2 * MCMCdraws + 1))
end


%% compute forecast error SV
% allocate memory
feSVdraws  = NaN(T,Nsurvey,MCMCdraws);
etaSVdraws = NaN(T,Nsurvey,MCMCdraws);

% progressbar(0)
for dd = 1 : MCMCdraws
    
    etaSVdraws(:,:,dd) = repmat(diag(SIGMA(:,:,dd))', [T , 1]);
    feSVdraws(:,:,dd)  = cumsum(etaSVdraws(:,:,dd),2);
end

%% simulate forecast errors
feDraws = zeros(Nsurvey,Nfedraws,MCMCdraws);

% model: eta = cholSigma* z
for n = 1 : MCMCdraws
    
    B = chol(SIGMA(:,:,n))';
    
    %     % generate z (old code, slower)
    %     z = randn(rndStream,Nsurvey,Nsurvey,Nfedraws);
    %     eta = NaN(Nsurvey,Nsurvey,Nfedraws);
    %     for m = 1 : Nfedraws
    %         eta(:,:,m) = B * z(:,:,m);
    %     end
    
    % generate z
    z = randn(rndStream,Nsurvey,Nsurvey * Nfedraws);    
    eta = B * z;
    eta = reshape(eta, Nsurvey, Nsurvey, Nfedraws);
    % note: 2nd col treated like "time when eta happens", 3rd column is the
    % horizon of eta
    % i.e. eta(t,h) is eta(t+h|t)
    
    % construct FE from eta
    for h = 0 : horizons(end)
        feDraws(h+1,:,n) = 0;
        for j = 0 : h
            feDraws(h+1,:,n) = feDraws(h+1,:,n) + permute(eta(1+h-j,j+1,:), [1 3 2]);
        end
    end
    
end

feDraws = reshape(feDraws, Nsurvey, Nfedraws * MCMCdraws);


