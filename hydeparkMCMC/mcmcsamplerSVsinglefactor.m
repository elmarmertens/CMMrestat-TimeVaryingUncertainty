function [feSVdraws, feDraws, etaSVdraws, hbar, bbeta] = mcmcsamplerSVsinglefactor(eta, etaNanny, T, horizons, MCMCdraws, Nfedraws, doZeroSlopes, doDiffuseSlopes, rndStream, showProgress)

% NOTE: SVdraws are *variance* draws (take sqrt to convert into Stoch*Vol*)

eta         = eta(1:T,:);
etaNanny    = etaNanny(1:T,:);

Nsurvey     = size(eta,2);
% [KSC, KSCt] = getKSCvalues(T, Nsurvey);
[KSC, KSCt] = getKSC10values(T, Nsurvey); % using 10-point grid as in Omori et al (2007, JoE)

lambdavarDof   = 2 + 1;
lambdavarT     = 0.2^2 * (lambdavarDof - 2);

beta0  = ones(Nsurvey, 1);
beta0V = .5 * ones(Nsurvey,1);
beta0V(1) = 0; % normalization beta(1) = 1;

Nshockslopes  = Nsurvey * (Nsurvey - 1) / 2;
E0shockslopes = zeros(Nshockslopes, 1);
V0shockslopes = eye(Nshockslopes);


%% set up SV



Eh0  = repmat(log(.5) * 2, Nsurvey,1);
Vh0  = repmat(10,Nsurvey,1);

bbeta     = NaN(Nsurvey, MCMCdraws);
hbar      = NaN(Nsurvey, MCMCdraws);
hh        = NaN(Nsurvey, T, MCMCdraws);
lambda    = NaN(T, MCMCdraws);
lambdavar = NaN(1,MCMCdraws);


SHOCKSLOPES = NaN(Nshockslopes,MCMCdraws);

SVinno = NaN(Nsurvey, T);

%% draw initial values
lambdavarPREV = igammadraw(lambdavarT, lambdavarDof, 1, rndStream);

betaPREV    = beta0 + sqrt(beta0V) .* randn(rndStream, Nsurvey, 1); 
betaPREV(1) = 1; % to be safe

lambda0    = 0;
lambdaPREV = sqrt(lambdavarPREV) * randn(rndStream, 1, T);
lambdaPREV = lambda0 + cumsum(lambdaPREV);
hbarPREV   = Eh0 + sqrt(Vh0) .* randn(rndStream, Nsurvey, 1);

hPREV      = bsxfun(@plus, betaPREV * lambdaPREV, hbarPREV);

if doZeroSlopes
    shockslopesPREV = zeros(Nshockslopes,1);
else
    shockslopesPREV = E0shockslopes + chol(V0shockslopes)' * randn(rndStream, Nshockslopes, 1);
end

% diffuse priors for hVar:
% hvarDof = 0;
% hvarT   = zeros(Nsurvey,1);

% burnin
if showProgress
    progressbar(0)
end
burnin = MCMCdraws;
for n = -burnin : MCMCdraws
    
    %% fill in missing values
    B      = eye(Nsurvey);
    offset = 0;
    for j = 2 : Nsurvey
        ndx        = offset + (1:j-1);
        B(j,1:j-1) = shockslopesPREV(ndx);
        offset     = ndx(end);
    end
    
    if any(etaNanny(:))
        %% draw missing values
        for t = find(any(etaNanny,2))'
            thisNdx = ~etaNanny(t,:);
            Nmissing = sum(~thisNdx);
            
            sigmaChol = B * diag(exp(hPREV(:,t) * .5));
            
            if Nmissing < Nsurvey
                sigma = sigmaChol * sigmaChol';
                
                etaData = eta(t,thisNdx);
                varData = sigma(thisNdx, thisNdx);
                covData = sigma(~thisNdx, thisNdx);
                beta    = covData / varData;
                
                varMissing  = sigma(~thisNdx, ~thisNdx);
                residVar    = varMissing - beta * covData';
                missingData = beta * etaData' + chol(residVar)' * randn(rndStream, Nmissing, 1);
                eta(t, ~thisNdx) = missingData;
            else
                missingData = randn(rndStream, 1, Nsurvey) * sigmaChol';
            end
            eta(t, ~thisNdx) = missingData;
        end
    end
    
    %% remainder of MCMC
    
    if doZeroSlopes
        SVinno = eta';
    else
        % run regressions to create Choleski of eta
        SV = exp(hPREV * 0.5);
        SVinno(1,:) = eta(:,1)';
        offset = 0;
        for m = 2 : Nsurvey
            ndx = offset + (1:m-1);
            rhs = bsxfun(@rdivide, SVinno(1:m-1,:), SV(m,:))';
            lhs = eta(:,m) ./ SV(m,:)';
            if doDiffuseSlopes
                [shockslopesPREV(ndx),resid] = bayesRegressionSlopesGibbsDrawDiffuse(lhs, rhs, 1, 1, rndStream);
            else
                [shockslopesPREV(ndx),resid] = bayesRegressionSlopesGibbsDraw(lhs, rhs, 1, E0shockslopes(ndx), V0shockslopes(ndx,ndx), 1, rndStream);
            end
            SVinno(m,:) = resid' .* SV(m,:);
            offset = ndx(end);
        end
    end
    
    % KSC draw
    logy2           = log(SVinno.^2 + 0.0001); % the original KSC practice was to add 0.001; we are using however 10-point grid as in Omori et al. (2007, JoE)
    
    
    [hPREV, lambdaPREV, lambdaResid, hbarPREV, kai2states] = StochVolSingleFactor(logy2, ...
        hPREV, betaPREV, lambdavarPREV, Eh0, Vh0, KSC, KSCt, Nsurvey, T, rndStream);

    
    % draw lambdavar
    lambdavarPREV =  bayesVCVgibbsDraw1(lambdavarT, lambdavarDof, lambdaResid', rndStream, true);
    
    % draw beta
    % - adjust obs
    obs = bsxfun(@minus, logy2, hbarPREV);
    obs = (obs - KSC.mean(kai2states)) ./ KSC.vol(kai2states);
    lambdatilde = bsxfun(@rdivide, lambdaPREV, KSC.vol(kai2states));
    for nn = 2 : Nsurvey
        betaPREV(nn) = bayesRegressionSlopesGibbsDraw(obs(nn,:)', lambdatilde(nn,:)', ...
            1, beta0(nn), beta0V(nn), 1, rndStream);
    end
    
    %% store draws post burnin
    if n > 0
        hh(:,:,n)             = hPREV;
        hbar(:,n)             = hbarPREV;
        lambda(:,n)           = lambdaPREV;
        lambdavar(n)          = lambdavarPREV;   
        bbeta(:,n)            = betaPREV;
        SHOCKSLOPES(:,n)      = shockslopesPREV;
    end
    
    if showProgress
        progressbar((n + burnin + 1) / (burnin +  MCMCdraws + 1))
    end
end

%% compute forecast error SV
% allocate memory
feSVdraws  = zeros(T,Nsurvey,MCMCdraws);
etaSVdraws = NaN(T,Nsurvey,MCMCdraws);
SVhorizons = horizons + 1; % We use t-1 data to project the SV's for t, t+1 etc

% progressbar(0)
for dd = 1 : MCMCdraws
    
        
    % construct B
    B      = eye(Nsurvey);
    offset = 0;
    for j = 2 : Nsurvey
        ndx        = offset + (1:j-1);
        B(j,1:j-1) = SHOCKSLOPES(ndx,dd);
        offset     = ndx(end);
    end
    
    % construct forecasts of hbar up to maxhorizon steps ahead
    VarianceForecast  = NaN(T,Nsurvey,Nsurvey);
    for j = 1 : Nsurvey
        k = SVhorizons(j);
        condmean = hh(:,:,dd);
        condvar  = bbeta(:,dd).^2 * k * lambdavar(dd);
        for n = 1 : Nsurvey
            VarianceForecast(:,n,j) = exp(condmean(n,:) + .5 * condvar(n));
        end
    end
    
    %% more generic code:
    % 1): forecast Omega(t+j)
    Omega = NaN(Nsurvey,Nsurvey,Nsurvey,T);
    for t = 1 : T
        for j = 1 : Nsurvey
            Omega(:,:,j,t) = B * diag(VarianceForecast(t,:,j)) * B';
        end
    end
    
    % feSV
    for t = 1 : T
        for n = 1 : Nsurvey
            %             feSVdraws(t,n,dd) = 0; % redundant since array already initialized to zero
            for j = 1 : n
                feSVdraws(t,n,dd) = feSVdraws(t,n,dd) + Omega(j,j,n - j + 1,t);
            end
        end
    end
    
    % one step ahead forecasts of etaSV
    j = 1;
    for t = 1 : T
        etaSVdraws(t,:,dd) = diag(Omega(:,:,j,t));
    end
    
    % progressbar(dd/MCMCdraws)
end

%% simulate forecast errors
feDraws = zeros(Nsurvey,Nfedraws,MCMCdraws);

% model: eta = B * Vol * z
for dd = 1 : MCMCdraws
    
    % construct B
    B      = eye(Nsurvey);
    offset = 0;
    for j = 2 : Nsurvey
        ndx        = offset + (1:j-1);
        B(j,1:j-1) = SHOCKSLOPES(ndx,dd);
        offset     = ndx(end);
    end
    
    % generate SV (note the second dimension are the simulated forecasts
    lambdasim = lambda(end,dd) + cumsum(sqrt(lambdavar(dd)) * randn(rndStream,Nsurvey,Nfedraws));
    hsim      = NaN(Nsurvey, Nsurvey, Nfedraws);
    for m = 1 : Nfedraws
        hsim(:,:,m) = bbeta(:,dd) * lambdasim(:,m)' + hbar(:,dd);
    end
    
    % generate z and scale by SV
    z = randn(rndStream,Nsurvey,Nsurvey,Nfedraws);
    z = exp(hsim * .5) .* z;
    z   = reshape(z, Nsurvey, Nsurvey * Nfedraws);
    eta = B * z;
    eta = reshape(eta, Nsurvey, Nsurvey, Nfedraws);
    
    
    % note: 1st col treated like "time when eta happens", 2nd column is the
    % horizon of eta
    % i.e. eta(t,h) is eta(t+h|t)
    
    % construct FE from eta
    for h = 0 : horizons(end)
        feDraws(h+1,:,dd) = 0; % this line should be redundant
        for j = 0 : h
            feDraws(h+1,:,dd) = feDraws(h+1,:,dd) + permute(eta(1+h-j,j+1,:), [1 3 2]);
        end
    end
    
end

feDraws = reshape(feDraws, Nsurvey, Nfedraws * MCMCdraws);
