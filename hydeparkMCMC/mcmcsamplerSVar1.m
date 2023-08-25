function [feSVdraws, feDraws, etaSVdraws, hhVCV, RHO, HBAR] = mcmcsamplerSVar1(eta, etaNanny, T, horizons, MCMCdraws, Nfedraws, doZeroSlopes, doDiffuseSlopes, rndStream, showProgress)

% NOTE: SVdraws are *variance* draws (take sqrt to convert into Stoch*Vol*)


maxShake    = 1e3;

eta         = eta(1:T,:);
etaNanny    = etaNanny(1:T,:);

Nsurvey     = size(eta,2);
% [KSC, KSCt] = getKSCvalues(T, Nsurvey);
[KSC, KSCt] = getKSC10values(T, Nsurvey); % using 10-point grid as in Omori et al (2007, JoE)

hvarDof   = 8 + Nsurvey + 1;
hvarT     = 0.2^2 * (hvarDof - Nsurvey - 1) * eye(Nsurvey);
hvarTsqrt = sqrt(0.2^2 * (hvarDof - Nsurvey - 1)) * eye(Nsurvey);

Nshockslopes  = Nsurvey * (Nsurvey - 1) / 2;
E0shockslopes = zeros(Nshockslopes, 1);
V0shockslopes = eye(Nshockslopes);


%% set up SV



Eh0  = repmat(log(.5) * 2, Nsurvey,1);
Vh0  = repmat(10,Nsurvey,1);

hh   = NaN([Nsurvey T MCMCdraws]);

hhVCV = NaN(Nsurvey, Nsurvey, MCMCdraws);
% generate initial values
hvcv  = iwishdraw(hvarT, hvarDof, 1, hvarTsqrt, rndStream);
h0    = Eh0 + sqrt(Vh0) .* randn(rndStream, Nsurvey, 1);

rho0       = repmat(0.8, Nsurvey, 1);
rhoV0i     = (1 / 0.2^2) * eye(Nsurvey);
rhosqrtV0i = (1 / 0.2) * eye(Nsurvey);

RHO         = NaN(Nsurvey,MCMCdraws);
HBAR        = NaN(Nsurvey,MCMCdraws);
SHOCKSLOPES = NaN(Nshockslopes,MCMCdraws);

SVinno = NaN(Nsurvey, T);

%% draw initial values
hPREV = chol(hvcv)' * randn(rndStream, Nsurvey, T);
hPREV = cumsum(hPREV);
hPREV = bsxfun(@plus, hPREV, h0);

if doZeroSlopes
    shockslopesPREV = zeros(Nshockslopes,1);
else
    shockslopesPREV = E0shockslopes + chol(V0shockslopes)' * randn(rndStream, Nshockslopes, 1);
end

rhoPREV = repmat(10,Nsurvey, 1);
ndx = abs(rhoPREV) >= 1;
while any(ndx)
    rhodraw       = rho0 + transpose(rhosqrtV0i) \ randn(rndStream, Nsurvey, 1);
    rhoPREV(ndx)  = rhodraw(ndx);
    ndx           = abs(rhoPREV) >= 1;
end
clear rhodraw


% burnin
if showProgress
    progressbar(0)
end
burnin = MCMCdraws;
n = -burnin;
while n <= MCMCdraws
    
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
                residVar    = varMissing - beta * varData * beta'; % bit clumys but numerically more stable than beta * covData
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
    logy2           = log(SVinno.^2 + 0.0001); % the original KSC practice was to add 0.001; also used at NBER
    hLAST   = hPREV;
    rhoLAST = rhoPREV;
    [hPREV, hbar, ~, htilde] = StochVolKSCcorrAR1(logy2, hPREV, rhoPREV, hvcv, Eh0, Vh0, KSC, KSCt, ...
        Nsurvey, T, rndStream); 
    
    htilde = htilde';
    
    
    rhodraws = bayesAR1SURdraw(htilde(2:end,:), htilde(1:end-1,:), hvcv, rho0, rhoV0i, maxShake, rndStream);
    shake    = 0;
    OK       = false;
    while ~OK && shake < maxShake
        shake = shake + 1;
        OK = all(abs(rhodraws(:,shake)) < 1);
    end
    
    rhoPREV = rhodraws(:,shake);
    
    if shake == maxShake
        warning('maxShake exhausted T=%d', T)
        rhoPREV = rhoLAST;
        hPREV   = hLAST;
        continue % redo entire MCMC step
    end
    
    hresid  = htilde(2:end,:) - bsxfun(@times, htilde(1:end-1,:), rhoPREV');

    % draw hvcv
    hvcv =  bayesVCVgibbsDraw1(hvarT, hvarDof, hresid, rndStream, false);
    
    
    %% store draws post burnin
    if n > 0
        hh(:,:,n)             = hPREV;
        hhVCV(:,:,n)          = hvcv;
        SHOCKSLOPES(:,n)      = shockslopesPREV;
        RHO(:,n)              = rhoPREV;
        HBAR(:,n)             = hbar;
    end
    
    
    % prepare next loop
    n = n + 1;
    if showProgress
        progressbar((n + burnin + 1) / (burnin + MCMCdraws + 1))
    end
end % while

%% compute forecast error SV
% allocate memory
feSVdraws  = zeros(T,Nsurvey,MCMCdraws);
etaSVdraws = NaN(T,Nsurvey,MCMCdraws);
SVhorizons = horizons + 1; % We use t-1 data to project the SV's for t, t+1 etc

% progressbar(0)
for dd = 1 : MCMCdraws
    
        
    % collect MCMC parameter draws
    hVAR          = diag(hhVCV(:,:,dd));
    hbar          = HBAR(:,dd);
    htilde0       = bsxfun(@minus, hh(:,:,dd), hbar);
    rho           = RHO(:,dd);
    %     hintercept    = (1 - rho) .* hbar;
    
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
        condmean = hbar + rho.^k .* htilde0;
        rho2     = rho.^2;
        condvar  = (1 - rho2.^k) ./ (1 - rho2) .* hVAR;
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
    hshocks = randn(rndStream,Nsurvey,Nsurvey,Nfedraws);
    hhVCVsqrt = chol(hhVCV(:,:,dd))';
    for m = 1 : Nfedraws
        hshocks(:,:,m) = hhVCVsqrt * hshocks(:,:,m);
    end

    h = NaN(Nsurvey,Nsurvey,Nfedraws);
    hintercept = (1 - RHO(:,dd)) .* HBAR(:,dd);
    j = 1;
    h(:,j,:) = hintercept + RHO(:,dd) .* hh(:,end,dd) + hshocks(:,j,:);
    for j = 2 : Nsurvey
        h(:,j,:) = hintercept + RHO(:,dd) .* h(:,j-1,:) + hshocks(:,j,:);
    end
    
    % generate z and scale by SV
    z = randn(rndStream,Nsurvey,Nsurvey,Nfedraws);
    z = exp(h * .5) .* z;
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
