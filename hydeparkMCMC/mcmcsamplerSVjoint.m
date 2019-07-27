function [feSVdraws, feDraws, etaSVdraws, hhVCV] = mcmcsamplerSVjoint(eta, etaNanny, T, horizons, MCMCdraws, burnin, Nfedraws, doZeroSlopes, doDiffuseSlopes, rndStream, showProgress)


eta         = eta(1:T,:,:);
etaNanny    = etaNanny(1:T,:,:);

Nsurvey     = size(eta,2);
Ny          = size(eta,3);
Nsv         = Ny * Nsurvey;


eta      = reshape(eta, T, Nsv);
etaNanny = reshape(etaNanny, T, Nsv);


[KSC, KSCt] = getKSC10values(T, Nsv); % using 10-point grid as in Omori et al (2007, JoE)

hvarDof   = 8 + Nsv + 1; 
hvarT     = 0.2^2 * (hvarDof - Nsv - 1) * eye(Nsv);
hvarTsqrt = sqrt(0.2^2 * (hvarDof - Nsv - 1)) * eye(Nsv);

Nshockslopes  = Nsv * (Nsv - 1) / 2;
E0shockslopes = zeros(Nshockslopes, 1);
V0shockslopes = eye(Nshockslopes);



Eh0  = repmat(log(.5) * 2, Nsv,1);
Vh0  = repmat(10,Nsv,1);

hh   = NaN([Nsv T MCMCdraws]);

hhVCV = NaN(Nsv, Nsv, MCMCdraws);
% generate initial values
hvcv  = iwishdraw(hvarT, hvarDof, 1, hvarTsqrt, rndStream);
h0    = Eh0 + sqrt(Vh0) .* randn(rndStream, Nsv, 1);

SHOCKSLOPES = NaN(Nshockslopes,MCMCdraws);


SVinno = NaN(Nsv, T);

%% draw initial values
hPREV = chol(hvcv)' * randn(rndStream, Nsv, T);
hPREV = cumsum(hPREV);
hPREV = bsxfun(@plus, hPREV, h0);

if doZeroSlopes
    shockslopesPREV = zeros(Nshockslopes,1);
else
    shockslopesPREV = E0shockslopes + chol(V0shockslopes)' * randn(rndStream, Nshockslopes, 1);
end



% burnin
if showProgress
    progressbar4xterm(0)
end
if isempty(burnin)
    burnin = MCMCdraws;
end
for n = -burnin : MCMCdraws
    
    %% fill in missing values
    B      = eye(Nsv);
    offset = 0;
    for j = 2 : Nsv
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
            
            if Nmissing < Nsv
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
                missingData = randn(rndStream, 1, Nsv) * sigmaChol';
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
        for m = 2 : Nsv
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
     
     [hPREV, hPREV0] = StochVolKSCcorr(logy2, hPREV, hvcv, Eh0, Vh0, KSC, KSCt, ...
         Nsv, T, rndStream);
     
    
    % draw hvcv
    hvcv =  bayesVCVgibbsDraw1(hvarT, hvarDof, diff([hPREV0 hPREV],[],2)', rndStream, false);
    
    
    %% store draws post burnin
    if n > 0
        hh(:,:,n)        = hPREV;
        hhVCV(:,:,n)     = hvcv;
        SHOCKSLOPES(:,n) = shockslopesPREV;
    end
    
    if showProgress
        progressbar4xterm((n + burnin + 1) / (burnin +  MCMCdraws + 1), sprintf('Estimaton for T=%d:', T))
    end
end

%% compute forecast error SV
% allocate memory
etaSVdraws = NaN(T,Nsv,MCMCdraws);
feSVdraws  = zeros(T,Nsurvey,Ny,MCMCdraws);
SVhorizons = horizons + 1; % We use t-1 data to project the SV's for t, t+1 etc

% progressbar(0)
for dd = 1 : MCMCdraws
    
    % collect MCMC parameter draws
    hVAR          = diag(hhVCV(:,:,dd));
    hbar          = hh(:,:,dd)';
    
    % construct B
    B      = eye(Nsv);
    offset = 0;
    for j = 2 : Nsv
        ndx        = offset + (1:j-1);
        B(j,1:j-1) = SHOCKSLOPES(ndx,dd);
        offset     = ndx(end);
    end
    
    % construct forecasts of hbar up to maxhorizon steps ahead
    VarianceForecast  = NaN(T,Nsv,Nsurvey);
    for j = 1 : Nsurvey
        k = SVhorizons(j); 
        for n = 1 : Nsv
            VarianceForecast(:,n,j) = exp(hbar(:,n) + .5 * k * hVAR(n));
        end
    end
    
    %% more generic code:
    % 1): forecast Omega(t+j)
    Omega = NaN(Nsv,Nsv,Nsurvey,T);
    for t = 1 : T
        for j = 1 : Nsurvey
            Omega(:,:,j,t) = B * diag(VarianceForecast(t,:,j)) * B';
        end
    end
    

    % feSV
    % note: the array dimensions are not optimally ordered for performance here ...
    for t = 1 : T
        for m = 1 : Ny
            for n = 1 : Nsurvey
                %  feSVdraws(t,n,m,dd) = 0; % redundant since array already initialized to zero
                for j = 1 : n
                    this = (m-1) * Nsurvey + j;
                    feSVdraws(t,n,m,dd) = feSVdraws(t,n,m,dd) + Omega(this,this,n - j + 1,t);
                end
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
    

%% reshuffle output dimensions
etaSVdraws = reshape(etaSVdraws, [T,Nsurvey,Ny,MCMCdraws]);

%% simulate forecast errors
feDraws = NaN(Nsurvey,Nfedraws, Ny, MCMCdraws);

% model: eta = B * Vol * z
for dd = 1 : MCMCdraws
    
    % construct B
    B      = eye(Nsv);
    offset = 0;
    for j = 2 : Nsv
        ndx        = offset + (1:j-1);
        B(j,1:j-1) = SHOCKSLOPES(ndx,dd);
        offset     = ndx(end);
    end
    
    % generate SV (note the second dimension are the simulated forecasts
    h = randn(rndStream,Nsv,Nsurvey,Nfedraws);
    hhVCVsqrt = chol(hhVCV(:,:,dd))';
    for m = 1 : Nfedraws
        h(:,:,m) = hhVCVsqrt * h(:,:,m);
    end
    h = cumsum(h,2);
    h = bsxfun(@plus, h, hh(:,end,dd));
    
    % generate z and scale by SV
    z   = randn(rndStream,Nsv,Nsurvey,Nfedraws);
    z   = exp(h * .5) .* z;
    z   = reshape(z, Nsv, Nsurvey * Nfedraws);
    eta = B * z;
    eta = reshape(eta, Nsurvey, Ny, Nsurvey, Nfedraws);
    eta = permute(eta, [1 3 4 2]);
    
    % note: 1st col treated like "time when eta happens", 2nd column is the
    % horizon of eta
    % i.e. eta(t,h) is eta(t+h|t)
    
    % construct FE from eta
    for n = 1 : Ny
        for h = 0 : horizons(end)
            feDraws(h+1,:,n,dd) = 0; % this line should be redundant
            for j = 0 : h
                feDraws(h+1,:,n,dd) = feDraws(h+1,:,n,dd) + permute(eta(1+h-j,j+1,:,n), [1 3 2 4]);
            end
        end
    end
end

feDraws = permute(feDraws, [1 3 2 4]);
feDraws = reshape(feDraws, Nsurvey, Ny, Nfedraws * MCMCdraws);
