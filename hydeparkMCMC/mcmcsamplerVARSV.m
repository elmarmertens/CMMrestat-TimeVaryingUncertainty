function [feDraws, etaDraws, Fdraws, maxlambdaDraws, hhVCV] = mcmcsamplerVARSV(eta, etaNanny, T, horizons, MCMCdraws, Nfedraws, doDiffuseSlopes, rndStream)

% NOTE: SVdraws are *variance* draws (take sqrt to convert into Stoch*Vol*)

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

SHOCKSLOPES = NaN(Nshockslopes,MCMCdraws);


%% setup VARSV priors
p    = 1;
Nreg = 1 + Nsurvey * p;
Nf   = Nsurvey * Nreg;

% tightness of own slope priors
lambda1 = 0.2;
% tightness of other slope priors
lambda2 = 0.5;
% tightness of steady state prior (Sims' unit root prior)
lambda0 = 1;


%% get estimates of residual variances
sigma2hat = NaN(Nsurvey,1);
for n = 1 : Nsurvey
    y = eta(~etaNanny(:,n),n);
    reggae = ols(y(2:end), [y(1:end-1) ones(length(y)-1,1)]);
    sigma2hat(n) = reggae.sige;
end

%% Minnesota prior for VAR
if p > 1
    error('need to implement prior for higher lags') % see mcmcsamplerFEVARSV.m for a solution
end
Ef0 = zeros(Nf,1);
Vf0 = zeros(Nf); % allocate memory
indexmat = reshape(1:Nf, Nreg, Nsurvey);

% see (13) of Carriero/Clark/Marcellino (2015, JAE)
% using same (i,j) notation for indices

% first lag coeff
for j = 1 : Nsurvey % rhs / regressor index, actually Nreg - 1; but so far only first lag implemented
    for i = 1 : Nsurvey % lhs index
        ndx = indexmat(j,i);
        if i == j 
            Vf0(ndx,ndx) = lambda1^2; %  * sigma2hat(i) / sigma2hat(j) = 1
        else
            Vf0(ndx,ndx) = lambda1^2 * lambda2^2 * sigma2hat(i) / sigma2hat(j);
        end
    end
end
% intercept
for i = 1 : Nsurvey
    ndx = indexmat(end,i);
    Vf0(ndx,ndx) = lambda0^2 * sigma2hat(i);
end

% note: the choleski in the next step will also tell whether all diagonals
% of Vf0 have been set to non-zero values
cholVf0 = chol(Vf0)';


Fdraws = NaN(Nf, MCMCdraws);
maxlambdaDraws = NaN(1,MCMCdraws);

%% draw initial values
hPREV = chol(hvcv)' * randn(rndStream, Nsurvey, T);
hPREV = cumsum(hPREV);
hPREV = bsxfun(@plus, hPREV, h0);

shockslopesPREV = E0shockslopes + chol(V0shockslopes)' * randn(rndStream, Nshockslopes, 1);

fPREV = Ef0 + cholVf0 * randn(rndStream, Nf, 1);

%% set up state space to sample missing values
Nx = Nsurvey + 1;
A  = zeros(Nx,Nx);
A(end,end) = 1;
C  = repmat(eye(Nsurvey,Nx), 1 , 1, T);
etaObs = eta;
etaObs(etaNanny) = 0;
etaObs = etaObs';
etaNanny = etaNanny';

for t = 1 : T
    C(etaNanny(:,t),etaNanny(:,t),t) = 0;
end

etaStateTdraws = NaN(Nx, MCMCdraws); % to storejump-off values; need to be drawn in case of missing obs

Iy = eye(Nsurvey);

%% prepare loop
A(1:Nsurvey,:) = reshape(fPREV, Nx, Nsurvey)';
maxlambda      = max(abs(eig(A(1:Nx-1,1:Nx-1))));    

%% loop over MCMCdraws
for n = -MCMCdraws : MCMCdraws
    
    
    %% update VARSV transition
    A(1:Nsurvey,:) = reshape(fPREV, Nx, Nsurvey)';
    
    %% update B
    B      = eye(Nx,Nsurvey);
    offset = 0;
    for j = 2 : Nsurvey
        ndx        = offset + (1:j-1);
        B(j,1:j-1) = shockslopesPREV(ndx);
        offset     = ndx(end);
    end
    
    SV = exp(hPREV * 0.5);
    
        
    
    %% update priors about initial eta
    eta0               = zeros(Nx,1);
    sqrtVeta0          = 10 * eye(Nx);
    eta0(end)          = 1;
    sqrtVeta0(end,end) = 0;
    if maxlambda < 1
        eta0(1:Nsurvey)                 = (Iy - A(1:Nsurvey,1:Nsurvey)) \ A(1:Nsurvey,end);
        BSV0 = B(1:Nsurvey,:) * diag(SV(:,1));
        sqrtVeta0(1:Nsurvey,1:Nsurvey)  = chol(dlyapdoubling(A(1:Nsurvey,1:Nsurvey), BSV0 * BSV0'))'; 
    end
    
    %% call Durbin-Koopman sampler
    [etaState, ~, etaState0] = abcDisturbanceSmoothingSamplerNaN1draw(A, B, C, etaObs, etaNanny, eta0, sqrtVeta0, ...
        [], SV, rndStream);
    
    %% remainder of MCMC
    % VARSV draws
    
    Yobs = cat(2, etaState0(1:Nsurvey), etaState(1:Nsurvey,:))';
    
    % prepare residual variance
    iB = Iy / B(1:Nsurvey,:);
    iSigmaResid = NaN(Nsurvey,Nsurvey,T);
    for t = 1 : T
        iOmega = diag(exp(-hPREV(:,t)));
        iSigmaResid(:,:,t) = iB' * iOmega * iB;
        % checkdiff(iSigmaResid(:,:,t), inv(B(1:Nsurvey,:) * diag(exp(hPREV(:,t))) * B(1:Nsurvey,:)'));
    end
    [fPREV, SVinno, ~, maxlambda]  = bayesVARSVgibbsDraw(Yobs, 1, true, iSigmaResid, Ef0, Vf0, rndStream);
    
    if (maxlambda > 1) && (n > 0)
        warning('unstable VAR')
    end
    
    SVinno = SVinno';
    
    %% VAR resid VCV
    
    % run regressions to create Choleski of eta
    offset = 0;
    for m = 2 : Nsurvey
        ndx = offset + (1:m-1);
        rhs = bsxfun(@rdivide, SVinno(1:m-1,:), SV(m,:))';
        lhs = (SVinno(m,:) ./ SV(m,:))';
        if doDiffuseSlopes
            [shockslopesPREV(ndx),resid] = bayesRegressionSlopesGibbsDrawDiffuse(lhs, rhs, 1, 1, rndStream);
        else
            [shockslopesPREV(ndx),resid] = bayesRegressionSlopesGibbsDraw(lhs, rhs, 1, E0shockslopes(ndx), V0shockslopes(ndx,ndx), 1, rndStream);
        end
        SVinno(m,:) = resid' .* SV(m,:);
        offset = ndx(end);
    end
    
    
    % KSC draw
    logy2           = log(SVinno.^2 + 0.0001); % the original KSC practice was to add 0.001; we are using however 10-point grid as in Omori et al. (2007, JoE)
    
    [hPREV, hPREV0] = StochVolKSCcorr(logy2, hPREV, hvcv, Eh0, Vh0, KSC, KSCt, ...
        Nsurvey, T, rndStream);
    
    
    % draw hvcv
    hvcv =  bayesVCVgibbsDraw1(hvarT, hvarDof, diff([hPREV0 hPREV],[],2)', rndStream, false);
    
    
    %% store draws post burnin
    if n > 0
        hh(:,:,n)             = hPREV;
        hhVCV(:,:,n)          = hvcv;
        SHOCKSLOPES(:,n)      = shockslopesPREV;
        Fdraws(:,n)           = fPREV;
        maxlambdaDraws(n)     = maxlambda;
        etaStateTdraws(:,n)   = etaState(:,end);
    end
    
    % progressbar((n + MCMCdraws + 1) / (2 * MCMCdraws + 1))
end


%% simulate forecast errors
feDraws  = zeros(Nsurvey,Nfedraws,MCMCdraws);
NforecastHorizons = 12;
etaDraws = zeros(Nsurvey,Nfedraws,NforecastHorizons,MCMCdraws); % dimensions will be reordered at the end


% model: eta = B * Vol * z
for dd = 1 : MCMCdraws
    
    % construct A
    A  = reshape(Fdraws(:,dd), Nx, Nsurvey)';
    
    % construct B
    B      = eye(Nsurvey);
    offset = 0;
    for j = 2 : Nsurvey
        ndx        = offset + (1:j-1);
        B(j,1:j-1) = SHOCKSLOPES(ndx,dd);
        offset     = ndx(end);
    end
    
    % generate SV (note the second dimension are the simulated forecasts
    h = randn(rndStream,Nsurvey,NforecastHorizons,Nfedraws);
    hhVCVsqrt = chol(hhVCV(:,:,dd))';
    for m = 1 : Nfedraws
        h(:,:,m) = hhVCVsqrt * h(:,:,m);
    end
    h = cumsum(h,2);
    h = bsxfun(@plus, h, hh(:,end,dd));
    
    h = permute(h, [1 3 2]);
    SV = exp(h * 0.5);
    
    % generate shocks z and scale by SV;
    z = SV .* randn(rndStream, Nsurvey, Nfedraws, NforecastHorizons);
    
    % simulate VAR
    
    m = 1;
    etaDraws(:,:,m,dd) = bsxfun(@plus, A * etaStateTdraws(:,dd), B * z(:,:,m));
    for m = 2 : NforecastHorizons
        etaDraws(:,:,m,dd) = A * cat(1, etaDraws(:,:,m-1,dd), ones(1,Nfedraws)) + B * z(:,:,m);
    end
    
    % construct FE from eta
    for h = 0 : horizons(end)
        feDraws(h+1,:,dd) = 0; % this line should be redundant
        for j = 0 : h
            feDraws(h+1,:,dd) = feDraws(h+1,:,dd) + etaDraws(1+h-j,:,j+1,dd);
        end
    end
    
end

feDraws  = reshape(feDraws, Nsurvey, Nfedraws * MCMCdraws);

etaDraws = permute(etaDraws, [1 3, 2, 4]); % turning into [Nsurvey, NforecastHorizons, Nfedraws, MCMCdraws]
etaDraws = reshape(etaDraws, Nsurvey, NforecastHorizons, Nfedraws * MCMCdraws);
