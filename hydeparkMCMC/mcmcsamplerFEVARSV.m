function [feDraws, maxlambdaDraws, hhVCV] = mcmcsamplerFEVARSV(eta, etaNanny, T, p, horizons, MCMCdraws, Nfedraws, doDiffuseSlopes, rndStream, showProgress)

% NOTE: SVdraws are *variance* draws (take sqrt to convert into Stoch*Vol*)

% NOTE: "eta" variable in this script stores actually FE ...

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
Ef0 = zeros(Nf,1);
Vf0 = zeros(Nf); % allocate memory
indexmat = reshape(1:Nf, Nreg, Nsurvey);

% see (13) of Carriero/Clark/Marcellin (2015, JAE)
% using same (i,j) notation for indices

% first lag coeff
for k = 1 : p
    for j = 1 : Nsurvey % rhs / regressor index, actually Nreg - 1; but so far only first lag implemented
        for i = 1 : Nsurvey % lhs index
            ndx = indexmat((k-1) * Nsurvey + j,i);
            if i == j
                Vf0(ndx,ndx) = lambda1^2 / k; %  * sigma2hat(i) / sigma2hat(j) = 1
            else
                Vf0(ndx,ndx) = lambda1^2 * lambda2^2 * sigma2hat(i) / sigma2hat(j) / k;
            end
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
Nx = Nsurvey * p + 1;
A  = zeros(Nx,Nx);
A(Nsurvey+1:end-1,1:Nsurvey*(p-1)) = eye(Nsurvey * (p-1));
A(end,end) = 1;
Cobs = eye(Nsurvey,Nx);
C    = repmat(Cobs, 1 , 1, T);
etaObs = eta;
etaObs(etaNanny) = 0;
etaObs = etaObs';
etaNanny = etaNanny';

for t = 1 : T
    C(etaNanny(:,t),etaNanny(:,t),t) = 0;
end

etaStateTdraws = NaN(Nx, MCMCdraws); % to storejump-off values; need to be drawn in case of missing obs

Iy = eye(Nsurvey);
Ix = eye(Nx-1);

%% prepare loop
A(1:Nsurvey,:) = reshape(fPREV, Nx, Nsurvey)';
maxlambda      = max(abs(eig(A(1:Nx-1,1:Nx-1))));    

%% loop over MCMCdraws
if showProgress
    progressbar4xterm(0)
end
burnin = MCMCdraws;
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
        eta0(1:Nx-1)                 = (Ix - A(1:Nx-1,1:Nx-1)) \ A(1:Nx-1,end);
        BSV0 = B(1:Nx-1,:) * diag(SV(:,1));
        sqrtVeta0(1:Nx-1,1:Nx-1)  = chol(dlyapdoubling(A(1:Nx-1,1:Nx-1), BSV0 * BSV0'))'; % faster than dlyapchol(A(1:Nx-1,1:Nx-1), BSV0)';
    end
    
    %% call Durbin-Koopman sampler
    [etaState, ~, etaState0] = abcDisturbanceSmoothingSamplerNaN1draw(A, B, C, etaObs, etaNanny, eta0, sqrtVeta0, ...
        [], SV, rndStream);
    
    %% remainder of MCMC
    % VARSV draws
    
    Yobs = cat(2, fliplr(reshape(etaState0(1:Nsurvey*p), Nsurvey, p)), ...
        etaState(1:Nsurvey,:))';
    
    % prepare residual variance
    iB = Iy / B(1:Nsurvey,:);
    iSigmaResid = NaN(Nsurvey,Nsurvey,T);
    for t = 1 : T
        iOmega = diag(exp(-hPREV(:,t)));
        iSigmaResid(:,:,t) = iB' * iOmega * iB;
        % checkdiff(iSigmaResid(:,:,t), inv(B(1:Nsurvey,:) * diag(exp(hPREV(:,t))) * B(1:Nsurvey,:)'));
    end
    [fPREV, SVinno, ~, maxlambda]  = bayesVARSVgibbsDraw(Yobs, p, true, iSigmaResid, Ef0, Vf0, rndStream);
    
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
    
    if showProgress
        progressbar4xterm((n + burnin + 1) / (burnin +  MCMCdraws + 1), sprintf('Estimaton for T=%d:', T))
    end
end


%% simulate forecast errors
NforecastHorizons = Nsurvey;
etaDraws          = zeros(Nsurvey,Nfedraws,NforecastHorizons,MCMCdraws); % dimensions will be reordered at the end


% model: eta = B * Vol * z
for dd = 1 : MCMCdraws
    
    % construct A
    A(1:Nsurvey,:)  = reshape(Fdraws(:,dd), Nx, Nsurvey)';
    
    % construct B
    B      = eye(Nx,Nsurvey);
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
    thisState = repmat(etaStateTdraws(:,dd), 1, Nfedraws);
    for m = 1 : NforecastHorizons
        thisState = A * thisState + B * z(:,:,m);
        etaDraws(:,:,m,dd) = Cobs * thisState;
    end
    
end


etaDraws = permute(etaDraws, [1 3, 2, 4]); % turning into [Nsurvey, NforecastHorizons, Nfedraws, MCMCdraws]
etaDraws = reshape(etaDraws, Nsurvey, NforecastHorizons, Nfedraws * MCMCdraws);

% feDraws  = reshape(feDraws, Nsurvey, Nfedraws * MCMCdraws);
feDraws = NaN(Nsurvey, Nfedraws * MCMCdraws);
for n = 1 : Nsurvey
    feDraws(n,:) = squeeze(etaDraws(n,1+horizons(n),:));
end