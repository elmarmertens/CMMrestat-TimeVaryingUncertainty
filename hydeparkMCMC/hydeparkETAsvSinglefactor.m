%%  ETA-SV model but w/single-factor SV specification
% uses MDS assumption
% computes also FE-SIMPLE alternative;
% ... results for the latter are labeled "RT..." as the underlying computations
% ... intend to mimic the approach of Reifschneider and Tulip (2007, 2017)

%% load toolboxes
path(pathdef)

addpath ../toolbox/emtools/
addpath ../toolbox/emtexbox/
addpath ../toolbox/emgibbsbox/
addpath ../toolbox/emeconometrics/


%% clear workspace
clear variables
clear global
close all
fclose all;
clc

%% some parameters

Nstreams    = max(1,getparpoolsize);
rndStreams  = initRandStreams(Nstreams, [], 0);


quicky     = false; % if TRUE: very short MCMC chains, no looping across variables,
%  useful for testing code execution, see below for specific settings
doStore    = true; % to store result files in directory specified by "localstore" m-file
doStoreXXL = false; % to store extended result files (with MCMC draws)
DATALABELS = {'CPI', 'TBILL', 'UNRATE', 'PGDP', 'RGDP'};


Tstart     = 60;
RTwindow   = Tstart;
evalTstart = Tstart;

doDiffuseSlopes = true;
doZeroSlopes = false;

tic

%% loop over variables
for d =  1 : length(DATALABELS)
    
    
    
    close all
    datalabel = DATALABELS{d};
    
    if quicky
        MCMCdraws    = 1e2;
        Nfedraws     = 10;
    else
        MCMCdraws    = 1e4;
        Nfedraws     = 100;
    end
    
    % construct a model-label indicating dataset and important parameters, to
    % append to picture names
    modellabel = datalabel;
    %#ok<*UNRCH>
    if doZeroSlopes
        modellabel = strcat(modellabel, '00');
    else
        modellabel = strcat(modellabel, '01');
    end
    
    modellabel = strcat(modellabel,'QRTdensitySVsinglefactor');
    
    modellabel = strcat(modellabel, sprintf('TSTART%dEVALSTART%d', Tstart, evalTstart));
    
    
    wrap = [];
    titlename = modellabel;
    if ~quicky
        initwrap
    end
    
    fprintf('Processing %s ... \n', modellabel)
    
    %% load data
    matfile = sprintf('hydepark%sdata', upper(datalabel));
    load(matfile, 'Yforecast', 'Yrealized', 'eta', 'FE', 'FEobserved', 'Nsurvey', 'dates', 'horizons', 'maxHorizon', 'Ylabel')
    
    
    %% record location of missing values (NaN)
    etaNanny = isnan(eta);
    
    samStart = find(~(etaNanny(:,1)),1);
    
    eta        = eta(samStart:end, :);
    etaNanny   = etaNanny(samStart:end, :);
    FE         = FE(samStart:end, :);
    FEobserved = FEobserved(samStart:end, :);
    dates      = dates(samStart:end);
    T          = length(dates);
    
    
    %% plot eta elements
    etavar = nanvar(eta);
    for n = 1 : Nsurvey
        figure
        hold on
        plot(dates,eta(:,n))
        plothorzline(2 * sqrt(etavar(n)), [], 'r--')
        plothorzline(-2 * sqrt(etavar(n)), [], 'r--')
        
        plothorzline(nanmean(eta(:,n)), [], 'r-.')
        
        plothorzline(0, [], 'k:')
        title(sprintf('eta h=%d', n-1))
        xtickdates(dates)
        wrapcf(sprintf('eta%d%s', n, modellabel), wrap)
        
    end
    
    %% plot FORECAST error
    FEvar   = nanvar(FE);
    FEetavar = cumsum(etavar) ;
    for n = 1 : Nsurvey
        figure
        hold on
        plot(dates,FE(:,n))
        plothorzline(2 * sqrt(FEvar(n)), [], 'r--')
        plothorzline(-2 * sqrt(FEvar(n)), [], 'r--')
        
        
        plothorzline(2 * sqrt(FEetavar(n)), [], 'b-.')
        plothorzline(-2 * sqrt(FEetavar(n)), [], 'b-.')
        
        plothorzline(nanmean(FE(:,n)), [], 'r-.')
        
        plothorzline(0, [], 'k:')
        title(sprintf('FE h=%d', n-1))
        xtickdates(dates)
        wrapcf(sprintf('fe%d%s', n, modellabel), wrap)
        
    end
    
    close all
    
    %% compute forecast error SV
    % allocate memory
    feSVdraws       = NaN(T,Nsurvey,MCMCdraws);
    etaSVdraws      = NaN(T,Nsurvey,MCMCdraws);
    feDraws         = NaN(T,Nsurvey,MCMCdraws*Nfedraws);
    feVol           = NaN(T,Nsurvey); % to store SE of feDraws
    hbar            = NaN(T, Nsurvey, MCMCdraws);
    betaSV          = NaN(T, Nsurvey, MCMCdraws);
    
    parfor thisT = Tstart : T - 1 % @parfor
        TID = parid;
        [mcmc_feSV, mcmc_fe, mcmc_etaSV, mcmc_hbar, mcmc_beta] = mcmcsamplerSVsinglefactor(eta, etaNanny, thisT, horizons, MCMCdraws, Nfedraws, doZeroSlopes, doDiffuseSlopes, rndStreams{TID}, false); %#ok<PFBNS>
        feSVdraws(thisT,:,:)             = mcmc_feSV(end,:,:);
        etaSVdraws(thisT,:,:)            = mcmc_etaSV(end,:,:);
        feDraws(thisT,:,:)               = mcmc_fe;
        feVol(thisT,:)                   = std(mcmc_fe,0,2);
        hbar(thisT,:,:)                  = mcmc_hbar;
        betaSV(thisT,:,:)                = mcmc_beta;
        fprintf('%s - QRT: done with t=%d\n', modellabel, thisT)
    end
    
    TID = 1;
    % doing the last step outside of the parfor loop (to grab final
    % draws)
    thisT                            = T;
    fprintf('%s - QRT proceeding to T=%d\n', modellabel, thisT)
    
    [mcmc_feSV, mcmc_fe, mcmc_etaSV, mcmc_hbar, mcmc_beta] = mcmcsamplerSVsinglefactor(eta, etaNanny, thisT, horizons, MCMCdraws, Nfedraws, doZeroSlopes, doDiffuseSlopes, rndStreams{TID}, true);
    feSVdraws(thisT,:,:)             = mcmc_feSV(end,:,:);
    etaSVdraws(thisT,:,:)            = mcmc_etaSV(end,:,:);
    feDraws(thisT,:,:)               = mcmc_fe;
    feVol(thisT,:)                   = std(mcmc_fe,0,2);
    hbar(thisT,:,:)                  = mcmc_hbar;
    betaSV(thisT,:,:)                = mcmc_beta;
    fprintf('%s - QRT: done with t=%d\n', modellabel, thisT)
    
    feSVFINALdraws                    = mcmc_feSV;
    etaSVFINALdraws                   = mcmc_etaSV;
    
    
    feSV      = sqrt(mean(feSVdraws,3));
    feSVmed   = sqrt(median(feSVdraws,3));
    feSVtails = sqrt(prctile(feSVdraws, [25 75], 3));
    
    
    feSVFINAL      = squeeze(sqrt(mean(feSVFINALdraws,3)));
    feSVFINALmed   = squeeze(sqrt(median(feSVFINALdraws,3)));
    feSVFINALtails = squeeze(sqrt(prctile(feSVFINALdraws, [25 75], 3)));
    
    etaSV      = sqrt(mean(etaSVdraws,3));
    etaSVmed   = sqrt(median(etaSVdraws,3));
    etaSVtails = sqrt(prctile(etaSVdraws, [25 75], 3));
    
    
    etaSVFINAL      = squeeze(sqrt(mean(etaSVFINALdraws,3)));
    etaSVFINALmed   = squeeze(sqrt(median(etaSVFINALdraws,3)));
    etaSVFINALtails = squeeze(sqrt(prctile(etaSVFINALdraws, [25 75], 3)));
    
    feMean = mean(feDraws, 3);
    feTails = prctile(feDraws, [.05 2.5 5 25 75 95 97.5 99.5], 3);
    
    clear mcmc_*
    
    %% set FE over estimation window to NaN
    FE(1:Tstart-1,:)      = NaN;
    eta(1:Tstart-1,:)     = NaN;
    % evalTstart should always be larger than Tstart, thus making the previous two lines redundant
    FE(1:evalTstart-1,:)  = NaN;
    eta(1:evalTstart-1,:) = NaN;

    

    %% compute coverage rates
    volScales = 1:2;
    [coverageTstat, coveragePval, coverageMean] = deal(NaN(Nsurvey,length(volScales)));
    for v = 1 : length(volScales)
        coverageNdx = double(abs(FE) < volScales(v) * feSV);
        coverageNdx(isnan(feSV)) = NaN;
        coverageNdx(isnan(FE))   = NaN;
        
        threshold   = normcdf(volScales(v)) - normcdf(-volScales(v));
        
        coverageMean(:,v) = nanmean(coverageNdx);
        
        for n = 1 : Nsurvey
            coverrate = coverageNdx(:,n) - threshold;
            nanny = ~isnan(coverrate);
            reggae = nwest(coverrate(nanny), ones(sum(nanny),1), ceil(1.5 + horizons(n)));
            %                                 hrulefill
            %                                 fprintf('Horizon %d, Coverage rate for %d times SV is %6.4f \n', horizons(n), volScales(v), coverageMean(n,v))
            %                                 prt(reggae)
            coverageTstat(n,v) = reggae.tstat;
            coveragePval(n,v)  = tdis_prb(reggae.tstat,reggae.nobs-reggae.nvar);
        end
    end
    
    % same for eta
    [ETAcoverageTstat, ETAcoveragePval, ETAcoverageMean] = deal(NaN(Nsurvey,length(volScales)));
    for v = 1 : length(volScales)
        ETAcoverageNdx = double(abs(eta) < volScales(v) * etaSV);
        ETAcoverageNdx(isnan(etaSV)) = NaN;
        ETAcoverageNdx(isnan(eta))   = NaN;
        
        threshold   = normcdf(volScales(v)) - normcdf(-volScales(v));
        
        ETAcoverageMean(:,v) = nanmean(ETAcoverageNdx);
        
        for n = 1 : Nsurvey
            coverrate = ETAcoverageNdx(:,n) - threshold;
            nanny = ~isnan(coverrate);
            reggae = nwest(coverrate(nanny), ones(sum(nanny),1), ceil(1.5 + horizons(n)));
            %                                 hrulefill
            %                                 fprintf('Horizon %d, ETA-Coverage rate for %d times SV is %6.4f \n', horizons(n), volScales(v), ETAcoverageMean(n,v))
            %                                 prt(reggae)
            ETAcoverageTstat(n,v) = reggae.tstat;
            ETAcoveragePval(n,v)  = tdis_prb(reggae.tstat,reggae.nobs-reggae.nvar);
        end
    end
    
    
    %% compute quantile-coverage-rates
    quantileP = [2.5 5 normcdf(-1) * 100 25 50 75 normcdf(1) * 100 95 97.5];
    criticalQuants = prctile(feDraws, quantileP, 3);
    
    fee                       = repmat(FE, [1 1 length(quantileP)]);
    quantcoverNdx             = double(fee < criticalQuants); % bsxfun would work, but fee needed for next step
    quantcoverNdx(isnan(fee) | isnan(criticalQuants)) = NaN;
    quantcoverMean            = squeeze(nanmean(quantcoverNdx,1));
    
    [quantcoverTstat, quantcoverPval] =deal(NaN(size(quantcoverMean)));
    
    % test for significance with Newey-West regressions
    for q = 1 : length(quantileP)
        
        for n = 1 : Nsurvey
            coverrate = quantcoverNdx(:,n,q) - quantileP(q) / 100;
            nanny = ~isnan(coverrate);
            reggae = nwest(coverrate(nanny), ones(sum(nanny),1), ceil(1.5 + horizons(n)));
            %                                 hrulefill
            %                                 fprintf('Horizon %d, Quantile-Coverage rate for %4.2f%% is %4.2f%% \n', horizons(n), quantileP(q), quantcoverMean(n,q) * 100);
            %                                 prt(reggae)
            quantcoverTstat(n,q) = reggae.tstat;
            quantcoverPval(n,q)  = tdis_prb(reggae.tstat,reggae.nobs-reggae.nvar);
        end
    end
    
    
    %% compute CRPS
    crps  = NaN(T,Nsurvey);
    feDraws = permute(feDraws, [3 1 2]); % should make the following loop quicker
    parfor n = 1 : Nsurvey
        for t = 1 : T % could just start at Tstart
            crps(t,n) = crpsDraws(FE(t,n), feDraws(:,t,n));
        end % t
    end
    if doStoreXXL
        feDraws = permute(feDraws, [2 3 1]);
    else
        clear feDraws
    end
    
    
    %% compute RT SV bands
    RTFEsv     = NaN(T,Nsurvey);
    
    for thisT = Tstart : T % @parfor
        theseFE           = FEobserved(thisT-RTwindow+1:thisT,:);
        RTFEsv(thisT,:)   = sqrt(nanmean(theseFE.^2));
        %                 fprintf('%s - RT : done with t=%d\n', modellabel, thisT)
    end
    
    RTFEdraws = randn(T,Nsurvey,MCMCdraws * Nfedraws);
    RTFEdraws = bsxfun(@times, RTFEdraws, RTFEsv);
    
    
    
    %% compute RT coverage rates
    volScales = 1:2;
    [RTcoverageTstat, RTcoveragePval, RTcoverageMean] = deal(NaN(Nsurvey,length(volScales)));
    for v = 1 : length(volScales)
        coverageNdx = double(abs(FE) < volScales(v) * RTFEsv);
        coverageNdx(isnan(feSV)) = NaN;
        coverageNdx(isnan(FE))   = NaN;
        
        if any(isnan(feSV) & ~isnan(FE))
            warning('mismatch between feSV and FE?')
        end
        threshold   = normcdf(volScales(v)) - normcdf(-volScales(v));
        
        RTcoverageMean(:,v) = nanmean(coverageNdx);
        
        for n = 1 : Nsurvey
            coverrate = coverageNdx(:,n) - threshold;
            nanny = ~isnan(coverrate);
            reggae = nwest(coverrate(nanny), ones(sum(nanny),1), ceil(1.5 + horizons(n)));
            %                                 hrulefill
            %                                 fprintf('Horizon %d, Coverage rate for %d times SV is %6.4f \n', horizons(n), volScales(v), RTcoverageMean(n,v))
            %                                 prt(reggae)
            RTcoverageTstat(n,v) = reggae.tstat;
            RTcoveragePval(n,v)  = tdis_prb(reggae.tstat,reggae.nobs-reggae.nvar);
        end
    end
    
    %% RT - quantiles
    quantileP = [2.5 5 normcdf(-1) * 100 25 50 75 normcdf(1) * 100 95 97.5];
    criticalQuants = prctile(RTFEdraws, quantileP, 3);
    
    fee                       = repmat(FE, [1 1 length(quantileP)]);
    quantcoverNdx             = double(fee < criticalQuants); % bsxfun would work, but fee needed for next step
    quantcoverNdx(isnan(fee) | isnan(criticalQuants)) = NaN;
    RTquantcoverMean        = squeeze(nanmean(quantcoverNdx,1));
    
    [RTquantcoverTstat, RTquantcoverPval] =deal(NaN(size(RTquantcoverMean)));
    
    % test for significance with Newey-West regressions
    for q = 1 : length(quantileP)
        
        for n = 1 : Nsurvey
            coverrate = quantcoverNdx(:,n,q) - quantileP(q) / 100;
            nanny = ~isnan(coverrate);
            reggae = nwest(coverrate(nanny), ones(sum(nanny),1), ceil(1.5 + horizons(n)));
            %                                 hrulefill
            %                                 fprintf('Horizon %d, Quantile-Coverage rate for %4.2f%% is %4.2f%% \n', horizons(n), quantileP(q), RTquantcoverMean(n,q) * 100);
            %                                 prt(reggae)
            RTquantcoverTstat(n,q) = reggae.tstat;
            RTquantcoverPval(n,q)  = tdis_prb(reggae.tstat,reggae.nobs-reggae.nvar);
        end
    end
    
    
    
    %% RT - CRPS
    RTcrps    = NaN(T,Nsurvey);
    RTFEdraws = permute(RTFEdraws, [3 1 2]); % should make the following loop quicker
    % progressbar(0)
    parfor n = 1 : Nsurvey
        for t = 1 : T % could just start at Tstart
            RTcrps(t,n) = crpsDraws(FE(t,n), RTFEdraws(:,t,n));
        end % t
    end % n
    clear RTFEdraws
    
    if ~isequal(isnan(RTcrps), isnan(FE))
        error('crps has NaN when there is data')
    end
    
    
    
    %% plot etaSV distro with abs eta
    for n = 1 : Nsurvey
        hanni = NaN(3,1);
        
        figure
        hold on
        barcol = .3 * [0 1 0];
        hanni(1) = bar(dates, abs(eta(:,n)), 'FaceColor', barcol, 'EdgeColor', barcol);
        
        hanni(2) = plot(dates, etaSV(:,n), 'k-', 'linewidth', 3);
        plot(dates, squeeze(etaSVtails(:,n,:)), 'k-', 'linewidth', 1)
        
        hanni(3) = plot(dates, etaSVFINAL(:,n), 'r--', 'linewidth', 3);
        plot(dates, squeeze(etaSVFINALtails(:,n,:)), 'r--', 'linewidth', 1)
        
        maxy = max(ylim);
        xtickdates(dates(evalTstart-1:end))
        ylim([0 maxy])
        
        legend(hanni, 'abs(eta)', 'QRT', 'Final')
        
        title(sprintf('ETA (h=%d): QRT Coverage Rate 70: %6.4f (p: %4.2f); 95: %6.4f (p: %4.2f)', horizons(n), ETAcoverageMean(n,1), coveragePval(n,1), coverageMean(n,2), coveragePval(n,2)))
        wrapcf(sprintf('etaSViqr%d%s', n, modellabel), wrap)
    end
    
    
    
    %% plot feSV distro with abs FE
    % close all
    for n = 1 : Nsurvey
        hanni = NaN(5,1);
        
        figure
        hold on
        barcol = .5 * [1 1 1];
        hanni(1) = bar(dates, abs(FE(:,n)), 'FaceColor', barcol, 'EdgeColor', barcol);
        
        hanni(2) = plot(dates, feSV(:,n), 'k-', 'linewidth', 3);
        plot(dates, squeeze(feSVtails(:,n,:)), 'k-', 'linewidth', 1)
        
        hanni(3) = plot(dates, feSVFINAL(:,n), 'r--', 'linewidth', 3);
        plot(dates, squeeze(feSVFINALtails(:,n,:)), 'r--', 'linewidth', 1)
        
        hanni(4) = plot(dates, feVol(:,n), 'y--', 'linewidth', 2);
        
        hanni(5) = plot(dates, RTFEsv(:,n), 'b-', 'linewidth', 3);
        
        maxy = max(ylim);
        xtickdates(dates(evalTstart-1:end))
        ylim([0 maxy])
        
        legend(hanni, 'abs(FE)', 'QRT', 'Final', 'QRT (montecarlo', 'Reifschneider-Tulip')
        
        title(sprintf('FE (h=%d): QRT Coverage Rate 70: %6.4f (p: %4.2f); 95: %6.4f (p: %4.2f)', horizons(n), coverageMean(n,1), coveragePval(n,1), coverageMean(n,2), coveragePval(n,2)))
        % title(sprintf('FEh=%d: QRT Coverage Rate 70: %6.4f (p: %4.2f); 95: %6.4f (p: %4.2f)', horizons(n), coverageMean(n,1), coveragePval(n,1), coverageMean(n,2), coveragePval(n,2)))
        title(sprintf('FEh (h=%d) -- Coverage Rate 70, SV/RT: %6.4f/%6.4f (p: %4.2f/%4.2f); 95: %6.4f/%6.4f (p: %4.2f/%4.2f)', ...
            horizons(n), coverageMean(n,1), RTcoverageMean(n,1), ...
            coveragePval(n,1), ...
            RTcoveragePval(n,1), ...
            coverageMean(n,2), RTcoverageMean(n,2), ...
            coveragePval(n,2), ...
            RTcoveragePval(n,2)))
        wrapcf(sprintf('feSViqr%d%s', n, modellabel), wrap)
    end
    
    %% compare CRPS
    [avgCRPS, RTavgCRPS] = deal(NaN(T,Nsurvey));
    for t = 1 : T
        avgCRPS(t,:) = nanmean(crps(1:t,:));
        RTavgCRPS(t,:) = nanmean(RTcrps(1:t,:));
    end
    
    for n = 1 : Nsurvey
        figure
        subplot(2,1,1)
        hanni = NaN(2,1);
        hold on
        hanni(1) = plot(dates, RTavgCRPS(:,n), 'k-', 'linewidth', 1);
        hanni(2) = plot(dates, avgCRPS(:,n), 'r--', 'linewidth', 2);
        plotOrigin
        xtickdates(dates(evalTstart-1:end))
        legend(hanni, 'RT', 'SV', 'location', 'best')
        title(sprintf('CRPS (RT - SV) / RT: %5.2f%% ', ...
            (RTavgCRPS(end,n) - avgCRPS(end,n))/RTavgCRPS(end,n) * 100))
        
        subplot(2,1,2)
        hanni = NaN(2,1);
        hold on
        hanni(1) = plot(dates, RTcrps(:,n), 'k-');
        hanni(2) = plot(dates, crps(:,n), 'r--');
        xtickdates(dates(evalTstart-1:end))
        legend(hanni, 'RT', 'SV')
        title('contributions')
        
        wrapcf(sprintf('crps%d%s', n, modellabel), wrap)
    end
    

    %% plot factor loadings
    figure
    plot(dates, median(betaSV,3), 'linewidth', 2)
    legend(Ylabel{:}, 'location', 'best')
    xtickdates(dates(evalTstart-1:end))
    wrapcf(sprintf('betaSV%s', modellabel), wrap)
    
    figure
    plot(dates, median(exp(hbar / 2),3), 'linewidth', 2)
    legend(Ylabel{:}, 'location', 'best')
    xtickdates(dates(evalTstart-1:end))
    ylimit = ylim;
    ylim([0 ylimit(2)])
    wrapcf(sprintf('barSV%s', modellabel), wrap)
    
    %% store results
    diary off
    if  doStore && ~quicky
        if doStoreXXL
            save(fullfile(localstore, modellabel), ...
                '-v7.3');
        end
        clear foo feSVdraws feSVFINALdraws feDraws etaSVdraws etaSVFINALdraws RTFEdraws
        save(fullfile(localstore, strcat('slim', modellabel)), ...
            '-v7.3');
    end
    finishwrap % if latex compilation fails, edit this script to turn the compilation off
    if quicky
        dockAllFigures
        toc
        return
    else
        close all
    end
end

toc

%% finish / clean up
finishscript
dockAllFigures
