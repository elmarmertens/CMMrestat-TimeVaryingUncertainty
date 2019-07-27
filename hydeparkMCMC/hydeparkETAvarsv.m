%% ETA-VAR-SV model
% etas modeled as VAR-SV without MDS assumption 
% computes also FE-SIMPLE alternative;
% ... results for the latter are labeled "RT..." as the underlying computations
% ... intend to mimic the approach of Reifschneider and Tulip (2007, 2017)

%% load toolboxes
path(pathdef)

addpath ../toolbox/emtools/
addpath ../toolbox/emtexbox/
addpath ../toolbox/emgibbsbox/
addpath ../toolbox/emstatespace/
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
    modellabel = strcat(modellabel, '01'); % legacy naming convention
    
    modellabel = strcat(modellabel,'QRTdensityVARSV');
    
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
    
    
    
    
    
    %% compute forecast error SV
    % allocate memory
    
    feDraws         = NaN(T,Nsurvey,MCMCdraws*Nfedraws);
    feRMSE          = NaN(T,Nsurvey); % to store SE of feDraws
    
    hVCVdraws       = NaN(Nsurvey,Nsurvey,MCMCdraws,T);
    maxlambdaDraws  = NaN(1,MCMCdraws,T);
    
    
    parfor thisT = Tstart : T % @parfor
        
        TID = parid;
        [mcmc_fe, ~, ~, mcmc_maxlambda, mcmc_hVCV] = mcmcsamplerVARSV(eta, etaNanny, thisT, horizons, MCMCdraws, Nfedraws, doDiffuseSlopes, rndStreams{TID}); %#ok<PFBNS>
        feDraws(thisT,:,:)               = mcmc_fe;
        feRMSE(thisT,:)                  = sqrt(nanmean(mcmc_fe.^2,2));
        hVCVdraws(:,:,:,thisT)           = mcmc_hVCV;
        maxlambdaDraws(:,:,thisT)        = mcmc_maxlambda;
        fprintf('%s - QRT: done with t=%d\n', modellabel, thisT)
    end
    
    
    feMean = mean(feDraws, 3);
    feTails = prctile(feDraws, [.05 2.5 5 25 75 95 97.5 99.5], 3);
    
    clear mcmc_*
    
    %% set FE over estimation window to NaN
    FE(1:Tstart-1,:)      = NaN;
    eta(1:Tstart-1,:)     = NaN;
    % evalTstart should always be larger than Tstart, thus making the previous two lines redundant
    FE(1:evalTstart-1,:)  = NaN;
    eta(1:evalTstart-1,:) = NaN;
    
    %% document correlation between SV shocks
    hVCVeigen = NaN(Nsurvey,MCMCdraws,T);
    for thisT = Tstart+ 1 : T
        for n = 1 : MCMCdraws
            hVCVeigen(:,n,thisT) = sort(eig(hVCVdraws(:,:,n,thisT)), 'descend');
        end
    end
    
    tails = prctile(hVCVeigen, [5 95], 2);
    figure
    for n = 1 : Nsurvey
        hold on
        plot(dates, squeeze(mean(hVCVeigen(n,:,:), 2)), '-', 'color', Colors4Plots(n), 'linewidth', 3)
        plot(dates, permute(tails(n,:,:), [3 2 1]), '-', 'color', Colors4Plots(n))
        ylim([0 max(ylim)])
        xtickdates(dates([Tstart T]))
    end
    wrapcf(sprintf('hVCVeigenvalues%s', modellabel), wrap)
    
    %% compute coverage rates
    volScales = 1:2;
    [coverageTstat, coveragePval, coverageMean] = deal(NaN(Nsurvey,length(volScales)));
    for v = 1 : length(volScales)
        coverageNdx = double(abs(FE) < volScales(v) * feRMSE);
        coverageNdx(isnan(feRMSE)) = NaN;
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
        coverageNdx(isnan(FE))   = NaN;
        
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
    
    
        
    %% just SV with 70 and 95 CI
    hanni = NaN(2,1);
    for n = 1 : Nsurvey
        figure
        %     subplot(2,1,1)
        hold on
        hanni(1) = plot(dates,FE(:,n), 'r--', 'linewidth', 1);
        
        hanni(2) = plot(dates,  1 * feRMSE(:,n),'k-', 'linewidth', 3);
        plot(dates, -1 * feRMSE(:,n),'k-', 'linewidth', 3)
        %             plot(dates,  2 * feRMSE(:,n),'k-', 'linewidth', 1);
        %             plot(dates, -2 * feRMSE(:,n),'k-', 'linewidth', 1)
        hanni(3) = plot(dates,  1 * RTFEsv(:,n),'b-', 'linewidth', 3);
        plot(dates, -1 * RTFEsv(:,n),'b-', 'linewidth', 3)
        
        plothorzline(0, [], 'k:')
        % title(sprintf('FEh=%d: QRT Coverage Rate 70: %6.4f (p: %4.2f); 95: %6.4f (p: %4.2f)', horizons(n), coverageMean(n,1), coveragePval(n,1), coverageMean(n,2), coveragePval(n,2)))
        title(sprintf('FEh (h=%d) -- Coverage Rate 70, VARSV/RT: %5.2f/%5.2f (p: %4.2f/%4.2f); 95: %5.2f/%5.2f (p: %4.2f/%4.2f)', ...
            horizons(n), coverageMean(n,1), RTcoverageMean(n,1), ...
            coveragePval(n,1), ...
            RTcoveragePval(n,1), ...
            coverageMean(n,2), RTcoverageMean(n,2), ...
            coveragePval(n,2), ...
            RTcoveragePval(n,2)))
        ylimmer = ylim;
        xtickdates(dates(evalTstart-1:end))
        ylim(ylimmer)
        
        %             legend(hanni, 'FE', 'QRT','FINAL')
        wrapcf(sprintf('fe7095SV%d%s', n, modellabel), wrap)
        
    end
    
    
    
    %% plot feSV distro with abs FE
    % close all
    for n = 1 : Nsurvey
        hanni = NaN(3,1);
        
        figure
        hold on
        barcol = .5 * [1 1 1];
        hanni(1) = bar(dates, abs(FE(:,n)), 'FaceColor', barcol, 'EdgeColor', barcol);
        
        hanni(2) = plot(dates, feRMSE(:,n), 'k-', 'linewidth', 3);
        
        
        hanni(3) = plot(dates, RTFEsv(:,n), 'b-', 'linewidth', 3);
        
        maxy = max(ylim);
        xtickdates(dates(evalTstart-1:end))
        ylim([0 maxy])
        
        legend(hanni, 'abs(FE)', 'VARSV', 'FE-CONST')
        
        % title(sprintf('FE (h=%d): QRT Coverage Rate 70: %6.4f (p: %4.2f); 95: %6.4f (p: %4.2f)', horizons(n), coverageMean(n,1), coveragePval(n,1), coverageMean(n,2), coveragePval(n,2)))
        % title(sprintf('FEh=%d: QRT Coverage Rate 70: %6.4f (p: %4.2f); 95: %6.4f (p: %4.2f)', horizons(n), coverageMean(n,1), coveragePval(n,1), coverageMean(n,2), coveragePval(n,2)))
        title(sprintf('FEh (h=%d) -- Coverage Rate 70, VARSV/RT: %5.2f/%5.2f (p: %4.2f/%4.2f); 95: %5.2f/%5.2f (p: %4.2f/%4.2f)', ...
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
    
    
    %% report maxlambda
    lambdaMean  = squeeze(mean(maxlambdaDraws,2));
    lambdaTails = squeeze(prctile(maxlambdaDraws,[5 25 75 95], 2));
    figure
    hold on
    plot(dates(Tstart:end), lambdaMean(Tstart:end), 'r-', 'linewidth', 4)
    plot(dates(Tstart:end), lambdaTails([1 4],Tstart:end), 'r--', 'linewidth', 1)
    ylim([0 1])
    xtickdates(dates(Tstart:end))
    
    wrapcf(sprintf('maxlambda%s', modellabel), wrap)
    
    %% store results
    diary off
    if  doStore && ~quicky
        if doStoreXXL
            save(fullfile(localstore, modellabel), ...
                '-v7.3');
        end
        clear foo feSVdraws feSVFINALdraws feDraws etaSVdraws etaSVFINALdraws RTFEdraws hVCVdraws
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
