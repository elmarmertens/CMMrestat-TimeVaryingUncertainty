%% JOINT ETA-SV model across multiple variables
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

Nstreams    = max(getparpoolsize, 1);
rndStreams  = initRandStreams(Nstreams, [], 0);


quicky     = false; % if TRUE: very short MCMC chains, no looping across variables,
%  useful for testing code execution, see below for specific settings
doStore    = true; % to store result files in directory specified by "localstore" m-file
doStoreXXL = false; % to store extended result files (with MCMC draws)
DATALABELS = {'UNRATE', 'RGDP', 'PGDP'};


Tstart     = 60;
RTwindow   = Tstart;
evalTstart = Tstart;

doDiffuseSlopes = true;
doZeroSlopes = false;


if quicky
    MCMCburnin   = 1e2;
    MCMCdraws    = 1e2;
    Nfedraws     = 5e1;
else
    MCMCburnin   = 1e5;
    MCMCdraws    = 1e3;
    Nfedraws     = 100;
end

% construct a model-label indicating dataset and important parameters, to
% append to picture names
modellabel = sprintf('joint%smodel', sprintf('%s', DATALABELS{:}));
%#ok<*UNRCH>
if doZeroSlopes
    modellabel = strcat(modellabel, '00');
else
    modellabel = strcat(modellabel, '01');
end

modellabel = strcat(modellabel,'QRTdensitySV');

modellabel = strcat(modellabel, sprintf('TSTART%dEVALSTART%d', Tstart, evalTstart));


wrap = [];
titlename = modellabel;
if ~quicky
    initwrap
end

fprintf('Processing %s ... \n', modellabel)

%% load data
Ny = length(DATALABELS);
Nsurvey = 5;
dates = genrQdates(1969,2018);
dates = dates(1:end-3);
T     = length(dates);

ETA         = NaN(T,Nsurvey,Ny);
FE          = NaN(T,Nsurvey,Ny);
FEobserved  = NaN(T,Nsurvey,Ny);
ETAnanny    = true(T, Nsurvey, Ny);

horizons = 0:4;

for d = 1 : Ny
    
    datalabel = DATALABELS{d};
    matfile = sprintf('hydepark%sdata', upper(datalabel));
    thisdata = load(matfile, 'Yforecast', 'Yrealized', 'eta', 'FE', 'FEobserved', 'T', 'Nsurvey', 'dates', 'horizons', 'maxHorizon', 'Ylabel');
    
    
    
    %% record location of missing values (NaN)
    eta        = thisdata.eta;
    etaNanny   = isnan(eta);
    fe         = thisdata.FE;
    feObserved = thisdata.FEobserved;
    
    ndx2all  = ismember(thisdata.dates, dates);
    ndx2this = ismember(dates, thisdata.dates);
    
    if thisdata.dates(end)  > dates(end)
        warning('individual data samples extend beyond common sample [ PRESS ANY KEY ]')
        pause
    end
    
    ETA(ndx2this,:,d)         = eta(ndx2all,:);
    ETAnanny(ndx2this,:,d)    = etaNanny(ndx2all,:);
    FE(ndx2this,:,d)          = fe(ndx2all,:);
    FEobserved(ndx2this,:,d)  = feObserved(ndx2all,:);
    
end

close all

%% compute forecast error SV
% allocate memory
feSVdraws       = NaN(T,Nsurvey,Ny,MCMCdraws);
etaSVdraws      = NaN(T,Nsurvey,Ny,MCMCdraws);
feDraws         = NaN(T,Nsurvey,Ny,MCMCdraws*Nfedraws);
feVol           = NaN(T,Nsurvey,Ny); % to store SE of feDraws

hVCVdraws         = NaN(Nsurvey*Ny,Nsurvey*Ny,MCMCdraws,T);

parfor thisT = Tstart : T - 1 % @parfor
    TID = parid;
    [mcmc_feSV, mcmc_fe, mcmc_etaSV, mcmc_hVCV] = mcmcsamplerSVjoint(ETA, ETAnanny, ...
        thisT, horizons, MCMCdraws, MCMCburnin, Nfedraws, doZeroSlopes, doDiffuseSlopes, rndStreams{TID}, TID == 1);  %#ok<*PFBNS>
    feSVdraws(thisT,:,:,:)             = mcmc_feSV(end,:,:,:);
    etaSVdraws(thisT,:,:,:)            = mcmc_etaSV(end,:,:,:);
    feDraws(thisT,:,:,:)               = mcmc_fe;
    feVol(thisT,:,:)                   = std(mcmc_fe,0,3);
    hVCVdraws(:,:,:,thisT)             = mcmc_hVCV;
    fprintf('%s - QRT: done with t=%d\n', modellabel, thisT)
end

TID = 1;
% doing the last step outside of the parfor loop (to grab final
% draws)
thisT                            = T;
fprintf('%s - QRT proceeding to T=%d\n', modellabel, thisT)

[mcmc_feSV, mcmc_fe, mcmc_etaSV, mcmc_hVCV] = mcmcsamplerSVjoint(ETA, ETAnanny, ...
    thisT, horizons, MCMCdraws, MCMCburnin, Nfedraws, doZeroSlopes, doDiffuseSlopes, rndStreams{TID}, true);
feSVdraws(thisT,:,:,:)             = mcmc_feSV(end,:,:,:);
etaSVdraws(thisT,:,:,:)            = mcmc_etaSV(end,:,:,:);
feDraws(thisT,:,:,:)               = mcmc_fe;
feVol(thisT,:,:)                   = std(mcmc_fe,0,3);
hVCVdraws(:,:,:,thisT)           = mcmc_hVCV;
fprintf('%s - QRT: done with t=%d\n', modellabel, thisT)


% [feSVFINALdraws, etaSVFINALdraws] = deal(NaN(T,Nsurvey,Ny,MCMCdraws));
feSVFINALdraws         = mcmc_feSV;
etaSVFINALdraws        = mcmc_etaSV;


feSV      = sqrt(mean(feSVdraws,4));
feSVmed   = sqrt(median(feSVdraws,4));
feSVtails = sqrt(prctile(feSVdraws, [25 75], 4));


feSVFINAL      = squeeze(sqrt(mean(feSVFINALdraws,4)));
feSVFINALmed   = squeeze(sqrt(median(feSVFINALdraws,4)));
feSVFINALtails = squeeze(sqrt(prctile(feSVFINALdraws, [25 75], 4)));

etaSV      = sqrt(mean(etaSVdraws,4));
etaSVmed   = sqrt(median(etaSVdraws,4));
etaSVtails = sqrt(prctile(etaSVdraws, [25 75], 4));


etaSVFINAL      = squeeze(sqrt(mean(etaSVFINALdraws,4)));
etaSVFINALmed   = squeeze(sqrt(median(etaSVFINALdraws,4)));
etaSVFINALtails = squeeze(sqrt(prctile(etaSVFINALdraws, [25 75], 4)));

feMean = mean(feDraws, 4);
feTails = prctile(feDraws, [.05 2.5 5 25 75 95 97.5 99.5], 4);

clear mcmc_*
clear etaSVFINALdraws* feSVFINALdraws*
clear etaSVdraws* feSVdraws*

%% set FE over estimation window to NaN
FE(1:Tstart-1,:)      = NaN;
eta(1:Tstart-1,:)     = NaN;
% evalTstart should always be larger than Tstart, thus making the previous two lines redundant
FE(1:evalTstart-1,:)  = NaN;
eta(1:evalTstart-1,:) = NaN;

%% document correlation between SV shocks
Nsv = Nsurvey*Ny;
hVCVeigen = NaN(Nsv,MCMCdraws,T);
for thisT = Tstart+ 1 : T
    for n = 1 : MCMCdraws
        hVCVeigen(:,n,thisT) = sort(eig(hVCVdraws(:,:,n,thisT)), 'descend');
    end
end

tails = prctile(hVCVeigen, [5 95], 2);
figure
for n = 1 : Nsv
    hold on
    plot(dates, squeeze(mean(hVCVeigen(n,:,:), 2)), '-', 'color', Colors4Plots(n), 'linewidth', 3)
    plot(dates, permute(tails(n,:,:), [3 2 1]), '-', 'color', Colors4Plots(n))
    ylim([0 max(ylim)])
    xtickdates(dates([Tstart T]))
end
wrapcf(sprintf('hVCVeigenvalues%s', modellabel), wrap)


%% compute coverage rates
volScales = 1:2;
[coverageTstat, coveragePval, coverageMean] = deal(NaN(Nsurvey,Ny,length(volScales)));
for v = 1 : length(volScales)
    coverageNdx = double(abs(FE) < volScales(v) * feSV);
    coverageNdx(isnan(feSV)) = NaN;
    coverageNdx(isnan(FE))   = NaN;
    
    threshold   = normcdf(volScales(v)) - normcdf(-volScales(v));
    
    coverageMean(:,:,v) = nanmean(coverageNdx);
    
    for m = 1 : Ny
        for n = 1 : Nsurvey
            coverrate = coverageNdx(:,n,m) - threshold;
            nanny = ~isnan(coverrate);
            reggae = nwest(coverrate(nanny), ones(sum(nanny),1), ceil(1.5 + horizons(n)));
            %                                 hrulefill
            %                                 fprintf('Horizon %d, Coverage rate for %d times SV is %6.4f \n', horizons(n), volScales(v), coverageMean(n,v))
            %                                 prt(reggae)
            coverageTstat(n,m,v) = reggae.tstat;
            coveragePval(n,m,v)  = tdis_prb(reggae.tstat,reggae.nobs-reggae.nvar);
        end
    end
end

% same for eta
[ETAcoverageTstat, ETAcoveragePval, ETAcoverageMean] = deal(NaN(Nsurvey, Ny,length(volScales)));
for v = 1 : length(volScales)
    ETAcoverageNdx = double(abs(ETA) < volScales(v) * etaSV);
    ETAcoverageNdx(isnan(etaSV)) = NaN;
    ETAcoverageNdx(isnan(eta))   = NaN;
    
    threshold   = normcdf(volScales(v)) - normcdf(-volScales(v));
    
    ETAcoverageMean(:,:,v) = nanmean(ETAcoverageNdx);
    
    for m = 1 : Ny
        for n = 1 : Nsurvey
            coverrate = ETAcoverageNdx(:,n,m) - threshold;
            nanny = ~isnan(coverrate);
            reggae = nwest(coverrate(nanny), ones(sum(nanny),1), ceil(1.5 + horizons(n)));
            %                                 hrulefill
            %                                 fprintf('Horizon %d, ETA-Coverage rate for %d times SV is %6.4f \n', horizons(n), volScales(v), ETAcoverageMean(n,v))
            %                                 prt(reggae)
            ETAcoverageTstat(n,m,v) = reggae.tstat;
            ETAcoveragePval(n,m,v)  = tdis_prb(reggae.tstat,reggae.nobs-reggae.nvar);
        end
    end
end


%% compute quantile-coverage-rates
quantileP = [2.5 5 normcdf(-1) * 100 25 50 75 normcdf(1) * 100 95 97.5];
criticalQuants = prctile(feDraws, quantileP, 4);

fee                       = repmat(FE, [1 1 1 length(quantileP)]);
quantcoverNdx             = double(fee < criticalQuants); % bsxfun would work, but fee needed for next step
quantcoverNdx(isnan(fee) | isnan(criticalQuants)) = NaN;
quantcoverMean            = squeeze(nanmean(quantcoverNdx,1));

[quantcoverTstat, quantcoverPval] =deal(NaN(size(quantcoverMean)));

% test for significance with Newey-West regressions
for q = 1 : length(quantileP)
    
    for m = 1 : Ny
        for n = 1 : Nsurvey
            coverrate = quantcoverNdx(:,n,m,q) - quantileP(q) / 100;
            nanny = ~isnan(coverrate);
            reggae = nwest(coverrate(nanny), ones(sum(nanny),1), ceil(1.5 + horizons(n)));
            %                                 hrulefill
            %                                 fprintf('Horizon %d, Quantile-Coverage rate for %4.2f%% is %4.2f%% \n', horizons(n), quantileP(q), quantcoverMean(n,q) * 100);
            %                                 prt(reggae)
            quantcoverTstat(n,m,q) = reggae.tstat;
            quantcoverPval(n,m,q)  = tdis_prb(reggae.tstat,reggae.nobs-reggae.nvar);
        end
    end
end

%% compute CRPS
crps     = NaN(T,Nsurvey,Ny);
feDraws = permute(feDraws, [4 1 2 3]); % should make the following loop quicker
for m = 1 : Ny
    for n = 1 : Nsurvey
        parfor t = 1 : T % could just start at Tstart
            crps(t,n,m) = crpsDraws(FE(t,n,m), feDraws(:,t,n,m));
        end % t
    end % n
end
clear feDraws

%% compute RT SV bands
RTFEsv     = NaN(T,Nsurvey,Ny);

for n = 1 : Ny
    for thisT = Tstart : T % @parfor
        theseFE           = FEobserved(thisT-RTwindow+1:thisT,:,n);
        RTFEsv(thisT,:,n)   = sqrt(nanmean(theseFE.^2));
        %                 fprintf('%s - RT : done with t=%d\n', modellabel, thisT)
    end
end

RTFEdraws = randn(T,Nsurvey,Ny,MCMCdraws * Nfedraws);
RTFEdraws = bsxfun(@times, RTFEdraws, RTFEsv);


%% compute RT coverage rates
volScales = 1:2;
[RTcoverageTstat, RTcoveragePval, RTcoverageMean] = deal(NaN(Nsurvey,Ny,length(volScales)));
for v = 1 : length(volScales)
    coverageNdx = double(abs(FE) < volScales(v) * RTFEsv);
    coverageNdx(isnan(feSV)) = NaN;
    coverageNdx(isnan(FE))   = NaN;
    
    threshold   = normcdf(volScales(v)) - normcdf(-volScales(v));
    
    RTcoverageMean(:,:,v) = nanmean(coverageNdx);
    
    for m = 1 : Ny
        for n = 1 : Nsurvey
            coverrate = coverageNdx(:,n,m) - threshold;
            nanny = ~isnan(coverrate);
            reggae = nwest(coverrate(nanny), ones(sum(nanny),1), ceil(1.5 + horizons(n)));
            %                                 hrulefill
            %                                 fprintf('Horizon %d, Coverage rate for %d times SV is %6.4f \n', horizons(n), volScales(v), RTcoverageMean(n,v))
            %                                 prt(reggae)
            RTcoverageTstat(n,m,v) = reggae.tstat;
            RTcoveragePval(n,m,v)  = tdis_prb(reggae.tstat,reggae.nobs-reggae.nvar);
        end
    end
end

%% RT - quantiles
quantileP = [2.5 5 normcdf(-1) * 100 25 50 75 normcdf(1) * 100 95 97.5];
criticalQuants = prctile(RTFEdraws, quantileP, 4);

fee                       = repmat(FE, [1 1 1 length(quantileP)]);
quantcoverNdx             = double(fee < criticalQuants); % bsxfun would work, but fee needed for next step
quantcoverNdx(isnan(fee) | isnan(criticalQuants)) = NaN;

RTquantcoverMean          = squeeze(nanmean(quantcoverNdx,1));

[RTquantcoverTstat, RTquantcoverPval] =deal(NaN(size(RTquantcoverMean)));

% test for significance with Newey-West regressions
for q = 1 : length(quantileP)
    
    for m = 1 : Ny
        for n = 1 : Nsurvey
            coverrate = quantcoverNdx(:,n,m,q) - quantileP(q) / 100;
            nanny = ~isnan(coverrate);
            reggae = nwest(coverrate(nanny), ones(sum(nanny),1), ceil(1.5 + horizons(n)));
            %                                 hrulefill
            %                                 fprintf('Horizon %d, Quantile-Coverage rate for %4.2f%% is %4.2f%% \n', horizons(n), quantileP(q), RTquantcoverMean(n,q) * 100);
            %                                 prt(reggae)
            RTquantcoverTstat(n,m,q) = reggae.tstat;
            RTquantcoverPval(n,m,q)  = tdis_prb(reggae.tstat,reggae.nobs-reggae.nvar);
        end
    end
end



%% RT - CRPS
RTcrps    = NaN(T,Nsurvey,Ny);
RTFEdraws = permute(RTFEdraws, [4 1 2 3]); % should make the following loop quicker
% progressbar(0)
for m = 1 : Ny
    for n = 1 : Nsurvey
        parfor t = 1 : T % could simply start at Tstart
            RTcrps(t,n,m) = crpsDraws(FE(t,n,m), RTFEdraws(:,t,n,m));
        end % t
    end % n
end % m


if ~isequal(isnan(RTcrps), isnan(FE))
    error('crps has NaN when there is data')
end

%% plot etaSV distro with abs eta
closeall
for m = 1 : Ny
    datalabel = DATALABELS{m};
    for n = 1 : Nsurvey
        hanni = NaN(3,1);
        
        figure
        hold on
        barcol = .3 * [0 1 0];
        hanni(1) = bar(dates, abs(ETA(:,n,m)), 'FaceColor', barcol, 'EdgeColor', barcol);
        
        hanni(2) = plot(dates, etaSV(:,n,m), 'k-', 'linewidth', 3);
        plot(dates, squeeze(etaSVtails(:,n,m,:)), 'k-', 'linewidth', 1)
        
        hanni(3) = plot(dates, etaSVFINAL(:,n,m), 'r--', 'linewidth', 3);
        plot(dates, squeeze(etaSVFINALtails(:,n,m,:)), 'r--', 'linewidth', 1)
        
        maxy = max(ylim);
        xtickdates(dates(evalTstart-1:end))
        ylim([0 maxy])
        
        legend(hanni, 'abs(eta)', 'QRT', 'Final')
        
        title(sprintf('%s ETA (h=%d): QRT Coverage Rate 70: %6.4f (p: %4.2f); 95: %6.4f (p: %4.2f)', ...
            datalabel, horizons(n), ETAcoverageMean(n,1), coveragePval(n,1), coverageMean(n,2), coveragePval(n,2)))
        wrapcf(sprintf('etaSViqr%d%s%s', n, modellabel, datalabel), wrap)
    end
end

%% plot feSV distro with abs FE
closeall
for m = 1 : Ny
    datalabel = DATALABELS{m};
    for n = 1 : Nsurvey
        hanni = NaN(5,1);
        
        figure
        hold on
        barcol = .5 * [1 1 1];
        hanni(1) = bar(dates, abs(FE(:,n,m)), 'FaceColor', barcol, 'EdgeColor', barcol);
        
        hanni(2) = plot(dates, feSV(:,n,m), 'k-', 'linewidth', 3);
        plot(dates, squeeze(feSVtails(:,n,m,:)), 'k-', 'linewidth', 1)
        
        hanni(3) = plot(dates, feSVFINAL(:,n,m), 'r--', 'linewidth', 3);
        plot(dates, squeeze(feSVFINALtails(:,n,m,:)), 'r--', 'linewidth', 1)
        
        hanni(4) = plot(dates, feVol(:,n,m), 'y--', 'linewidth', 2);
        
        hanni(5) = plot(dates, RTFEsv(:,n,m), 'b-', 'linewidth', 3);
        
        maxy = max(ylim);
        xtickdates(dates(evalTstart-1:end))
        ylim([0 maxy])
        
        legend(hanni, 'abs(FE)', 'QRT', 'Final', 'QRT (montecarlo)', 'Reifschneider-Tulip')
        
        title(sprintf('FE (h=%d): QRT Coverage Rate 70: %6.4f (p: %4.2f); 95: %6.4f (p: %4.2f)', horizons(n), coverageMean(n,1), coveragePval(n,1), coverageMean(n,2), coveragePval(n,2)))
        % title(sprintf('FEh=%d: QRT Coverage Rate 70: %6.4f (p: %4.2f); 95: %6.4f (p: %4.2f)', horizons(n), coverageMean(n,1), coveragePval(n,1), coverageMean(n,2), coveragePval(n,2)))
        title(sprintf('FEh (h=%d) -- Coverage Rate 70, SV/RT: %6.4f/%6.4f (p: %4.2f/%4.2f); 95: %6.4f/%6.4f (p: %4.2f/%4.2f)', ...
            horizons(n), coverageMean(n,1), RTcoverageMean(n,1), ...
            coveragePval(n,1), ...
            RTcoveragePval(n,1), ...
            coverageMean(n,2), RTcoverageMean(n,2), ...
            coveragePval(n,2), ...
            RTcoveragePval(n,2)))
        wrapcf(sprintf('feSViqr%d%s%s', n, modellabel, datalabel), wrap)
    end
end

%% compare CRPS
[avgCRPS, RTavgCRPS] = deal(NaN(T,Nsurvey,Ny));
for t = 1 : T
    avgCRPS(t,:,:)   = nanmean(crps(1:t,:,:),1);
    RTavgCRPS(t,:,:) = nanmean(RTcrps(1:t,:,:),1);
end

closeall
for m = 1 : Ny
    datalabel = DATALABELS{m};
    for n = 1 : Nsurvey
        figure
        subplot(2,1,1)
        hanni = NaN(2,1);
        hold on
        hanni(1) = plot(dates, RTavgCRPS(:,n,m), 'k-', 'linewidth', 1);
        hanni(2) = plot(dates, avgCRPS(:,n,m), 'r--', 'linewidth', 2);
        plotOrigin
        xtickdates(dates(evalTstart-1:end))
        legend(hanni, 'RT', 'SV', 'location', 'best')
        title(sprintf('%s -- CRPS (RT - SV) / RT: %5.2f%% ', ...
            datalabel, (RTavgCRPS(end,n,m) - avgCRPS(end,n,m))/RTavgCRPS(end,n,m) * 100))
        
        subplot(2,1,2)
        hanni = NaN(2,1);
        hold on
        hanni(1) = plot(dates, RTcrps(:,n,m), 'k-');
        hanni(2) = plot(dates, crps(:,n,m), 'r--');
        xtickdates(dates(evalTstart-1:end))
        legend(hanni, 'RT', 'SV')
        title('contributions')
        
        wrapcf(sprintf('crps%d%s%s', n, modellabel,datalabel), wrap)
    end
end

%% store results
diary off
if  doStore && ~quicky
    if doStoreXXL
        save(fullfile(localstore, modellabel), ...
            '-v7.3');
    end
    clear feDraws  RTFEdraws hVCVdraws
    save(fullfile(localstore, strcat('slim', modellabel)), ...
        '-v7.3');
end

%% finish / clean up
finishwrap % if latex compilation fails, edit this script to turn the compilation offfinishscript
finishscript
dockAllFigures
