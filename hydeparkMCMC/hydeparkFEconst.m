%% FE-SIMPLE estimates
% ... results for the latter are labeled "RT..." as the underlying computations
% ... intend to mimic the approach of Reifschneider and Tulip (2007, 2017)
% ... the same results are also computed as part of hydeparkETAsv etc;
% ... this scripts serves to store the FE-SIMPLE results separately as well

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

rng(0);
rndStream  = getDefaultStream;

quicky     = false;
doStore    = true;
doStoreXXL = false;
DATALABELS = {'CPI', 'TBILL', 'UNRATE', 'PGDP', 'RGDP'};

%% loop over window-lengths choices and variabels 
for Tstart = [40 60 80]
    RTwindow = Tstart;
    switch Tstart
        case {40,60}
            evalTstart   = 60;
        case 80
            evalTstart   = 80;
    end
    
    doDiffuseSlopes = true;
    doZeroSlopes    = false;
    tic
    
    for d =  1 : length(DATALABELS)
        
        close all
        datalabel = DATALABELS{d};
        
        if quicky
            Ndraws = 1e4;
        else
            Ndraws = 1e5;
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
        
        modellabel = strcat(modellabel, 'QRTfeconst');
        
        
        
        modellabel = strcat(modellabel, sprintf('TSTART%dEVALSTART%d', Tstart, evalTstart));
        
        
        wrap = [];
        
        
        fprintf('Processing %s ... \n', modellabel)
        
        %% load data
        matfile = sprintf('hydepark%sdata', upper(datalabel));
        load(matfile, 'Yforecast', 'Yrealized', 'eta', 'FE', 'FEobserved', 'T', 'Nsurvey', 'dates', 'horizons', 'maxHorizon', 'Ylabel')
        
        
        
        
        %% record location of missing values (NaN)
        etaNanny = isnan(eta);
        
        samStart = find(~(etaNanny(:,1)),1);
        
        eta        = eta(samStart:end, :);
        etaNanny   = etaNanny(samStart:end, :);
        FE         = FE(samStart:end, :);
        FEobserved = FEobserved(samStart:end, :);
        dates      = dates(samStart:end);
        T          = length(dates);
        
        
        
        %% compute RT SV bands
        RTFEsv     = NaN(T,Nsurvey);
        
        parfor thisT = Tstart : T % @parfor
            theseFE           = FEobserved(thisT-RTwindow+1:thisT,:); %#ok<PFBNS>
            RTFEsv(thisT,:)   = sqrt(nanmean(theseFE.^2));
            %                 fprintf('%s - RT : done with t=%d\n', modellabel, thisT)
        end
        
        RTFEdraws = randn(T,Nsurvey,Ndraws);
        RTFEdraws = bsxfun(@times, RTFEdraws, RTFEsv);
        
        %% set FE over evaluation window to NaN
        FE(1:Tstart-1,:)      = NaN;
        eta(1:Tstart-1,:)     = NaN;
        % evalTstart should always be larger than Tstart, thus making he previous two lines redundant
        FE(1:evalTstart-1,:)  = NaN;
        eta(1:evalTstart-1,:) = NaN;
        
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
        RTcrps  = NaN(T,Nsurvey);
        foodraws = permute(RTFEdraws, [3 1 2]); % should make the following loop quicker
        % progressbar(0)
        parfor n = 1 : Nsurvey
            for t = 1 : T % could just start at Tstart
                RTcrps(t,n) = crpsDraws(FE(t,n), foodraws(:,t,n));
            end % t
        end % n
        clear foodraws
        
        %%
        
        if ~isequal(isnan(RTcrps), isnan(FE))
            error('crps has NaN when there is data')
        end
        
        %% plot feSV distro with abs FE
        % close all
        for n = 1 : Nsurvey
            hanni = NaN(5,1);
            
            figure
            hold on
            barcol = .5 * [1 1 1];
            hanni(1) = bar(dates, abs(FE(:,n)), 'FaceColor', barcol, 'EdgeColor', barcol);
            
            hanni(5) = plot(dates, RTFEsv(:,n), 'b-', 'linewidth', 3);
            maxy = max(ylim);
            xtickdates(dates(evalTstart-1:end))
            ylim([0 maxy])
            
            wrapcf(sprintf('feSViqr%d%s', n, modellabel), wrap)
        end
        
        %% plot CRPS
        RTavgCRPS = deal(NaN(T,Nsurvey));
        for t = 1 : T
            RTavgCRPS(t,:) = nanmean(RTcrps(1:t,:));
        end
        
        for n = 1 : Nsurvey
            figure
            subplot(2,1,1)
            %                             hanni = NaN(2,1);
            hold on
            plot(dates, RTavgCRPS(:,n), 'k-', 'linewidth', 1);
            plotOrigin
            xtickdates(dates(evalTstart-1:end))
            % legend(hanni, 'RT', 'CONST', 'location', 'best')
            %                             title(sprintf('CRPS (RT - CONST) / RT: %5.2f%% ', ...
            %                                 (RTavgCRPS(end,n) - avgCRPS(end,n))/RTavgCRPS(end,n) * 100))
            
            subplot(2,1,2)
            hanni = NaN(2,1);
            hold on
            hanni(1) = plot(dates, RTcrps(:,n), 'k-');
            xtickdates(dates(evalTstart-1:end))
            % legend(hanni, 'RT', 'CONST')
            title('contributions')
            
            wrapcf(sprintf('crps%d%s', n, modellabel), wrap)
        end
        
        
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
        finishwrap
        if quicky
            dockAllFigures
            toc
            return
        else
            close all
        end
    end
end

toc

%% finish / clean up
finishscript
dockAllFigures
