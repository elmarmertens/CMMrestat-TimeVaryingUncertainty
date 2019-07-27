%% compare results based on CR70

% note: CR70 coverage rate regressions are already computed as part of the result files,
%       they are recomputed here for different sample choices

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


initwrap
wrap = diary2wrap(wrap, [], false);

%% define subsamples
sam(1).start = datenum(1969,1,1);
sam(1).stop  = datenum(2018,1,1);
sam(1).label = ''; % full sample

sam(2).start = datenum(1995,1,1);
sam(2).stop  = datenum(2018,1,1);
sam(2).label = 'since1995';


sam(3).start = datenum(1969,1,1);
sam(3).stop  = datenum(2007,10,1);
sam(3).label = 'preGreatRecession';


%% select choice of data

for COMPARISONS = {'SVvsRT', 'VARSVvsRT', 'CONSTvsRT',  ...
        'JOINT3SVvsRT','JOINT3SVsinglefactorvsRT', ...
        'FEVARSVlag2vsRT','FEVARSVlag5vsRT', ...
        'SVsinglefactorvsRT', 'SVAR1vsRT'}
    
    comparisonType = COMPARISONS{:};
    jointFlag = false; % default, unless overriden by one of the switch clauses below
    
    
    %% loop over samples
    for s = 1 : length(sam)
        
        % default directories
        datadir        = '../resultfiles/';
        datadirFECONST = '../resultfiles/';
        
        % default varlists
        varlist   = {'eta', 'FE', 'T', 'Nsurvey', ...
            'dates', 'horizons',  ...
            'feSV', ...
            'coverageTstat', 'coveragePval', 'coverageMean'};
        
        varlistRT   = {'dates', ...
            'RTFEsv', ...
            'RTcoverageTstat', 'RTcoveragePval', 'RTcoverageMean'};
        
        switch comparisonType
            case 'CONSTvsRT'
                DATALABELS = {'RGDP', 'UNRATE', 'PGDP', 'CPI', 'TBILL'};
                DATAPRETTYLABELS = {'RGDP', 'UNRATE', 'PGDP', 'CPI', 'TBILL'};
                baselineType  = 'const';
                baselineLabel = 'ETA-SIMPLE';
            case 'SVvsRT'
                DATALABELS = {'RGDP', 'UNRATE', 'PGDP', 'CPI', 'TBILL'};
                DATAPRETTYLABELS = {'RGDP', 'UNRATE', 'PGDP', 'CPI', 'TBILL'};
                baselineType  = 'densitySV';
                baselineLabel = 'SV';
            case 'JOINT3SVvsRT'
                DATALABELS = {'RGDP', 'UNRATE', 'PGDP'};
                DATAPRETTYLABELS = {'RGDP', 'UNRATE', 'PGDP'};
                baselineType  = 'densitySV';
                baselineLabel = 'JOINT-SV';
                jointFlag = true;
                jointLabel = 'UNRATERGDPPGDPmodel';
                jointOrder = [2 1 3];
            case 'JOINT3SVsinglefactorvsRT'
                DATALABELS = {'RGDP', 'UNRATE', 'PGDP'};
                DATAPRETTYLABELS = {'RGDP', 'UNRATE', 'PGDP'};
                baselineType  = 'densitySVsinglefactor';
                baselineLabel = 'JOINT-SV-Factor';
                jointFlag = true;
                jointLabel = 'UNRATERGDPPGDPmodel';
                jointOrder = [2 1 3];
            case 'SVsinglefactorvsRT'
                DATALABELS = {'RGDP', 'UNRATE', 'PGDP', 'CPI', 'TBILL'};
                DATAPRETTYLABELS = {'RGDP', 'UNRATE', 'PGDP', 'CPI', 'TBILL'};
                baselineType  = 'densitySVsinglefactor';
                baselineLabel = 'SV-Factor';
            case 'SVAR1vsRT'
                DATALABELS = {'RGDP', 'UNRATE', 'PGDP', 'CPI', 'TBILL'};
                DATAPRETTYLABELS = {'RGDP', 'UNRATE', 'PGDP', 'CPI', 'TBILL'};
                baselineType  = 'densitySVAR1';
                baselineLabel = 'SV-AR(1)';
            case 'VARSVvsRT'
                DATALABELS = {'RGDP', 'UNRATE', 'PGDP', 'CPI', 'TBILL'};
                DATAPRETTYLABELS = {'RGDP', 'UNRATE', 'PGDP', 'CPI', 'TBILL'};
                baselineType  = 'densityVARSV';
                baselineLabel = 'VAR-SV';
                varlist   = {'eta', 'FE', 'T', 'Nsurvey', ...
                    'dates', 'horizons',  ...
                    'feRMSE', ...
                    'coverageTstat', 'coveragePval', 'coverageMean'};
            case {'FEVARSVlag2vsRT', 'FEVARSVlag5vsRT'}
                DATALABELS = {'RGDP', 'UNRATE', 'PGDP', 'CPI', 'TBILL'};
                DATAPRETTYLABELS = {'RGDP', 'UNRATE', 'PGDP', 'CPI', 'TBILL'};
                baselineType  = 'densityFEVARSV';
                if strcmpi(comparisonType, 'FEVARSVlag2vsRT')
                    baselineLabel = 'FE-VAR(2)-SV';
                else
                    baselineLabel = 'FE-VAR(5)-SV';
                end
                varlist   = {'FE', 'T', 'Nsurvey', ...
                    'dates', 'horizons',  ...
                    'feRMSE', ...
                    'coverageTstat', 'coveragePval', 'coverageMean'};
                
            otherwise
                error('comparisonType %s not supported', comparisonType)
        end
        
        
        
        dummy    = cat(1, varlist, cell(1,length(varlist)));
        models   = repmat(struct(dummy{:}), length(DATALABELS), 1);
        dummy    = cat(1, varlistRT, cell(1,length(varlistRT)));
        modelsRT = repmat(struct(dummy{:}), length(DATALABELS), 1);
        
        
        
        TSTART = 60; % [41 61 81] % [21 41 61 81]
        if TSTART ~= 60 && strcmp(comparisonType, 'VARSVvsRT')
            continue
        end
        EVALTSTART = TSTART;
        
        for d = 1 : length(DATALABELS)
            
            datalabel = DATALABELS{d};
            
            modellabel = datalabel;
            %#ok<*UNRCH>
            
            modellabel = strcat(modellabel, '01');
            
            
            
            modellabelRT   = sprintf('%sQRTfeconst', modellabel);
            modellabelSV   = sprintf('%sQRT%s', modellabel, baselineType);
            
            
            modellabelRT = strcat(modellabelRT, sprintf('TSTART%dEVALSTART%d', TSTART, EVALTSTART));
            modellabelSV = strcat(modellabelSV, sprintf('TSTART%dEVALSTART%d', TSTART, EVALTSTART));
            
            if strcmpi(comparisonType, 'FEVARSVlag2vsRT')
                modellabelSV = strcat(modellabelSV, 'lags2');
            end
            if strcmpi(comparisonType, 'FEVARSVlag5vsRT')
                modellabelSV = strcat(modellabelSV, 'lags5');
            end
            
            %% load results
            if jointFlag
                matfilename = sprintf('slimjoint%s01QRT%sTSTART%dEVALSTART%d.mat', jointLabel, baselineType, TSTART, EVALTSTART);
                models(d)   = load(fullfile(datadir, matfilename), varlist{:});
                models(d).FE    = squeeze(models(d).FE(:,:,jointOrder(d)));
                models(d).feSV  = squeeze(models(d).feSV(:,:,jointOrder(d)));
                
                models(d).coverageMean  = squeeze(models(d).coverageMean(:,jointOrder(d),:));
                models(d).coverageTstat = squeeze(models(d).coverageTstat(:,jointOrder(d),:));
                models(d).coveragePval  = squeeze(models(d).coveragePval(:,jointOrder(d),:));
            else
                matfilename = strcat('slim', modellabelSV, '.mat');
                models(d)   = load(fullfile(datadir, matfilename), varlist{:});
            end
            
            matfilename  = strcat('slim', modellabelRT, '.mat');
            modelsRT(d)  = load(fullfile(datadirFECONST, matfilename), varlistRT{:});
            if ~isequal(modelsRT(d).dates, models(d).dates)
                error('date mismatch')
            end
            
            
            %% set sample
            % set FE outside sample to NaN
            % (which will then also limit the coverage regressions)
            sammy = models(d).dates >= sam(s).start & models(d).dates <= sam(s).stop;
            
            models(d).FE(~sammy,:)   = NaN;
 
            %% re-run coverage rate regressions
            Nsurvey     = models(1).Nsurvey;
            horizons    = models(1).horizons;
            switch comparisonType
                case {'VARSVvsRT', 'FEVARSVlag2vsRT', 'FEVARSVlag5vsRT'}
                    feSV = models(d).feRMSE;
                otherwise
                    feSV = models(d).feSV;
            end
            
            %% RECOMPUTE SV coverage rates
            volScales = 1:2;
            [coverageTstat, coveragePval, coverageMean] = deal(NaN(Nsurvey,length(volScales)));
            for v = 1 : length(volScales)
                coverageNdx = double(abs(models(d).FE) < volScales(v) * feSV);
                coverageNdx(isnan(feSV)) = NaN;
                coverageNdx(isnan(models(d).FE))   = NaN;
                
                threshold   = normcdf(volScales(v)) - normcdf(-volScales(v));
                
                coverageMean(:,v) = nanmean(coverageNdx);
                
                for n = 1 : Nsurvey
                    coverrate = coverageNdx(:,n) - threshold;
                    nanny = ~isnan(coverrate);
                    reggae = nwest(coverrate(nanny), ones(sum(nanny),1), ceil(1.5 + horizons(n)));
                    coverageTstat(n,v) = reggae.tstat;
                    coveragePval(n,v)  = tdis_prb(reggae.tstat,reggae.nobs-reggae.nvar);
                end
            end
            
            models(d).coverageMean  = coverageMean;
            models(d).coverageTstat = coverageTstat;
            models(d).coveragePval  = coveragePval;
            
            %% RECOMPUTE RT coverage rates
            volScales = 1:2;
            [RTcoverageTstat, RTcoveragePval, RTcoverageMean] = deal(NaN(Nsurvey,length(volScales)));
            for v = 1 : length(volScales)
                coverageNdx = double(abs(models(d).FE) < volScales(v) * modelsRT(d).RTFEsv);
                coverageNdx(isnan(models(d).FE))   = NaN;
                
                
                threshold   = normcdf(volScales(v)) - normcdf(-volScales(v));
                
                RTcoverageMean(:,v) = nanmean(coverageNdx);
                
                for n = 1 : Nsurvey
                    coverrate = coverageNdx(:,n) - threshold;
                    nanny = ~isnan(coverrate);
                    reggae = nwest(coverrate(nanny), ones(sum(nanny),1), ceil(1.5 + horizons(n)));
                    RTcoverageTstat(n,v) = reggae.tstat;
                    RTcoveragePval(n,v)  = tdis_prb(reggae.tstat,reggae.nobs-reggae.nvar);
                end
            end
            
            
            modelsRT(d).RTcoverageMean  = RTcoverageMean;
            modelsRT(d).RTcoverageTstat = RTcoverageTstat;
            modelsRT(d).RTcoveragePval  = RTcoveragePval;
            
        end % datalabel
        
        %% collect a few parameters
        % (Assuming model(1) is representative)
        Nsurvey     = models(1).Nsurvey;
        horizons    = models(1).horizons;
        
        
        %% CR70 comparison table
        [tableDataRT, tableDataRTpvalues, tableDataSV, tableDataSVpvalues] = deal(NaN(length(DATALABELS), Nsurvey + 2));
        for m = 1 : length(DATALABELS)
            
            
            tableDataRT(m,1:Nsurvey)         = modelsRT(m).RTcoverageMean(:,1) * 100;
            tableDataRTpvalues(m,1:Nsurvey)  = modelsRT(m).RTcoverageTstat(:,1);
            
            tableDataSV(m,1:Nsurvey)         = models(m).coverageMean(:,1) * 100;
            tableDataSVpvalues(m,1:Nsurvey)  = models(m).coverageTstat(:,1);
            
            
            ndx = find(any(~isnan(models(m).FE),2), 1);
            %             if ~isequal(ndx, EVALTSTART)
            %                 warning('Evaluation window seems to have started at t=%d, notEVALTSTART=%d', ndx, EVALTSTART)
            %             end
            tableDataRT(m,end-1) = models(m).dates(ndx);
            
            ndx = find(any(~isnan(models(m).FE),2), 1, 'last');
            tableDataRT(m,end) = models(m).dates(ndx);
        end
        
        filename = sprintf('%stableCR70TSTART%dEVAL%d%s.tex', comparisonType, TSTART, EVALTSTART, sam(s).label);
        
        
        fid = fopen(filename, 'wt');
        % header
        fprintf(fid, '\\begin{normalsize}\n\\begin{center}\n');
        fprintf(fid, '\\begin{tabular}{l%sll}\n', repmat('.4', 1, Nsurvey)');
        fprintf(fid, '\\toprule\n');
        
        % list horizons
        fprintf(fid, '  & \\multicolumn{%d}{c}{Forecast horizon}', Nsurvey);
        fprintf(fid, '  & \\multicolumn{2}{c}{eval. window}');
        fprintf(fid, '\\\\\n');
        fprintf(fid, '\\cmidrule(r){%d-%d}', 1+1,1+Nsurvey);
        fprintf(fid, '\\cmidrule(l){%d-%d}', 2+Nsurvey,3+Nsurvey);
        fprintf(fid, '\n');
        
        fprintf(fid, 'Variable ');
        fprintf(fid, ' & \\ccol{%d} ', horizons);
        fprintf(fid, ' & \\ccol{eval. begin} ');
        fprintf(fid, ' & \\ccol{eval. end} ');
        fprintf(fid, '\\\\\n'); % new line
        
        % PANEL A: SV
        fprintf(fid, '\\midrule\n');
        fprintf(fid, '\\multicolumn{%d}{c}{\\textbf{Panel A: %s}}\\\\\n', Nsurvey + 2, baselineLabel);
        fprintf(fid, '\\midrule\n');
        for m = 1 : size(tableDataRT,1)
            
            fprintf(fid, '%s  ', DATAPRETTYLABELS{m});
            
            for n = 1 : Nsurvey
                fprintf(fid, '& %5.2f%s ', tableDataSV(m,n), Zstar(tableDataSVpvalues(m,n)));
            end
            fprintf(fid, '& %s ', datestr(tableDataRT(m,end-1), 'YYYY:QQ'));
            fprintf(fid, '& %s ', datestr(tableDataRT(m,end), 'YYYY:QQ'));
            fprintf(fid, '\\\\\n');
        end
        % PANEL B: CONST
        fprintf(fid, '\\midrule\n');
        fprintf(fid, '\\multicolumn{%d}{c}{\\textbf{Panel B: FE-SIMPLE}}\\\\\n', Nsurvey + 2);
        fprintf(fid, '\\midrule\n');
        for m = 1 : size(tableDataRT,1)
            
            fprintf(fid, '%s  ', DATAPRETTYLABELS{m});
            
            for n = 1 : Nsurvey
                fprintf(fid, '& %5.2f%s ', tableDataRT(m,n), Zstar(tableDataRTpvalues(m,n)));
            end
            fprintf(fid, '& %s ', datestr(tableDataRT(m,end-1), 'YYYY:QQ'));
            fprintf(fid, '& %s ', datestr(tableDataRT(m,end), 'YYYY:QQ'));
            
            fprintf(fid, '\\\\\n');
        end
        
        
        % footer
        % evalStart = modelsCONST(m).dates(TSTART);
        evalStop  = models(m).dates(end);
        fprintf(fid, '\\bottomrule\n');
        fprintf(fid, '\\end{tabular}\n');
        fprintf(fid, '\\end{center}\n');
        
        % baseline labels
        H0model   = 'FE-SIMPLE';
        % adjust as needed
        switch comparisonType
            case {'SVvsRT', 'spfSVvsRT'}
                thisModel = 'SV';
            case {'VARSVvsRT', 'spfVARSVvsRT'}
                thisModel = 'VAR-SV';
            case 'CONSTvsRT'
                thisModel = 'ETA-SIMPLE';
            case 'JOINT3SVvsRT'
                thisModel = 'JOINT-SV';
            case 'JOINT3SVsinglefactorvsRT'
                thisModel = 'JOINT-SV-Factor';
            case 'SVsinglefactorvsRT'
                thisModel = 'SV-Factor';
            case 'SVAR1vsRT'
                thisModel = 'SV-AR(1)';
            case 'FEVARSVlag2vsRT'
                thisModel = 'FE-VAR(2)-SV';
            case 'FEVARSVlag5vsRT'
                thisModel = 'FE-VAR(5)-SV';
            otherwise
                error('comparisonType %s not defined', comparisonType)
        end
        
        fprintf(fid, 'Note: The table reports the empirical out-of-sample coverage rates of one-standard-deviation bands.\n');
        if any(contains(DATALABELS, 'Greenbook'))
            fprintf(fid, 'The sample uses predictions made from the date given in the second right-most column through the date in the right-most column (and realized forecast errors as far as available), in case of the SPF, and through 2011:Q4, in case of the Greenbook  (evaluated against realized data as far as available in both cases).\n');
        else
            fprintf(fid, 'The sample uses predictions made from the date given in the second right-most column through the date in the right-most column (and realized forecast errors as far as available).\n');
        end
        fprintf(fid, 'The upper panel provides results based on our proposed multi-horizon %s model.\n', thisModel);
        fprintf(fid, 'The lower panel provides results based on the %s model estimated over rolling windows with %d~quarterly observations.\n', H0model, TSTART);
        fprintf(fid, 'Statistically significant departures from a nominal coverage of 68\\%% (as predicted under a normal distribution) are indicated by *, **, or ***, corresponding to 10, 5, and 1 percent significance, respectively.\n');
        fprintf(fid, '\\end{normalsize}\n');
        
        fclose(fid);
        type(filename)
        display(filename)
        
        
        if ~isempty(wrap)
            movefile(filename, wrap.dir)
            latexwrapper(wrap, 'add', 'table', filename);
        end
        
        
    end % sample
end % comparisonType



%% finish / clean up
finishwrap
finishscript
dockAllFigures

