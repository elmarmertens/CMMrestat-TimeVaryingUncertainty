%% compare results based on CR70, compared against FE-SIMPLE using 80-quarter window
% thus using also a shorter evaluation window
% note: SV results are obtained from a separate run of the MCMC file starting only as of the 80th quarter (TSTART=80)

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


%% select choice of data

for COMPARISONS = {'SVvsRT'}
    
    comparisonType = COMPARISONS{:};
    jointFlag = false; % default, unless overriden by one of the switch clauses below
    
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
        case 'SVvsRT'
            DATALABELS = {'RGDP', 'UNRATE', 'PGDP', 'CPI', 'TBILL'};
            DATAPRETTYLABELS = {'RGDP', 'UNRATE', 'PGDP', 'CPI', 'TBILL'};
            baselineType  = 'densitySV';
            baselineLabel = 'ETA-SV';
        otherwise
            error('comparisonType %s not supported', comparisonType)
    end
    
    
    
    dummy    = cat(1, varlist, cell(1,length(varlist)));
    models   = repmat(struct(dummy{:}), length(DATALABELS), 1);
    dummy    = cat(1, varlistRT, cell(1,length(varlistRT)));
    modelsRT = repmat(struct(dummy{:}), length(DATALABELS), 1);
    
    
    
    for TSTART = 80
        
        for EVALTSTART = 80
            
            if TSTART > EVALTSTART
                continue
            end
            for d = 1 : length(DATALABELS)
                
                datalabel = DATALABELS{d};
                
                modellabel = datalabel;
                %#ok<*UNRCH>
                
                modellabel = strcat(modellabel, '01');
                
                
                
                modellabelRT   = sprintf('%sQRTfeconst', modellabel);
                modellabelSV   = sprintf('%sQRT%s', modellabel, baselineType);
                
                
                modellabelRT = strcat(modellabelRT, sprintf('TSTART%dEVALSTART%d', TSTART, EVALTSTART));
                modellabelSV = strcat(modellabelSV, sprintf('TSTART%dEVALSTART%d', EVALTSTART, EVALTSTART));
                
                
                %% load results
                matfilename = strcat('slim', modellabelSV, '.mat');
                models(d)   = load(fullfile(datadir, matfilename), varlist{:});
                
                matfilename  = strcat('slim', modellabelRT, '.mat');
                modelsRT(d)  = load(fullfile(datadirFECONST, matfilename), varlistRT{:});
                if ~isequal(modelsRT(d).dates, models(d).dates)
                    error('date mismatch')
                end
                
            end % datalabel
            
            %% collect a few parameters
            % (Assuming model(1) is representative)
            Nsurvey     = models(1).Nsurvey;
            horizons    = models(1).horizons;
            
            
            
            
            %% CR70 comparison table
            [tableDataRT, tableDataRTpvalues, tableDataSV, tableDataSVpvalues] = deal(NaN(length(DATALABELS), Nsurvey + 1));
            for m = 1 : length(DATALABELS)
                
                
                tableDataRT(m,1:Nsurvey)         = modelsRT(m).RTcoverageMean(:,1) * 100;
                tableDataRTpvalues(m,1:Nsurvey)  = modelsRT(m).RTcoverageTstat(:,1);
                
                tableDataSV(m,1:Nsurvey)         = models(m).coverageMean(:,1) * 100;
                tableDataSVpvalues(m,1:Nsurvey)  = models(m).coverageTstat(:,1);
                
                
                ndx = find(any(~isnan(models(m).FE),2), 1);
                if ~isequal(ndx, EVALTSTART)
                    warning('Evaluation window seems to have started at t=%d, notEVALTSTART=%d', ndx, EVALTSTART)
                end
                tableDataRT(m,end) = models(m).dates(ndx);
                
            end
            
            filename = sprintf('%stableCR70TSTART%dEVAL%d.tex', comparisonType, TSTART, EVALTSTART);
            
            fid = fopen(filename, 'wt');
            % header
            fprintf(fid, '\\begin{normalsize}\n\\begin{center}\n');
            fprintf(fid, '\\begin{tabular}{l%sl}\n', repmat('.4', 1, Nsurvey)');
            fprintf(fid, '\\toprule\n');
            
            % list horizons
            fprintf(fid, ' & \\multicolumn{%d}{c}{Forecast horizon}\\\\\n', Nsurvey);
            fprintf(fid, '\\cmidrule{%d-%d}\n', 1+1,1+Nsurvey);
            
            fprintf(fid, 'Variable ');
            fprintf(fid, ' & \\ccol{%d} ', horizons);
            fprintf(fid, ' & \\ccol{eval. begin} ');
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
                fprintf(fid, 'The sample uses predictions made from the date given in the right-most column through 2017:Q4 (and realized forecast errors as far as available), in case of the SPF, and through 2011:Q4, in case of the Greenbook  (evaluated against realized data as far as available in both cases).\n');
            else
                fprintf(fid, 'The sample uses predictions made from the date given in the right-most column through 2017:Q4 (and realized forecast errors as far as available).\n');
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
            
        end
        
    end
    
end


%% finish / clean up
finishwrap
finishscript
dockAllFigures

