%% compare results: ETA-SV vs ETA-SIMPLE based on CRPS

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


for COMPARISONS = {'SVvsCONST'}
    
    comparisonType = COMPARISONS{:};
    jointFlag = false; % default, unless overriden by one of the switch clauses below
    
    % default directories
    datadir        = '../resultfiles/';
    datadirFECONST = '../resultfiles/';
    
    varlist   = {'eta', 'FE', 'T', 'Nsurvey', ...
        'dates', 'horizons',  ...
        'crps'};
    
    varlistCONST = varlist;
    
    switch comparisonType
        case 'SVvsCONST'
            DATALABELS = {'RGDP', 'UNRATE', 'PGDP', 'CPI', 'TBILL'};
            DATAPRETTYLABELS = {'RGDP', 'UNRATE', 'PGDP', 'CPI', 'TBILL'};
            baselineType  = 'densitySV';
            baselineLabel = 'SV';
        otherwise
            error('comparisonType %s not supported', comparisonType)
    end
    
    
    
    dummy    = cat(1, varlist, cell(1,length(varlist)));
    models   = repmat(struct(dummy{:}), length(DATALABELS), 1);
    dummy    = cat(1, varlistCONST, cell(1,length(varlistCONST)));
    modelsCONST = repmat(struct(dummy{:}), length(DATALABELS), 1);
    
    
    TSTART = 60;
    EVALTSTART = TSTART;
    
    for d = 1 : length(DATALABELS)
        
        datalabel = DATALABELS{d};
        
        modellabel = datalabel;
        %#ok<*UNRCH>
        
        modellabel = strcat(modellabel, '01');
        
        
        
        modellabelCONST   = sprintf('%sQRTconst', modellabel);
        modellabelSV      = sprintf('%sQRT%s', modellabel, baselineType);
        
        
        modellabelCONST = strcat(modellabelCONST, sprintf('TSTART%dEVALSTART%d', TSTART, EVALTSTART));
        modellabelSV = strcat(modellabelSV, sprintf('TSTART%dEVALSTART%d', TSTART, EVALTSTART));
        
        
        %% load results
        matfilename = strcat('slim', modellabelSV, '.mat');
        models(d)   = load(fullfile(datadir, matfilename), varlist{:});
        
        matfilename  = strcat('slim', modellabelCONST, '.mat');
        modelsCONST(d)  = load(fullfile(datadirFECONST, matfilename), varlistCONST{:});
        if ~isequal(modelsCONST(d).dates, models(d).dates)
            error('date mismatch')
        end
        
        
        %% consistency check between CONST and SV models
        if ~isequal(isnan(modelsCONST(d).crps), isnan(models(d).crps))
            error('CRPS NaN mismatch')
        end
        
    end % datalabel
    
    %% collect a few parameters
    % (Assuming model(1) is representative)
    Nsurvey     = models(1).Nsurvey;
    horizons    = models(1).horizons;
    
    %% CRPS comparison table
    tableData   = NaN(length(DATALABELS), Nsurvey + 1);
    tableData2  = NaN(length(DATALABELS), Nsurvey + 1);
    [crpsDMmu, crpsDMtstat, crpsDMpval] = deal(NaN(length(DATALABELS), Nsurvey));
    for m = 1 : length(DATALABELS)
        
        
        for h = 1 : Nsurvey
            [crpsDMmu(m,h), crpsDMtstat(m,h), crpsDMpval(m,h)] = dmtest(modelsCONST(m).crps(:,h),models(m).crps(:,h),h+1); % note that h+! is the same as ceil(1.5 + horizons(h))
        end
        
        
        crpsCONST = nanmean(modelsCONST(m).crps);
        crpsSV    = nanmean(models(m).crps);
        tableData(m,1:Nsurvey)  = (crpsCONST - crpsSV) ./ crpsCONST * 100;
        tableData2(m,1:Nsurvey) = crpsCONST;
        
        
        checkdiff(crpsCONST - crpsSV, crpsDMmu(m,:));
        
        ndx = find(any(~isnan(models(m).FE),2), 1);
        if ~isequal(ndx, EVALTSTART)
            warning('Evaluation window seems to have started at t=%d, not EVALTSTART=%d', ndx, EVALTSTART)
        end
        tableData(m,end) = models(m).dates(ndx);
    end
    
    filename = sprintf('%stableCRPSTSTART%dEVAL%d.tex', comparisonType, TSTART, EVALTSTART);
    
    fid = fopen(filename, 'wt');
    % header
    fprintf(fid, '\\begin{normalsize}\n\\begin{center}\n');
    fprintf(fid, '\\begin{tabular}{l%sl}\n', repmat('.5', 1, Nsurvey)');
    fprintf(fid, '\\toprule\n');
    
    % list horizons
    fprintf(fid, '  & \\multicolumn{%d}{c}{Forecast horizon}\\\\\n', Nsurvey);
    fprintf(fid, '\\cmidrule{%d-%d}\n', 1+1,1+Nsurvey);
    
    fprintf(fid, 'Variable  ');
    fprintf(fid, ' & \\ccol{%d} ', horizons);
    fprintf(fid, ' & \\ccol{eval. begin} ');
    fprintf(fid, '\\\\\n'); % new line
    fprintf(fid, '\\midrule\n');
    
    for m = 1 : size(tableData,1)
        
        fprintf(fid, '\\textbf{%s}  \\\\\n ', DATAPRETTYLABELS{m});
        fprintf(fid, '\\quad (%s rel.) ', baselineLabel);
        
        for n = 1 : Nsurvey
            fprintf(fid, '& %5.2f\\%%%s ', tableData(m,n), Pstar(crpsDMpval(m,n)));
        end
        fprintf(fid, '& %s ', datestr(tableData(m,end), 'YYYY:QQ'));
        fprintf(fid, '\\\\\n');
        
        fprintf(fid, '\\quad (ETA-SIMPLE) ');
        for n = 1 : Nsurvey
            fprintf(fid, '& %5.2f ', tableData2(m,n));
        end
        fprintf(fid, '\\\\\n');
    end
    
    
    % footer
    % evalStart = modelsCONST(m).dates(TSTART);
    evalStop  = models(m).dates(end);
    fprintf(fid, '\\bottomrule\n');
    fprintf(fid, '\\end{tabular}\n');
    fprintf(fid, '\\end{center}\n');
    
    % baseline labels
    H0model   = 'ETA-SIMPLE';
    % adjust as needed
    switch comparisonType
        case {'SVvsCONST'}
            thisModel = 'ETA-SV';
        otherwise
            error('comparisonType %s not defined', comparisonType)
    end
    
    fprintf(fid, 'Note: The table reports CRPS results for out-of-sample density forecasts.\n');
    if any(contains(DATALABELS, 'Greenbook'))
        fprintf(fid, 'The sample uses predictions made from the date given in the right-most column through 2017:Q4 (and realized forecast errors as far as available), in case of the SPF, and through 2011:Q4, in case of the Greenbook  (evaluated against realized data as far as available in both cases).\n');
    else
        fprintf(fid, 'The sample uses predictions made from the date given in the right-most column through 2017:Q4 (and realized forecast errors as far as available).\n');
    end
    fprintf(fid, 'For each variable, the top row reports the relative CRPS calculated as the percentage decrease of the CRPS when using %s rather than %s; positive numbers indicate improvement of %s over the %s case.\n', ...
        thisModel, H0model,thisModel, H0model);
    fprintf(fid, 'The bottom row reports the CRPS for the %s case.\n', ...
        H0model);
    fprintf(fid, 'Both %s and %s have been estimated over rolling windows with %d~quarterly observations.\n', ...
        thisModel, H0model, TSTART);
    
    
    fprintf(fid, 'Statistical significance of the differences in average CRPS --- assessed with a Diebold and Mariano (1995) test --- is indicated by *, **, or ***, corresponding to 10, 5, and 1 percent significance, respectively.\n');
    fprintf(fid, '\\end{normalsize}\n');
    
    fclose(fid);
    type(filename)
    display(filename) %#ok<DSPS>
    
    
    if ~isempty(wrap)
        movefile(filename, wrap.dir)
        latexwrapper(wrap, 'add', 'table', filename);
    end
    
end



%% finish / clean up
finishwrap
finishscript
dockAllFigures

