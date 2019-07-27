%% construct data set for project "hydepark"
% all variables, except for GDP and PGDP (NIPA data matched against RTDSM)

%% load toolboxes
path(pathdef)

addpath ../toolbox/emtools/
% addpath ../toolbox/emtexbox/
% addpath ../toolbox/emgibbsbox/
% addpath ../toolbox/emeconometrics/

%% clear workspace
clear variables
clear global
close all
fclose all;
clc

%% select choice of data
DATALABELS = {'UNRATE', 'CPI', 'TBILL'};

for d = 1 : length(DATALABELS)
    
    datalabel = DATALABELS{d};
    
    %% prepare variables
    
    % quarters are time-stamped on the first day of the quarter (FRED convention)
    % note: the following could also be automated, but might be good to set a
    % few things manually (also forces us to check things whne updating the data)
    
    dates    = genrQdates(1968,2018);
    dates    = dates(4:end-3);     % these must be the dates as stated in the SPF file (1968:Q4 through 2018:Q1)
    
    
    T        = length(dates);
    
    maxHorizon = 4;
    horizons   = 0 : maxHorizon;     % SPF forecast horizons
    Nsurvey  = length(horizons);
    
    Ylabel = cell(Nsurvey,1);
    for n = 1 : Nsurvey
        Ylabel{n} = sprintf('h%d', horizons(n));
    end
    
    %% load SPF responses for the level of PGDP and convert into inflation rates
    % the results of this cell are to be stored in Ydata(:,2:6)
    
    switch upper(datalabel)
        case 'UNRATE'
            SPFimport = importdata('Mean_UNEMP_Level.xlsx'); % this is the original xls file converted into csv
        case 'TBILL'
            SPFimport = importdata('Mean_TBILL_Level.xlsx'); % this is the original xls file converted into csv
        case {'CPI'}
            SPFimport = importdata('Mean_CPI_Level.xlsx'); % this is the original xls file converted into csv
        otherwise
            error('datalabel %s not yet supported', upper(datalabel))
    end
    
    % col 1: year
    % col 2: quarter
    % col 3: "forecasts" for previous quarter
    % col 4: nowcast
    % col 5: forecast for next quarter
    % col 6: 2-quarter forecast
    % col 7: 3-quarter forecast
    % col 8: 4-quarter forecast
    % col 9: Annual forecast for current year
    % col 10: Annual forecast for next year
    
    
    % convert -999 into NaN
    SPFimport.data(SPFimport.data == -999) = NaN;
    
    % check dates
    Y = SPFimport.data(:,1);
    if ~isequal(Y, year(dates))
        error('inconsistent date vector in the SPF file (years)');
    end
    Q = SPFimport.data(:,2);
    if ~isequal(Q, quarter(dates))
        error('inconsistent date vector in the SPF file (quarters)');
    end
    
    % pull out SPF forecasts for quarter 0 through 4
    Yforecast = SPFimport.data(:,4:8);
    
    
    %% load realized data from current FRED vintage
    
    switch upper(datalabel)
        case 'UNRATE'
            FRED = fredreadcsv('UNRATE',[], 'm', 'q', 'avg');
            Yrealized = fredexpand(FRED, dates);
        case 'TBILL'
            FRED = fredreadcsv('TB3MS',[], 'm', 'q', 'avg');
            Yrealized = fredexpand(FRED, dates);
        case {'CPI'}
            FRED = fredreadcsv('CPIAUCSL',[], 'm', 'q', 'avg');
            pricelevel = fredexpand(FRED, cat(1, datenum(1968,7,1), dates));
            Yrealized  = ((pricelevel(2:end) ./ pricelevel(1:end-1)).^4 - 1) * 100;
            if dates(1) ~= datenum(1968,10,1)
                error('date mismatch')
            end
            if length(Yrealized) ~= T
                error('date mismatch')
            end
        otherwise
            error('datalabel %s not yet supported', upper(datalabel))
    end
    
    
    % RAW data is now complete
    
    
    %% construct Forecast errors --- timing: future realizations matched to current forecast!
    % timing: congruent with realized values
    FE = NaN(T, Nsurvey);
    
    for t = 1 : T
        for h = 0 : min(T-t,maxHorizon)
            hcol = h + 1;
            FE(t,hcol) = Yrealized(t+h) - Yforecast(t,hcol);
        end
    end
    
    FEobserved = NaN(T, Nsurvey);
    
    for t = 2 : T
        for h = 0 : maxHorizon
            hcol = h + 1;
            tForecast = t - 1 - h;
            if tForecast > 0
                FEobserved(t,hcol) = Yrealized(t-1) - Yforecast(tForecast,hcol);
            end
        end
    end
    %% construct eta vector of expectational updates
    eta = NaN(T,Nsurvey);
    % *lagged* nowcast error
    for t = 2 : T
        eta(t,1) = Yrealized(t-1) - Yforecast(t-1,1);
    end
    
    % survey updates
    for j = 0 : maxHorizon-1
        for t = 2 : T
            eta(t,2+j) = Yforecast(t,1+j) - Yforecast(t-1,2+j);
        end
    end
    
    %% verify: FE replicated from eta
    FE2 = NaN(size(FE));
    
    % % h = 0
    % for t = 1 : T-1
    %     FE2(t,1) = eta(t+1,1);
    % end
    % % h = 1
    % for t = 1 : T-2
    %     FE2(t,2) = eta(t+2,1) + eta(t+1,2);
    % end
    % % h = 2
    % for t = 1 : T-3
    %     FE2(t,3) = eta(t+3,1) + eta(t+2,2) + eta(t+1,3);
    % end
    
    for h = 0 : maxHorizon
        for t = 1 : T-h-1
            FE2(t,h+1) = 0;
            for j = 0 : h
                FE2(t,h+1) = FE2(t,h+1) + eta(t+1+h-j,j+1);
            end
        end
    end
    checkdiff(FE,FE2);
    
    
    
    %% store data
    matfile = sprintf('hydepark%sdata', upper(datalabel));
    
    
    save(matfile, 'Yforecast', 'Yrealized', 'eta', 'FE', 'FEobserved', ...
        'Nsurvey', 'T', 'dates', 'horizons', 'maxHorizon', 'Ylabel')
    
end