%% construct data set (NIPA variables) for project "hydepark"

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

%% select choice of data: Inflation (GDPdeflator) or GDP growth
DATALABELS = {'PGDP', 'RGDP'};

for d = 1 : length(DATALABELS)
    
    datalabel = DATALABELS{d};
    
    
    releaseCol = 1; % first, second, third or latest (4th)
    
    %% prepare variables
    
    % quarters are time-stamped on the first day of the quarter (FRED convention)
    % note: the following could also be automated, but might be good to set a
    % few things manually (also forces us to check things whne updating the data)
    
    dates    = genrQdates(1968,2018);
    dates    = dates(4:end-3);     % these must be the dates as stated in the SPF file (1968:Q4 through 2018:Q1)
    
    
    % major-release dates; line up with *first_second_third*csv
    mrdates = genrQdates(1965,2017);
    mrdates = mrdates(3:end); % through 2017Q4
    
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
        case 'PGDP'
            SPFimport = importdata('Mean_PGDP_Level.xlsx'); % this is the original xls file converted into csv
        case 'RGDP'
            SPFimport = importdata('Mean_RGDP_Level.xlsx'); % this is the original xls file converted into csv
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
    
    % convert SPF level forecasts into forecasts of *quarterly* inflation
    
    Yforecast = 400 * (log(SPFimport.data(:,4:8)) - log(SPFimport.data(:,3:7)));
    
    
    %% load second revision data and use it as measure of realized inflation
    % the csv files below are obtained by converting the corresponding xls
    % files from the website of FRB-PHIL
    % in particular:
    %   -- the CSV files correspond to the sheet "DATA",
    %   -- entries of "#NA" need to be changed to "-999"
    %   -- remove headers
    
    switch upper(datalabel)
        case 'PGDP'
            majorreleases = importdata('p_first_second_third.csv');
        case 'RGDP'
            majorreleases = importdata('routput_first_second_third.csv');
        otherwise
            error('datalabel %s not yet supported', upper(datalabel))
    end
    
    % check dates
    checkdates = datenum(majorreleases.textdata(1:end,1), 'YYYY:QQ');
    
    if ~isequal(checkdates, mrdates)
        error('mismatch in release dates')
    end
    
    realizeddata = majorreleases.data(:,releaseCol);
    realizeddata(realizeddata == -999) = NaN;
    
    Yrealized = NaN(T,1);
    Yrealized(ismember(dates, mrdates),1) = realizeddata(ismember(mrdates,dates));
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
    switch upper(datalabel)
        case 'PGDP'
            matfile = 'hydeparkPGDPdata';
        case 'RGDP'
            matfile = 'hydeparkRGDPdata';
        otherwise
            error('datalabel %s not yet supported', upper(datalabel))
    end
    
    save(matfile, 'Yforecast', 'Yrealized', 'eta', 'FE', 'FEobserved', ...
        'Nsurvey', 'T', 'dates', 'horizons', 'maxHorizon', 'Ylabel')
    % what
    
end