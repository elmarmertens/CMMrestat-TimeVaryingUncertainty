%% plot Fan Charts estimated from a given Hydeparkmodel


%% load toolboxes
path(pathdef)

addpath ../toolbox/emtools/
addpath ../toolbox/emtexbox/
% addpath ../toolbox/emgibbsbox/
% addpath ../toolbox/emeconometrics/

initscript
initwrap

DATALABELS = {'CPI', 'TBILL', 'UNRATE', 'PGDP', 'RGDP'};

for d = 1 : length(DATALABELS)
    
    datalabel = DATALABELS{d};
    
    %% model parameters
    
    evalTstart = 60;
    
    datadir        = '../resultfiles/'; % point to location of directory that stores mat-files generted by code in hydeparkMCMC

    modellabel = strcat(datalabel, 'sv');
    
    matfile = sprintf('slim%s01QRTdensitySVTSTART60EVALSTART60', datalabel);
    
    %% load estimates
    
    
    varlist    = {'dates', 'feSV', 'RTFEsv', ...
        'horizons', 'samStart'};
    resultData = load(fullfile(datadir, matfile), varlist{:});
    feSV       = resultData.feSV;
    RTFEsv     = resultData.RTFEsv;
    dates      = resultData.dates;
    horizons   = resultData.horizons;
    
    T = length(dates);
    
    %% load Yforecast and Yrealized
    % need to be loaded separately to match dates of result file
    ydatafile = sprintf('hydepark%sdata', upper(datalabel));
    ydata     = load(ydatafile, 'Yforecast', 'Yrealized', 'dates');
    
    ndx       = ismember(ydata.dates, dates);
    Yrealized = ydata.Yrealized(ndx,:);
    Yforecast = ydata.Yforecast(ndx,:);
        
    %% correct Yforecast sample
    % the estimation code did not truncate initial obs from Yforecast that were outside the sample
    if size(Yforecast,1) ~= T
        if size(Yforecast,1) == T + resultData.samStart - 1
            Yforecast = Yforecast(resultData.samStart:end,:);
        else
            error('sample mismatch with Yforecast?')
        end
    end
    if size(Yrealized,1) ~= T
        if size(Yrealized,1) == T + resultData.samStart - 1
            Yrealized= Yrealized(resultData.samStart:end,:);
        else
            error('sample mismatch with Yforecast?')
        end
    end
    
    
    %% padding Yrealized with blanks (for out of sample fan charts)
    Yrealized = cat(1, Yrealized, NaN(max(horizons),1));
    
    %% construct fans
    YSVupper = Yforecast + feSV;
    YSVlower = Yforecast - feSV;
    
    YRTupper = Yforecast + RTFEsv;
    YRTlower = Yforecast - RTFEsv;
    
    
    %% fan chart for a given forecast origin
    if ~isempty(wrap)
        close all
    end
    
    fontsize = 15;
    
    panelStart = find(dates == datenum(1987,1,1));
    hanni = NaN(5,1);
    
    for t = panelStart : 40 : T
        histlags = -4 : -1;
        
        if isempty(wrap)
            figure
        else
            clf reset
        end
        set(gca, 'fontsize', fontsize)
        hold on
        % forecast
        hanni(2) = plot(horizons, Yforecast(t,:), 'k:', 'linewidth', 3);
        
        % future realizations
        plot(horizons, Yrealized(t+horizons), 'kx', 'linewidth', 1, 'markersize', 8);
        hanni(3) = plot(horizons, Yrealized(t+horizons), 'ko', 'linewidth', 2, 'markersize', 8);
        
        % history
        hanni(1) = plot(histlags, Yrealized(t+histlags), 'k-', 'linewidth', 3);
        
        % SV bands
        hanni(4) = plot(horizons, YSVlower(t,:), 'r--', 'linewidth', 3);
        plot(horizons, YSVupper(t,:), 'r--', 'linewidth', 3)
        
        % RT bands
        hanni(5) = plot(horizons, YRTlower(t,:), 'b-.', 'linewidth', 3);
        plot(horizons, YRTupper(t,:), 'b-.', 'linewidth', 3)
        
        ylimits = ylim;
        switch datalabel
            case {'CPI'}
                if all(ylimits > 0)
                    ylimits(1) = 0;
                elseif all(ylimits < 0)
                    ylimits(2) = 0;
                end
        end
        ylim(ylimits)
        xlim([min(histlags) horizons(end) + 1])
        xlabel('Quarters')
        
        plotvertline(0, [], 'k:', 'linewidth', 1)
        
        set(gca, 'xtick', [histlags horizons])
        
        hl = legend(hanni, 'History', 'Forecast', 'Realized', 'ETA-SV bands', 'FE-SIMPLE bands', ...
            'location', 'best');
        wrapcf(sprintf('FAN%s%sWITHLEGEND', datestr(dates(t), 'yyyyqq'), modellabel), wrap)
        delete(hl)
        wrapcf(sprintf('FAN%s%s', datestr(dates(t), 'yyyyqq'), modellabel), wrap)
    end
    
end
%% finish
dockAllFigures
finishwrap
finishscript