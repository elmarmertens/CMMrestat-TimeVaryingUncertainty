%% plot ETAs

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

% in order to store graphics as files, define wrap as a structure with field "dir"
% graphics will then be stored in directory wrap.dir
% wrap.dir = pwd; % in order to store graphs in current directory

%% loop over data

wrap = [];
initwrap % comment this line out to avoid storing figures
diary2wrap(wrap, 'screen.log')

DATALABELS = {'CPI', 'TBILL', 'UNRATE', 'PGDP', 'RGDP'};

NSURVEY = 5;

for d = 1  : length(DATALABELS)
    
    close all
    datalabel = DATALABELS{d};
    
    
    fprintf('Processing %s ... \n', datalabel)
    
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
    
    hrulefill
    fprintf('%s: \n', datalabel)
    
    fprintf('Data from %s to %s\n', datestr(dates(1), 'YYYY:QQ'), datestr(dates(end), 'YYYY:QQ'))
    
    estimationStart = find(~(etaNanny(:,1)),1);
    fprintf('Estimation sample from %s to %s\n', datestr(dates(estimationStart), 'YYYY:QQ'), datestr(dates(end), 'YYYY:QQ'))
    for RTwindow = [40 60 80]
        fprintf('RTwindow %d: Evaluation sample from %s to %s\n', RTwindow, datestr(dates(estimationStart+RTwindow) , 'YYYY:QQ'), datestr(dates(end), 'YYYY:QQ'))
    end
    
    fprintf('The number of missing obs (for each horizon) is:\n')
    fprintf('\t %d ', 0:size(etaNanny,2)-1)
    fprintf('\n');
    fprintf('\t %d ', sum(etaNanny))
    fprintf('\n');
    
    hrulefill
    
    
    
    %% plot individual figures
    hanni = NaN(2,1);
    for n = 1 : NSURVEY
        
        % FE and ETA
        figure
        set(gca, 'FontSize', 18)
        hold on
        hanni(1) = plot(dates, eta(:,n), 'r-', 'linewidth', 2);
        hanni(2) = plot(dates, FEobserved(:,n), 'k-.', 'linewidth', 2);
        
        nbershades(dates)
        hl = legend(hanni, '\eta', 'FE', 'location', 'best');
        wrapcf(sprintf('dataFETA%sh%dwithlegend', datalabel,n), wrap)
        delete(hl)
        wrapcf(sprintf('dataFETA%sh%d', datalabel,n), wrap)
        
        % ETA
        figure
        set(gca, 'FontSize', 18)
        hold on
        hanni(1) = plot(dates, eta(:,n), 'r-', 'linewidth', 2);
        nbershades(dates)
        wrapcf(sprintf('dataETA%sh%d', datalabel,n), wrap)
        
        % FE
        figure
        set(gca, 'FontSize', 18)
        hold on
        hanni(1) = plot(dates, FEobserved(:,n), 'k-.', 'linewidth', 2);
        nbershades(dates)
        wrapcf(sprintf('dataFE%sh%d', datalabel,n), wrap)
        
    end
    
    
    %% multi-panel plot
    tickvec = datenum(1965:10:2025,1,1);
    
    ndx = 2:5; %[1 2 3 5];
    
    
    hanni = NaN(length(ndx));
    figure
    for n = 1 : length(ndx)
        
        subplot(2,2,n)
        set(gca, 'FontSize', 14)
        hold on
        plot(dates, eta(:,ndx(n)), 'r-', 'linewidth', 2);
        hanni(n) = plot(dates, FEobserved(:,ndx(n)), 'k-.', 'linewidth', 1);
        nbershades(dates,[],[],tickvec)
        title(sprintf('h = %d', ndx(n) - 1))
    end
    orient('landscape')
    wrapcf(sprintf('dataFETAall%s', datalabel), wrap)
    
    figure
    for n = 1 : length(ndx)
        
        subplot(2,2,n)
        set(gca, 'FontSize', 14)
        hold on
        plot(dates, eta(:,ndx(n)), 'r-', 'linewidth', 2);
        nbershades(dates,[],[],tickvec)
        title(sprintf('h = %d', ndx(n) - 1))
    end
    orient('landscape')
    wrapcf(sprintf('dataETAall%s', datalabel), wrap)
    
    
    ndx = 2:5;
    figure
    for n = 1 : length(ndx)
        
        subplot(2,2,n)
        set(gca, 'FontSize', 14)
        hold on
        plot(dates, FEobserved(:,ndx(n)), 'k-.', 'linewidth', 2);
        nbershades(dates,[],[],tickvec)
        title(sprintf('h = %d', ndx(n) - 1))
    end
    orient('landscape')
    wrapcf(sprintf('dataFEall%s', datalabel), wrap)
    
    
    
end

%% finish / clean up
finishwrap
finishscript
dockAllFigures

