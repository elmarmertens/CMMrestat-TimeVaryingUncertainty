%% load data set and perform estimates using "ETA" notation
% 01 version: assumes that elements of eta are serially uncorrelated;
%             estimates contemporaneous correlation (as in "Extension 1")

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

% in order to store graphics as files, define wrap as a structure with field "dir"
% graphics will then be stored in directory wrap.dir
% wrap.dir = pwd; % in order to store graphs in current directory

%% select choice of data: Inflation (GDPdeflator) or GDP growth

Nstreams    = max(1,getparpoolsize);
rndStreams  = initRandStreams(Nstreams, [], 0);

DATALABELS = {'CPI', 'TBILL', 'UNRATE', 'PGDP', 'RGDP'};

tic

for d =   1  : length(DATALABELS)
    
    
    
    
    
    close all
    datalabel = DATALABELS{d};
    
    
    % construct a model-label indicating dataset and important parameters, to
    % append to picture names
    modellabel = datalabel;
    %#ok<*UNRCH>
    modellabel = strcat(modellabel, '01'); % legacy naming convention
    
    modellabel = strcat(modellabel,'QRTdensityFEVARSV');
        
    
    wrap = [];
    titlename = modellabel;

    
    fprintf('Processing %s ... \n', modellabel)
    
    %% load data
    matfile = sprintf('hydepark%sdata', upper(datalabel));
    load(matfile, 'Yforecast', 'Yrealized', 'FE', 'FEobserved', 'T', 'Nsurvey', 'dates', 'horizons', 'maxHorizon', 'Ylabel')
    FEnanny = isnan(FEobserved);
    

    
    %% compute IC for different lag length choices
    
    % drop initial obs that have NaN
    ndx = any(isnan(FEobserved),2);
    samStart = find(ndx,1,'last') + 1;
    maxlag = 8;
    [AIC,SIC, HQIC] = VARlagIC(FEobserved(samStart:end,:), 1:maxlag);
    
    [~, minAIC]  = min(AIC);
    [~, minSIC]  = min(SIC);
    [~, minHQIC] = min(HQIC);
    
    hrulefill
    fprintf('%s\n', datalabel)
    fprintf('Using continuous sample since %s\n', datestr(dates(samStart)))
    fprintf('%10s \t', 'lags', 'AIC', 'BIC', 'HQIC')
    fprintf('\n')
    for k = 1 : maxlag
        fprintf('%10d', k)
        fprintf('\t%10.4f', AIC(k));
        if k == minAIC
            fprintf('*')
        end
        fprintf('\t%10.4f', SIC(k));
        if k == minSIC
            fprintf('*')
        end
        fprintf('\t%10.4f', HQIC(k));
        if k == minHQIC
            fprintf('*')
        end
        fprintf('\n')
    end
    
    fprintf('\n')
    hrulefill
end

toc

%% finish / clean up
finishscript
dockAllFigures
