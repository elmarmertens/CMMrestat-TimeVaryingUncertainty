%% plot SV figures for Baseline model (ETA SV)

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

wrap = [];
initwrap


%% select choice of data: Inflation (GDPdeflator) or GDP growth


SVTSTART     = 60;
EVALTSTART   = 60;


CONSTTSTART     = 60;
CONSTEVALTSTART = EVALTSTART;

DATALABELS = {'CPI', 'TBILL', 'UNRATE', 'PGDP', 'RGDP'};

modelType    = 'QRTdensitySV';
% modelType    = 'QRTdensityVARSV';

datadir        = '../resultfiles/'; % point to location of directory that stores mat-files generted by code in hydeparkMCMC

fontsize = 22;

%% loop over data lables
for d = 1  : length(DATALABELS)
    
    close all
    datalabel = DATALABELS{d};
    
    matfile = sprintf('hydepark%sdata', upper(datalabel));
    basedata   = load(matfile, 'Yforecast', 'Yrealized', 'eta', 'FE', 'T', 'Nsurvey', 'dates', 'horizons', 'maxHorizon', 'Ylabel');
    
    maxHorizon = basedata.maxHorizon;
    Ylabel     = basedata.Ylabel;
    
    
    modellabel = datalabel;
    modellabel = strcat(modellabel, '01');
    
    RTmodellabel   = sprintf('%s%s', modellabel, 'QRTfeconst');
    modellabel     = sprintf('%s%s', modellabel, modelType);
    
    
    modellabel        = strcat(modellabel, sprintf('TSTART%dEVALSTART%d', SVTSTART, EVALTSTART));
    RTmodellabel      = strcat(RTmodellabel, sprintf('TSTART%dEVALSTART%d', CONSTTSTART, CONSTEVALTSTART));
    
    
    fprintf('Processing %s ... \n', modellabel)
    
    %% load  estimates
    % define variables to be loaded from prior estimation run
    if strcmpi(modelType, 'QRTdensityVARSV')
        varlist   = {'eta', 'FE', 'T', 'Nsurvey', 'dates', 'horizons', 'feRMSE', ...
            'RTwindow', 'RTFEsv'};
        load(fullfile(datadir, strcat('slim', modellabel, '.mat')), varlist{:});
        feSV   = feRMSE;
    else
        varlist   = {'eta', 'FE', 'T', 'Nsurvey', 'dates', 'horizons', 'feSV', 'feSVtails', ...
            'etaSV', 'etaSVFINAL', 'etaSVtails', 'etaSVFINALtails', 'RTwindow', 'RTFEsv'};
        load(fullfile(datadir, strcat('slim', modellabel, '.mat')), varlist{:});
    end
    
    checkRT = RTFEsv;
    
    %% load RTFEsv
    RTdata = load(fullfile(datadir, strcat('slim', RTmodellabel, '.mat')), 'RTFEsv', 'dates');
    if ~isequal(RTdata.dates, dates)
        error('date mismatch')
    end
    if ~isequal(isnan(RTdata.RTFEsv), isnan(checkRT))
        error('RT NaN mismatch')
    end
    
    RTFEsv = RTdata.RTFEsv;
    clear RTdata
    
    
    %% ETA data
    % the estimation code has set eta to NaN prior to the evaluation window
    % for the plots, we would like to show the full eta matrix however
    
    base2model    = ismember(basedata.dates, dates);
    fullsampleEta = basedata.eta(base2model,:);
    
    checkdiff(eta, fullsampleEta); % automatically ignoring the NaN
    
    %% report sample
  
    hrulefill
    fprintf('%s: \n', datalabel)
    
    fprintf('Estimation sample from %s to %s\n', datestr(dates(1), 'YYYY:QQ'), datestr(dates(end), 'YYYY:QQ'))
    fprintf('RTwindow %d: Evaluation sample from %s to %s\n', RTwindow, datestr(dates(RTwindow+1) , 'YYYY:QQ'), datestr(dates(end), 'YYYY:QQ'))
    
    
    %% abs ETA vs SV
    
    if ~strcmpi(modelType, 'QRTdensityVARSV')
        for n = 1 : Nsurvey
            hanni = NaN(3,1);
            
            figure
            set(gca, 'FontSize', fontsize)
            
            % orient landscape
            hold on
            barcol = .7 * [1 1 1];
            hanni(1) = bar(dates, abs(fullsampleEta(:,n)), 1, 'FaceColor', barcol, 'EdgeColor', barcol);
            
            % hanni(2) = plot(dates, etaSV(:,n), 'k-', 'linewidth', 3);
            % plot(dates, squeeze(etaSVtails(:,n,:)), 'k-', 'linewidth', 1)
            
            hanni(3) = plot(dates, etaSVFINAL(:,n,end), 'r-.', 'linewidth', 5);
            hanni(2) = plot(dates, etaSV(:,n), 'k-', 'linewidth', 3);
            %         plot(dates, squeeze(etaSVFINALtails(:,n,:,end)), 'k--', 'linewidth', 1)
            
            % xtickdates(dates(CONSTEVALTSTART:end))
            xtickdates(dates, 'keepticks')
            maxy = max(ylim);
            ylim([0 maxy])
            
            hl = legend(hanni, '|\bf\eta|', 'SV: QRT', 'SV: Final', 'location', 'best');
            wrapcf(sprintf('etaSV%d%swithlegend', n, modellabel), wrap)
            delete(hl)
            wrapcf(sprintf('etaSV%d%s', n, modellabel), wrap)
        end
        
    end
    
    %% just SV with 70 and 95 CI
    hanni = NaN(3,1);
    for n = 1 : Nsurvey
        figure
        set(gca, 'FontSize', fontsize)
        
        hold on
        
        % plot baseline
        hanni(1) = plot(dates, FE(:,n), 'k:', 'linewidth', 2);
        
        hanni(2) = plot(dates,  1 * feSV(:,n),'r-', 'linewidth', 2);
        plot(dates, -1 * feSV(:,n),'r-', 'linewidth', 2)
        
        hanni(3) = plot(dates,  1 * RTFEsv(:,n),'b--', 'linewidth', 2);
        plot(dates, -1 * RTFEsv(:,n),'b--', 'linewidth', 2)
        
        %         plot(dates,  1 * checkRT(:,n),'k-', 'linewidth', 2);
        %         plot(dates, -1 * checkRT(:,n),'k-', 'linewidth', 2)
        
        
        plothorzline(0, [], 'k:')
        %         title(sprintf('FEh=%d', horizons(n)))
        xtickdates(dates(CONSTEVALTSTART-1:end), 'keepticks')
        
        hl = legend(hanni, 'FE', 'SV-QRT', 'FE-SIMPLE-QRT', 'location', 'best');
        wrapcf(sprintf('feSV%d%swithlegend', n, modellabel), wrap)
        delete(hl)
        wrapcf(sprintf('feSV%d%s', n, modellabel), wrap)
        
    end
    
    
    if ~isempty(wrap)
        close all
    end
end

%% finish / clean up
finishwrap
finishscript
dockAllFigures

