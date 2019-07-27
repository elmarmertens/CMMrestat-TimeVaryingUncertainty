%% plot SV figures for Baseline model (ETA SV) vs an alternative

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


%% select choice of data: Inflation (GDPdeflator) or GDP growth


SVTSTART     = 60;
EVALTSTART   = 60;


CONSTTSTART     = 60;
CONSTEVALTSTART = EVALTSTART;

DATALABELS = {'CPI', 'TBILL', 'UNRATE', 'PGDP', 'RGDP'};

modelType    = 'QRTdensitySV'; % QRTdensity
% modelType    = 'QRTdensityVARSV';

% use AR1 as alternative
% altType      = 'QRTdensitySVAR1'; % QRTdensity
% ALTlabel     = 'SV-AR1';

% use SV1 as alternative
% altType      = 'QRTdensitySVsinglefactor'; % QRTdensity
% ALTlabel     = 'SV-1FCTR';

% use SV1 as alternative
altType      = 'QRTconst'; % QRTdensity
ALTlabel     = 'ETA-SIMPLE';

titlename = sprintf('HydeparkFiguresSVvs%s', ALTlabel);
initwrap

datadir        = '../resultfiles/';
datadirFECONST = '../resultfiles/';
datadirALT     = '../resultfiles/';


%% loop over data lables
for d = 1  : length(DATALABELS)
    
    datalabel = DATALABELS{d};
    
    matfile = sprintf('hydepark%sdata', upper(datalabel));
    basedata   = load(matfile, 'Yforecast', 'Yrealized', 'eta', 'FE', 'T', 'Nsurvey', 'dates', 'horizons', 'maxHorizon', 'Ylabel');
    
    maxHorizon = basedata.maxHorizon;
    Ylabel     = basedata.Ylabel;
    
    
    modellabel = datalabel;
    modellabel = strcat(modellabel, '01');
    
    RTmodellabel   = sprintf('%s%s', modellabel, 'QRTfeconst');
    ALTmodellabel  = sprintf('%s%s', modellabel, altType);
    modellabel     = sprintf('%s%s', modellabel, modelType);
    
    
    modellabel        = strcat(modellabel, sprintf('TSTART%dEVALSTART%d', SVTSTART, EVALTSTART));
    ALTmodellabel     = strcat(ALTmodellabel, sprintf('TSTART%dEVALSTART%d', SVTSTART, EVALTSTART));
    RTmodellabel      = strcat(RTmodellabel, sprintf('TSTART%dEVALSTART%d', CONSTTSTART, CONSTEVALTSTART));
    
    
    fprintf('Processing %s ... \n', modellabel)
    
    %% define variables to be loaded
    varlist   = {'eta', 'FE', 'T', 'Nsurvey', 'dates', 'horizons', 'feSV', 'feSVtails', ...
        'etaSV', 'etaSVFINAL', 'etaSVtails', 'etaSVFINALtails', 'RTwindow', 'RTFEsv'};
    
    %% load  estimates
    load(fullfile(datadir, strcat('slim', modellabel, '.mat')), varlist{:});
    
    checkRT = RTFEsv;
    
    %% load RTFEsv
    RTdata = load(fullfile(datadirFECONST, strcat('slim', RTmodellabel, '.mat')), 'RTFEsv', 'dates');
    if ~isequal(RTdata.dates, dates)
        error('date mismatch')
    end
    if ~isequal(isnan(RTdata.RTFEsv), isnan(checkRT))
        error('RT NaN mismatch')
    end
    
    RTFEsv = RTdata.RTFEsv;
    clear RTdata

    %% load ALTERNATIVE
    ALTdata = load(fullfile(datadirALT, strcat('slim', ALTmodellabel, '.mat')), varlist{:});
    if ~isequal(ALTdata.dates, dates)
        error('date mismatch')
    end
    if ~isequal(isnan(ALTdata.RTFEsv), isnan(checkRT))
        error('RT NaN mismatch')
    end
    
    
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
    
    
    for n = 1 : Nsurvey
        hanni = NaN(3,1);
        
        figure
        set(gca, 'FontSize', 18)
        
        
        % orient landscape
        hold on
        barcol = .7 * [1 1 1];
        hanni(1) = bar(dates, abs(fullsampleEta(:,n)), 1, 'FaceColor', barcol, 'EdgeColor', barcol);
        
        hanni(2) = plot(dates, etaSV(:,n), 'k-', 'linewidth', 5);
        plot(dates, squeeze(etaSVtails(:,n,:)), 'k-', 'linewidth', 1)
        
        hanni(3) = plot(dates, ALTdata.etaSV(:,n), 'r--', 'linewidth', 3);
        plot(dates, squeeze(ALTdata.etaSVtails(:,n,:)), 'r--', 'linewidth', 1)
        
        xtickdates(dates(CONSTEVALTSTART:end))
        maxy = max(ylim);
        ylim([0 maxy])
        
        hl = legend(hanni, '|\bf\eta|', 'SV', sprintf('%s', ALTlabel));
        wrapcf(sprintf('etaSV%d%swithlegend', n, ALTmodellabel), wrap)
        delete(hl)
        wrapcf(sprintf('etaSV%d%s', n, ALTmodellabel), wrap)
    end
    
    
    
    %% just SV with 70 and 95 CI
    hanni = NaN(4,1);
    for n = 1 : Nsurvey
        figure
        set(gca, 'FontSize', 18)
        
        hold on
        
        % plot baseline
        hanni(1) = plot(dates, FE(:,n), 'k:', 'linewidth', 2);
        
        hanni(2) = plot(dates,  1 * feSV(:,n),'r-', 'linewidth', 2);
        plot(dates, -1 * feSV(:,n),'r-', 'linewidth', 2)
        
        hanni(3) = plot(dates,  1 * RTFEsv(:,n),'b--', 'linewidth', 2);
        plot(dates, -1 * RTFEsv(:,n),'b--', 'linewidth', 2)
        
        hanni(4) = plot(dates,  1 * ALTdata.feSV(:,n),'k-.', 'linewidth', 3);
        plot(dates, -1 * ALTdata.feSV(:,n),'k-.', 'linewidth', 3)
        
        %         plot(dates,  1 * checkRT(:,n),'k-', 'linewidth', 2);
        %         plot(dates, -1 * checkRT(:,n),'k-', 'linewidth', 2)
        
        
        plothorzline(0, [], 'k:')
        %         title(sprintf('FEh=%d', horizons(n)))
        xtickdates(dates(CONSTEVALTSTART-1:end))
        
        hl = legend(hanni, 'FE', 'SV-QRT', 'FE-SIMPLE-QRT',  sprintf('%s-QRT', ALTlabel), 'location', 'best');
        wrapcf(sprintf('feSV%d%swithlegend', n, ALTmodellabel), wrap)
        delete(hl)
        wrapcf(sprintf('feSV%d%s', n, ALTmodellabel), wrap)
        
    end
    
    
    if ~isempty(wrap)
        close all
    end
end

%% finish / clean up
finishwrap
finishscript
dockAllFigures

