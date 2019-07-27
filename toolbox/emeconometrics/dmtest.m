function [mu, tstat, pvalue] = dmtest(loss1, loss2, nlag)
% DMTEST ... 
%  
%   ... 

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 28-Jul-2017 20:32:48 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 9.2.0.556344 (R2017a) 
% FILENAME  : dmtest.m 


nanny1 = isnan(loss1);
nanny2 = isnan(loss2);

if ~isequal(nanny1, nanny2)
    error('mismatch in obs')
end

delta = loss1(~nanny1) - loss2(~nanny2);
Nobs  = length(delta);

if nargin < 3 || isempty(nlag)
    nlag = floor( 4 * (Nobs / 100)^(2/9) );
end


reggae = nwest(delta, ones(Nobs,1), nlag);

mu      = reggae.beta;
tstat   = reggae.tstat;
pvalue  = tdis_prb(tstat,Nobs-1);

% fprintf('Doing DM test with %d NW lags\n', nlag)
