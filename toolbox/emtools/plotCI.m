function [h, hanni] = plotCI(y, tails, x, varargin)
% function plotCI(y,tails,x, yLinespec)
% plots series y against x plus confidence intervals
%   (e.g.: y is IRF and x are the lags) 
% confidence intervals are given in tails
% - each column of tails is supposed to contain fractiles of bootstraped distribution
% - tail can contain several fractiles, e.g. [2.5; 97.5; 50; 95; 16; 84];
% - the columns of tails are assumed to be ordered by fractiles (otherwise they will be sorted)
% - if tails has a single column, will interpret it as 1 std of y and construct a 95% symmetric CI
% - number of tails can odd (median), LAST column will then be treated as "median"
% CI will be shaded
% see plotCIdemo

%   Coded by  Elmar Mertens, em@elmarmertens.com

[T, N] = size(tails);
if nargin < 3
	x = 1 : T;
end

if isempty(varargin)
    yLinespec = {'k-', 'Linewidth', 2};
else
    yLinespec = varargin;
end
y = y(:);
x = x(:);
if length(y) ~= T || length(x) ~= T
	error('inconsistent input dimensions')
end

if N == 1
	tails = [y - 2 * tails; y + 2 * tails];
end

if isodd(N) % if last fractile is, say, mean/median
	single = tails(:,end);
	tails  = tails(:,1:end-1);
end
	
tails = sort(tails, 2); % note: this is just a crude swap of columns. it relies on tails being sortable

% if size(unique(i, 'rows'), 1) > 1
% 	error('tails sort not simply swapping columns')
% end

cla % CHECKME: really never needed?
p = plot(x, [y tails]);
YLIM = ylim;
delete(p);

% denan
nanny = ~any(isnan([y tails]), 2);
y     = y(nanny);
tails = tails(nanny,:);
x     = x(nanny);

hold on

hanni = area(x, [tails(:,1) diff(tails, 1, 2)], min(YLIM), 'EdgeColor', 'none');

set(hanni(1), 'facecolor', ones(1,3));

switch (length(hanni) - 1) 
   case 3
      areacolors = [.8 .4 .8];
   case 7
      areacolors = [.8 .6 .4 .2 .4 .6 .8];
   case 5
      areacolors = [.8 .6 .4 .6 .8];
      % areacolors = [.75 .5 0 .5 .75];
   case 1
      areacolors = .8;
   otherwise
      error('unprepared for this number of tails ...')
end

for n = 2 : length(hanni)
   set(hanni(n), 'facecolor', repmat(areacolors(n - 1),1,3));
end

if isodd(N)
	plot(x,single, 'w-', 'MarkerSize', 3)
end

xhanni = plot(x, y, yLinespec{:});
set(gca,'Layer','top')

if nargout > 0
    h = xhanni;
end
