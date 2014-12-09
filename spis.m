function [xx, ys] = spis(ogee, tol, m)
%SPIS SPline Interpolation (or extrapolation) and Smoothing
%

narginchk(1, 3);

if nargin < 3, m = 3; end
if nargin < 2, tol = 1.e-2; end

n = numel(ogee);

x = (1:n)'-0.5; % center points
xx = [0; x; n]; % add boundaries

% exclude NANs
notNan = ~isnan(ogee);
x = x(notNan);
ogee = ogee(notNan);

% extrapolation to envelop NAN points and two boundaries
yy = spline(x, ogee, xx);

% spline smoothing
[~, ys] = spaps(xx, yy, tol, m);
ys = ys'; % row to column

