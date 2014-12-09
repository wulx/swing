function [locs, ptype] = findpoints(ys, opts, popts, vopts)
%FINDPOINTS find target points, including peaks and valleys, and special
%inflection points
%
% p.s. to divide the spline into S shaped curves, we should find all local 
% maxima(peaks), local minima(valleys) and special inflection points. The 
% inflection points refer specifically to those inflection points at ogee 
% curves.
%


if nargin < 2
    opts = {'MINPEAKHEIGHT', -1, 'MINPEAKDISTANCE', 15, 'NPEAKS', 5};
end
if nargin < 3
    popts = {'MINPEAKDISTANCE', 15, 'NPEAKS', 5};
end
if nargin < 4
    vopts = {'MINPEAKDISTANCE', 15, 'NPEAKS', 5};
end

% find all target points
y1 = -abs(diff(ys, 1));
[~, locs] = findpeaks(y1, opts{:});

% % exclude the points out of mean-std
% locs(pks < mean(pks)-std(pks,1)) = [];

% find peaks
[~, plocs] = findpeaks(ys, popts{:});

% plot(xx(plocs), ys(plocs), 'ro')

% find valleys
[~, vlocs] = findpeaks(-ys, vopts{:});

% plot(xx(vlocs), ys(vlocs), 'bo')

% others are the infletion points at ogees
% point type:
%   1   --  peaks
%   0   --  inflection points at ogees
%   -1  --  valleys
ptype = zeros(size(locs));
d = 1.5;
for i = 1:numel(locs)
    if any((plocs < locs(i)+d) & (plocs > locs(i)-d))
        ptype(i) = 1;
    elseif any((vlocs < locs(i)+d) & (vlocs > locs(i)-d))
        ptype(i) = -1;
    end
end
