function dwellTime = timecount(projWidths, strokeTime, scaleDivs, margins)
%TIMECOUNT dwell time counts
%
% varargin:
%   projWidths   --  projected widths
%   strokeTime  --  elapsed time per stroke
%   scaleDivs    --  scale divisions on the grating substrate
%   margins      --  head and tail margins
%
% varargout:
%   dwellTime  --  dwell time
%

% copyright (c) wulx, gurdy.woo@gmail.com
% last modified by wulx, 2013/11/18


if any(~isfinite(projWidths) | projWidths<0) 
    error('projected width should be finite non-negtive.')
end

if ~isfinite(strokeTime) || strokeTime<0
    error('elapsed time should be finite non-negtive.')
end

% if any(~isfinite(scaleDivs) | scaleDivs<0)
%     error('substrate scale divisions should be finite non-negtive.')
% end

if nargin < 4, margins = [0, 0]; end % default no margins

nScan = numel(projWidths);

firstPos = scaleDivs(1) - margins(1);
lastPos = scaleDivs(end) + margins(2);

% stroke = lastPos - firstPos;
% vLeaf = stroke / strokeTime; % constant speed

leafCenterPos = linspace(firstPos, lastPos, nScan);
counts = zeros(size(scaleDivs));
for i = 1:nScan
    pos = leafCenterPos(i);
    halfWidth = 0.5*projWidths(i);
    hits = scaleDivs>(pos-halfWidth) & scaleDivs<=(pos+halfWidth);
    counts = counts + hits;
end

timePerScan = strokeTime / nScan;
dwellTime = timePerScan * counts;

