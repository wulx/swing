function [r, dwellTime] = dwellopt(wa, wd, wf, w1, w2, crt, strokeTime, timeStep, stepAngleDeg, leafWidth)
%DWELLOPT DWELL time OPTimization
% varargin:
%   @params to be optimized
%   wa   --  weight factors of number of acceleration steps, (0, 1)
%   wd   --  weight factors of number of deceleration steps, (0, 1)
%   wf   --  weight factors of initial frequency, (0, 1)
%   w1   --  to enlarge amplitude of peaks or valleys
%   w2   --  to braoden step differences
%   @params to set
%   ogee          --  initial dwell time data
%   strokeTime    --  time per a stroke [s]
%   timeStep      --  time step
%   stepAngleDeg  --  step angle in degree
%   leafWidth     --  width of the dynamic leaf
% varargout:
%   r           --  RMSD of dwell time
%   dwellTime   --  calculated dwell time

% by wulx, 2014/7/14

% global r dwellTime

[stepNums, timeSeqs] = stepsamp(wa, wd, wf, w1, w2, crt, strokeTime);

% timeStep = 0.001; % 1 ms
nScan = ceil(strokeTime / timeStep);

stepList = zeros(nScan, 1);
timeline = linspace(0, strokeTime, nScan)';

nStep = numel(timeSeqs);
for i = 1:nStep
    stepList(timeline>=sum(timeSeqs(1:i-1)) & timeline<=sum(timeSeqs(1:i))) = stepNums(i);
end

% projected withs of the dynamic leaf
projWidths = step2width(stepList, stepAngleDeg, leafWidth);

% caculated dwell times
dwellTime = timecount(projWidths, strokeTime, crt.scaleDivs);

% root-mean-square deviation of dwell time
dwellTime = dwellTime(2:end-1);
% dwellTime(isnan(crt.ogee)) = nan;

r = rmsd(crt.ogee, dwellTime);

