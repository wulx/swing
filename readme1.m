%% Trajectory optimization using pattern search
%
% Back to the <matlab:webExp('html/readme.html') main> page.

%%
% *initialization*
close all; clear all; clc


%% Dwell time optimization
% Conduct dwell time optimization using Pattern Search method
%
%   [r, dwellTime] = dwellopt(wa, wd, wf, w1, w2, crt, strokeTime, timeStep, stepAngleDeg, leafWidth)
%
%  varargin:
%    @params to be optimized
%      wa   --  weight factors of number of acceleration steps, (0, 1)
%      wd   --  weight factors of number of deceleration steps, (0, 1)
%      wf   --  weight factors of initial frequency, (0, 1)
%      w1   --  to enlarge amplitude of peaks or valleys
%      w2   --  to braoden step differences
%    @params to set
%      ogee          --  initial dwell time data
%      strokeTime    --  time per a stroke [s]
%      timeStep      --  time step
%      stepAngleDeg  --  step angle in degree
%      leafWidth     --  width of the dynamic leaf
%  varargout:
%    r           --  RMSD of dwell time
%    dwellTime   --  calculated dwell time
%
% Source code:
%
% * <matlab:openInEditor('dwelloptRun.m') *dwelloptRun.m*>
% * + <matlab:openInEditor('spdiv.m') spdiv.m>
% * ++ <matlab:openInEditor('spis.m') spis.m>
% * ++ <matlab:openInEditor('findpoints.m') findpoints.m>
% * + <matlab:openInEditor('spdivplot.m') spdivplot.m>
% * <matlab:openInEditor('dwellopt.m') *dwellopt.m*>
% * + <matlab:openInEditor('stepsamp.m') stepsamp.m>
% * ++ <matlab:openInEditor('timeperstepd.m') timeperstepd.m>
% * + <matlab:openInEditor('step2width.m') step2width.m>
% * + <matlab:openInEditor('timecount.m') timecount.m>
% * + <matlab:openInEditor('rmsd.m') rmsd.m>
%

leafWidth = 60; % 60 mm
timeStep = 0.001; % 1ms
stepAngleDeg = 1.8 / 8; % degree

% load the latest edition of DWELL_TIME_*.mat
latestFile = selatest('data', 'DWELL_TIME_*.mat');
D = load(latestFile);

% call the script for dividing splines into S-curves
useDefaultSettings = false; spdivRun

% append C to DWELL_TIME_*.mat
save(latestFile, 'C', '-append')
disp(['... append C to ' latestFile])

% Profile the trajectory optimization
% profile on

% call the script for generating the ion beam current density distribution
useDefaultSettings = false; dwelloptRun

% profile off
% profsave(profile('info'),'dwellopt_profile')


fileName = ['data/OPT_RESULT_' num2str(now, 12) '.mat'];
save(fileName, 'C', 'params', 'rmsds', 'optDwellTime')
disp(['... saved as ' fileName])

