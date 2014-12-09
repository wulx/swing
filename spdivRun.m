
close all; clear all; clc

useDefaultSettings = true;
if useDefaultSettings
    stepAngleDeg = 1.8 / 8;
    leafWidth = 60; % 60 mm
    timeStep = 0.001; % 1ms
    
    % load the latest edition of DWELL_TIME_*.mat
    latestFile = selatest('data', 'DWELL_TIME_*.mat');
    D = load(latestFile);
end

dwellTime = D.dwellTime;
strokeTime = D.strokeTime;
strkSet = D.strkSet;
% Note that
% leafWidth/vStroke = 2*maxEtchDepthContrast/(1+maxEtchDepthContrast)
dwellTimeUpperLimit = 2 * strokeTime * D.maxEtchDepthContrast / (1 + D.maxEtchDepthContrast);

% divide splines into S-curves
opts = {'MINPEAKHEIGHT', -0.001, 'MINPEAKDISTANCE', 30, 'NPEAKS', 5};
C = spdiv(dwellTime, dwellTimeUpperLimit, strkSet, stepAngleDeg, opts);

% % append C to DWELL_TIME_*.mat
% save(latestFile, 'C', '-append')
% disp(['... append C to ' latestFile])


% %SPDIVPLOT plot spdiv results
% % by wulx, 2014/7/12
% 
% nStrks = sum([strkSet.nStrks]);
% 
% headsAndTails = nan(2, nStrks);
% 
% for j = 1:numel(C)
%     x = C(j).scaleDivs;
%     y = C(j).ogee;
%     sp = C(j).spline;
% 
%     hFig = figure;
%     hold on;
%     hOg = plot(x(2:end-1), y, 'k:');
%     hSp = plot(x, sp, 'k-');
% 
%     % mark all critical points at the spline
%     for ni = 1:numel(C(j).locs)
%         li = C(j).locs(ni);
%         ti = C(j).types(ni);
%         
%         segmarker('make', [x(li), sp(li)], ti, hSp)
%     end
%     
%     segmarker on
%     hCtrl = uicontrol('Position',[20 20 80 20],'String','Continue',...
%         'Callback','uiresume(gcbf)');
%     uiwait(gcf);
%     
%     hSeg = findobj(gcf, '-regexp', ...
%         'Tag', 'segmarker_(foothill|valley|peak)');
%     
%     segmarker off
%     close(hFig);
%     
%     crtDwellTimes = sp(C(ind).locs);
%     crtDegs = rad2deg( asin(crtDwellTimes / dwellTimeUpperLimit) );
%     C(ind).steps = round(crtDegs / stepAngleDeg); % estimates
%     
%     headsAndTails(:, ind) = C(ind).steps([1 end]);
% 
% end
% 
% %#!IMPORTANT
% % average and splice heads and tails
% headsAndTails(1, 2:end) = round((headsAndTails(1, 2:end) + headsAndTails(2, 1:end-1)) / 2);
% headsAndTails(2, 1:end-1) = headsAndTails(1, 2:end);
% 
% for j = 1:nStrks
%     % revise the heads and tails
%     C(j).steps([1 end]) = headsAndTails(:, j);
% end

