function C = spdiv(dwellTime, dwellTimeUpperLimit, strkSet, stepAngleDeg, opts)
%SPDIV DIVide the SPlines into S curves and extract critical points of the S curves

% by wulx, 2014/7/10

if nargin < 5
    opts = {'MINPEAKHEIGHT', -1, 'MINPEAKDISTANCE', 15, 'NPEAKS', 5};
end

% struct to store critical points
C = struct('ogee', {}, ...
    'spline', {}, ...
    'scaleDivs', {}, ...
    'locs', {}, ...
    'types', {}, ...
    'steps', {});

nLayers = numel(strkSet);
nStrks = sum([strkSet.nStrks]);

headsAndTails = nan(2, nStrks);

for i = 1:nLayers
    for ind = strkSet(i).indices
        % reverse processing directions for down stroks
        if strcmp(strkSet(i).strkDir, 'DOWN')
            C(ind).ogee = dwellTime(end:-1:1, ind);
        else
            C(ind).ogee = dwellTime(:, ind);
        end
        
        % spline interpolation (or extrapolation) and smoothing
        [xs, ys] = spis(C(ind).ogee);
               
        % low pass filtering
        lowPass = ys>dwellTimeUpperLimit;
        if any(lowPass)
            ys(lowPass) = dwellTimeUpperLimit;
        end
        
        % high pass filtering
        highPass = ys<0.01;
        if any(highPass)
            ys(highPass) = 0.01;
        end
        
        C(ind).scaleDivs = xs;
        C(ind).spline = ys;
        
        % find critical points
        [locs, ptype] = findpoints(ys, opts);
        
        C(ind).locs = [1; locs; numel(ys)]; % all critical locations
        C(ind).types = [0; ptype; 0]; % treate boundary points as infletion points
      
        hFig = figure;
        hold on;
        plot(xs(2:end-1), C(ind).ogee, 'k:')
        hSp = plot(xs, ys, 'k-');
        
        % mark all critical points at the spline
        for ni = 1:numel(C(ind).locs)
            li = C(ind).locs(ni);
            ti = C(ind).types(ni);
            
            segmarker('make', [xs(li), ys(li)], ti, hSp);
        end

        segmarker on
        uicontrol('Position',[5 5 80 30],'String','Continue',...
            'Callback','uiresume(gcbf)')
        uiwait(gcf)
        
        hSegs = findobj(gcf, '-regexp', ...
            'Tag', '^segmarker_(foothill|valley|peak)$');
        
        nSeg = numel(hSegs);
        
        xd = cell2mat(get(hSegs, 'XData'));
        yd = cell2mat(get(hSegs, 'YData'));
        [xb, xi] = sort(xd);
        yd = yd(xi);
        hSegs = hSegs(xi);
        
        [~, ~, iSegs] = unique([xb; xs]);
        C(ind).locs = iSegs(1:nSeg);
        
        types = nan(nSeg, 1);
        for ns = 1:nSeg
            tagName = get(hSegs(ns), 'Tag');
            
            switch tagName
                case 'segmarker_foothill'
                    sType = 0;
                case 'segmarker_peak'
                    sType = 1;
                case 'segmarker_valley'
                    sType = -1;
            end
            
            types(ns) = sType;
        end
        C(ind).types = types;
        
        segmarker off
        close(hFig);
        
        crtDwellTimes = yd; %ys(C(ind).locs);
        crtDegs = rad2deg( asin(crtDwellTimes / dwellTimeUpperLimit) );
        C(ind).steps = round(crtDegs / stepAngleDeg); % estimates
        
        headsAndTails(:, ind) = C(ind).steps([1 end]);
    end
end

%#!IMPORTANT
% average and splice heads and tails
headsAndTails(1, 2:end) = round((headsAndTails(1, 2:end) + headsAndTails(2, 1:end-1)) / 2);
headsAndTails(2, 1:end-1) = headsAndTails(1, 2:end);

for j = 1:nStrks
    % revise the heads and tails
    C(j).steps([1 end]) = headsAndTails(:, j);
end

