function outDwellTime = data2dwell(strokeTime, strkSet, C, DataTable, leafWidth, stepAngleDeg, timeStep)
%DATA2DWELL reversely transform or reflect the data table into dwell times

% (c) wulx, 2014/12/5
% 


nScan = ceil(strokeTime / timeStep);

steps = nan(1, nScan);
timeline = linspace(0, strokeTime, nScan);

nsum = @(n, list) sum(list(1:n));

outDwellTime = nan(numel(C(1).ogee), numel(C));

for j = 1:numel(strkSet)
    for k = strkSet(j).indices
        pulseNums = [DataTable(k).pulseNums{:}];
        pulseFreqs = [DataTable(k).pulseFreqs{:}];
        
        timeSeqs = abs(pulseNums) ./ pulseFreqs;
        sumSteps = DataTable(k).crtSteps(1) + arrayfun(@(n) nsum(n, pulseNums), 1:numel(pulseNums));
        
        nStair = numel( pulseNums );
        for i = 1:(nStair-1)
            inds = timeline>=nsum(i-1, timeSeqs) & timeline<=nsum(i, timeSeqs);
            steps(inds) = sumSteps(i);
        end
        steps(timeline>=nsum(nStair-1, timeSeqs)) = sumSteps(nStair);
        
        projWidths = step2width(steps, stepAngleDeg, leafWidth);
        
        dwTime = timecount(projWidths, strokeTime, C(k).scaleDivs);
        dwTime = dwTime(2:end-1); % remove boundaries

        mask = isnan(C(k).ogee);
        
        if strcmp(strkSet(j).strkDir, 'DOWN')
            dwTime = dwTime(end:-1:1);
            mask = mask(end:-1:1);
        end

        dwTime(mask) = nan;
        
        % figure, plot(C(k).scaleDivs(2:end-1), dwTime);
        
        outDwellTime(:, k) = dwTime;
    end
end

