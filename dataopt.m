function [pulseFreqs, pulseNums, nRamps, crtSteps] = dataopt(wa, wd, wf, w1, w2, crt, strokeTime)
%DATAOPT data table optimization

% by wulx, 2014/7/11, 2014/7/13, 2014/11/5


% variable shortcuts
scaleDivs = crt.scaleDivs;
locs = crt.locs;
types = crt.types;
crtSteps = crt.steps;

% steps difference should not be less than 4
stepDiffs = crtSteps(2:end) - crtSteps(1:end-1);
narrows = stepDiffs<4 | -stepDiffs<4 ;
broadenMask = [narrows(1); narrows(1:end)] | [narrows(1:end); narrows(end)];

% broaden differences of neighbors
amp1 = 10 * w1 * (types>0); % pull up peaks
amp2 = -10 * w2 * (types<0); % pull down valeys
amp3 = 2 * broadenMask .* types; % broaden the narrow gaps

crtSteps = round(crtSteps + amp1 + amp2 + amp3);

nRamps = numel(crt.steps) - 1;

pulseNums = cell(1, nRamps);
pulseFreqs = cell(1, nRamps);

timeDiff = 0;
for j = 1:nRamps
    if crtSteps(j+1)>crtSteps(j)
        nSteps = crtSteps(j+1) - crtSteps(j);
    else
        nSteps = crtSteps(j) - crtSteps(j+1);
    end
    
    % distance
    targDistance = scaleDivs(locs(j+1)) - scaleDivs(locs(j));
    strkLength = scaleDivs(end) - scaleDivs(1);
    targTime = strokeTime*targDistance/strkLength - timeDiff;
    
    % numbers of steps, [nAcc nConst nDec]
    nAcc = round(2 + wa(j)*(nSteps-4));
    nDec = round(2 + wd(j)*(nSteps-nAcc-2));
    nConst = nSteps - nAcc - nDec;

    % time per step algorithm
    [F, ~, timeDiff] = timeperstepd(nAcc, nConst, nDec, targTime, wf(j));

    freqs = [F(1:nAcc-1); ones(nConst+2,1)*F(nAcc); F(nAcc+1:end)];
    tseqs = 1 ./ freqs;
    
    %! using image processing techniques to reduce datatable
    freqs_rshift = [0; freqs(1:end-1)];
    
    neighbors = (freqs == freqs_rshift);
    
    [L, nLabels] = bwlabel(neighbors);
    
    L_lshift = [L(2:end); 0];
    
    labels = max(L, L_lshift);
    
    reducedFreqs = freqs;
    reducedTimeSeqs = tseqs;
    for i = 1:nLabels
        inds = find(labels == i);
        
        reducedFreqs(inds(2:end)) = nan;
        reducedTimeSeqs(inds(1)) = sum(tseqs(inds));
    end
    
    % clear all NaNs
    nans = isnan(reducedFreqs);
    reducedFreqs(nans) = [];
    reducedTimeSeqs(nans) = [];
    
    pulseFreqs{j} = reducedFreqs'; % in row
    
    if crtSteps(j+1)>crtSteps(j)
        pulseNums{j} = round(reducedFreqs .* reducedTimeSeqs)';
    else
        pulseNums{j} = -round(reducedFreqs .* reducedTimeSeqs)';
    end
end

