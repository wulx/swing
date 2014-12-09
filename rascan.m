function [r, dwellTime, strokeTime, maxEtchDepthContrast, strkSet] = rascan(depth, nTiers, ionBeamWidth, leafWidth, etchRate, showImage)
%RASCAN RASter SCANning algorithm

% by wulx, 2014/7/6, 2014/7/7, 2014/7/10, 2014/11/30

if nargin < 6, showImage = false; end

% dimensions of the substrate, reference unit (default: mm)
[subHeight, subWidth] = size(depth);

% horizontally and vertically flipped depth map
depthFlipped = depth(end:-1:1, end:-1:1); % default in nm

% the ion beam width should be integer multiples of the number of tiers
halfPeriod = round(ionBeamWidth/nTiers);
ionBeamWidth = nTiers * halfPeriod;

vStroke = leafWidth + subHeight; % vertical stroke

% stroke sets
%   data: swept depth map
%   padding: [top right bottom left]'
%   indices: stroke No. indices
%   nStrks: number of strokes
%   ribbons: etch time profile
%   shear: shear along X-axis
%     [ 1   0   0
%      shx  1   0
%       0   0   1 ];
%   strkDir: bottom-UP or top-DOWN
%  
strkSet = struct('data', {}, ...
    'padding', {}, ...
    'indices', {}, ...
    'nStrks', {}, ...
    'ribbons', {}, ...
    'shear', {}, ...
    'strkDir', {});

avgDepth = depthFlipped / nTiers;
shx = halfPeriod / vStroke;

% strokes: 1..nTiers, nTiers+1..nTotalStrks
%   start from nTiers, the width of max-tier-area increases halfPeriod with
%   every new stroke adding.
nTotalStrks = nTiers + ceil(subWidth / halfPeriod);

% #1. grouping
nPairs = round(nTiers/2);
for i = 1:nPairs
    strkSet(2*i-1).strkDir = 'UP';
    strkSet(2*i-1).shear = -shx;

    strkSet(2*i-1).indices = (2*i-1):nTiers:nTotalStrks;
    strkSet(2*i-1).nStrks = numel(strkSet(2*i-1).indices);
    
    leftPadding = (nTiers - (2*i-1)) * halfPeriod;
    rightPadding = strkSet(2*i-1).nStrks*ionBeamWidth - halfPeriod - subWidth  - leftPadding;

    strkSet(2*i-1).padding = [
        0.5*leafWidth;
        rightPadding;
        0.5*leafWidth;
        leftPadding];
    
    strkSet(2*i).strkDir = 'DOWN';
    strkSet(2*i).shear = shx;

    strkSet(2*i).indices = (2*i):nTiers:nTotalStrks;
    strkSet(2*i).nStrks =  numel(strkSet(2*i).indices);
    
    leftPadding2 = (nTiers - 2*i) * halfPeriod;
    rightPadding2 = strkSet(2*i).nStrks*ionBeamWidth - halfPeriod - subWidth - leftPadding2;
    
    strkSet(2*i).padding = [
        0.5*leafWidth;
        rightPadding2;
        0.5*leafWidth;
        leftPadding2];
end

depthSup = zeros(subHeight, subWidth); % superposed depth map

for j = 1:nTiers
    % padding
    preLeft = round([strkSet(j).padding(1), strkSet(j).padding(4)]);
    prePad = padarray(avgDepth, preLeft, nan, 'pre');
    postRight = round([strkSet(j).padding(3), strkSet(j).padding(2)]);
    postPad = padarray(prePad, postRight, nan, 'post');
    
    % affine transform
    xform = [1, 0, 0; strkSet(j).shear, 1, 0; 0, 0, 1];
    tform = maketform('affine', xform);
    %tform = affine2d(xform);
    shPad = imtransform(postPad, tform, 'nearest', 'FillValues', nan);
    %shPad = imwarp(postPad, tform, 'nearest', 'FillValues', nan);
    
    % divide sheared padding area into ribbons
    [shHeight, shWidth] = size(shPad);
    startIdx = 1:ionBeamWidth:shWidth;
    endIdx = [startIdx(2:end)-1, shWidth];
    
    rasterStrips = nan(shHeight, strkSet(j).nStrks);
    rasterMap = nan(shHeight, shWidth);

    % #2. averaging
    for n = 1:strkSet(j).nStrks
        colIdx = startIdx(n):endIdx(n);
        for r = 1:shHeight
            row_r = shPad(r, colIdx);
            rowIdx = ~isnan( row_r );
            
            if any(rowIdx) > 0 % any() is better than sum() logically
                rasterStrips(r, n) = mean( row_r(rowIdx) );
                rasterMap(r, colIdx(rowIdx)) = rasterStrips(r, n);
            end
        end
    end
    
    if isscalar(etchRate)
        strkSet(j).ribbons = rasterStrips / etchRate;
    elseif isvector(etchRate) % etch time profile along the vertical
        strkSet(j).ribbons = cell2mat(arrayfun(@(n) rasterStrips(:,n)./etchRate(:), 1:strkSet(j).nStrks, 'UniformOutput', false));
    end
    
    xform2 = [1, 0, 0; -strkSet(j).shear, 1, 0; 0, 0, 1];
    tform2 = maketform('affine', xform2);
    %tform2 = affine2d(xform2);
    
    xdata = strkSet(j).padding(4) + [1 subWidth];
    if strcmp(strkSet(j).strkDir, 'UP')
        xdata = xdata + halfPeriod;
    end
    
    ydata = strkSet(j).padding(1) + [1 subHeight];
    
    strkSet(j).data = imtransform(rasterMap, tform2, 'nearest', 'FillValues', nan, ...
        'XData', xdata, 'YData', ydata);
    %strkSet(j).data = imwarp(rasterMap, tform2, 'nearest', 'FillValues', nan, ...
    %    'XData', xdata, 'YData', ydata);
     
    % #3. superposing
    depthSup = depthSup + strkSet(j).data;
end

% RMSD of the flipped depth map and the superposed depth map
r = rmsd(depthFlipped, depthSup);

fprintf('Ion beam width: %d mm, Number of tiers: %d, RMSD: %f\n', ionBeamWidth, nTiers, r)

% convert depth map to time profiles -------------------------------------%
ribbons = [strkSet.ribbons];
ribbonsWithoutNan = ribbons(~isnan(ribbons));

meanEtchTimeG = mean(ribbonsWithoutNan); % global mean etch time
% medEtchTimeG = median(stripNums);
maxEtchTimeG = max(ribbonsWithoutNan);
minEtchTimeG = min(ribbonsWithoutNan);

% define the maximal etch depth contrast
maxEtchDepthContrast = leafWidth / (2*vStroke - leafWidth);

% the necessary and sufficient condition for fully realizing the fine
% adjustment of etch depths.
% P.S. The step-by-step deduction is well described on the notebook.
etchTimeBase = (minEtchTimeG + maxEtchTimeG) / 2;
if maxEtchDepthContrast > (maxEtchTimeG-minEtchTimeG)/(maxEtchTimeG+minEtchTimeG)
    %etchTimeBase = (minEtchTimeG + maxEtchTimeG) / 2;
    
    % the additional condition for taking the meaETchTimeG as the etch time
    % baseline.
    if ((meanEtchTimeG > maxEtchTimeG/(1+maxEtchDepthContrast)) && (meanEtchTimeG<=minEtchTimeG/(1-maxEtchDepthContrast)))
        etchTimeBase = meanEtchTimeG;
        disp('The mean etch time is selected as the baseline of etch times')
    end
else
    warning('The range of etch times cannot be fully covered.')
end

strokeTime = (1 + maxEtchDepthContrast) * etchTimeBase;

% convert etch time map into dwell time map ------------------------------%
% dwellTime = strokeTime - etchTime

% maximium tunable dwell time
maxDwellTime = strokeTime * leafWidth / vStroke;

dwellTime = nan(size(ribbons));
for k = 1:nTiers
    dwellTime(:, strkSet(k).indices) = strokeTime - strkSet(k).ribbons;
end


if showImage
    figure('Name','flipped etch depth map')
    imshow( mat2gray(depthFlipped) )
    set(gca, 'XDir', 'reverse', 'YDir', 'normal', 'Visible', 'on', 'YAxisLocation','right')
    colormap jet

    for k = 1:nTiers
        figure, imshow(mat2gray(strkSet(k).data));
        title([strkSet(k).strkDir ' STROKES'])
        set(gca, 'XDir', 'reverse', 'YDir', 'normal', 'Visible', 'on', 'YAxisLocation', 'right')
        colormap jet
    end
    
    figure('Name', ['superposed depth map (RMSD: ' num2str(r) ')'])
    imshow( mat2gray(depthSup) )
    set(gca, 'XDir', 'reverse', 'YDir', 'normal', 'Visible', 'on', 'YAxisLocation', 'right')
    colormap jet
    
    figure, hold on;
    plot(ribbons)
    title('opt stroke time')
    
    plot(etchTimeBase*ones(vStroke, 1), 'k--')
    
    plot(strokeTime*ones(vStroke, 1), 'r-')
    plot((2*etchTimeBase-strokeTime)*ones(vStroke, 1), 'b:')
    
    axis tight
    
    figure, hold on;
    plot(dwellTime)
    plot(maxDwellTime*ones(vStroke, 1), 'r:')
    title('dwell time of all strokes')
    
    axis tight
end
