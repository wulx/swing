function outEtchDepth = dwell2depth(strkSet, inDwellTime, strokeTime, etchRate, ionBeamWidth)
%DWELL2DEPTH general convertor for transforming dwell times into etch
% depths

% (c) wulx, 2014/12/5


% init or preallocation
subSize = size(strkSet(1).data);

outEtchDepth = zeros(subSize);

nanPad = nan(subSize);

nLayers = numel(strkSet);
halfPeriod = ionBeamWidth / nLayers;

for i = 1:nLayers
    if isscalar(etchRate)
        ribbons = etchRate * (strokeTime - inDwellTime(:, strkSet(i).indices)); % etch rate can be a scalar
    else
        ribbons = repmat(etchRate, 1, strkSet(i).nStrks) .* ...
            (strokeTime - inDwellTime(:, strkSet(i).indices)); % or a vector
    end
    
    % padding --------------------------------------------------------%
    preLeft = round([strkSet(i).padding(1), strkSet(i).padding(4)]);
    prePad = padarray(nanPad, preLeft, nan, 'pre');
    postRight = round([strkSet(i).padding(3), strkSet(i).padding(2)]);
    postPad = padarray(prePad, postRight, nan, 'post');
    
    % affine transform -----------------------------------------------%
    xform = [1, 0, 0; strkSet(i).shear, 1, 0; 0, 0, 1];
    tform = maketform('affine', xform);
    shPad = imtransform(postPad, tform, 'nearest', 'FillValues', nan);
    
    % divide sheared padding area into stroke strips -----------------%
    [shHeight, shWidth] = size(shPad);
    startIdx = 1:ionBeamWidth:shWidth;
    endIdx = [startIdx(2:end)-1, shWidth];
    
    rasterMap = nan(shHeight, shWidth);
    
    for n = 1:strkSet(i).nStrks
        for r = 1:shHeight
            rasterMap(r, startIdx(n):endIdx(n)) = ribbons(r, n);
        end
    end
    
    xform2 = [1, 0, 0; -strkSet(i).shear, 1, 0; 0, 0, 1];
    tform2 = maketform('affine', xform2);
    
    xdata = strkSet(i).padding(4) + [1 subSize(2)];
    if strcmp(strkSet(i).strkDir, 'UP')
        xdata = xdata + halfPeriod;
    end
    
    ydata = strkSet(i).padding(1) + [1 subSize(1)];
    
    outData = imtransform(rasterMap, tform2, 'nearest', 'FillValues', nan, ...
        'XData', xdata, 'YData', ydata);
    
    outEtchDepth = outEtchDepth + outData;
end
