function dwell2depth(strkSet, strokeTime, optDwellTime, runDwellTime)
%DWELL2DEPTH general convertor for transforming dwell times into etch
% depths

% (c) wulx, 2014/12/5




% init or preallocation
subSize = size(strkSet(1).data);

nanPad = nan(subSize);

supDepth = zeros(subSize);
optSupDepth = supDepth;
runSupDepth = supDepth;

for i = 1:numel(strkSet)
    optRibbons = meanEtchRate * (strokeTime - optDwellTime(:, strkSet(i).indices));
    runRibbons = meanEtchRate * (strokeTime - runDwellTime(:, strkSet(i).indices));
    
    % padding ------------------------------------------------------------%
    preLeft = round([strkSet(i).padding(1), strkSet(i).padding(4)]);
    prePad = padarray(nanPad, preLeft, nan, 'pre');
    postRight = round([strkSet(i).padding(3), strkSet(i).padding(2)]);
    postPad = padarray(prePad, postRight, nan, 'post');
    
    % affine transform ---------------------------------------------------%
    xform = [1, 0, 0; strkSet(i).shear, 1, 0; 0, 0, 1];
    tform = maketform('affine', xform);
    shPad = imtransform(postPad, tform, 'nearest', 'FillValues', nan);
    
    % divide sheared padding area into stroke strips ---------------------%
    [shHeight, shWidth] = size(shPad);
    startIdx = 1:ionBeamWidth:shWidth;
    endIdx = [startIdx(2:end)-1, shWidth];
    
    optRasterMap = nan(shHeight, shWidth);
    runRasterMap = optRasterMap;
    
    for n = 1:strkSet(i).nStrks
        for r = 1:shHeight
            optRasterMap(r, startIdx(n):endIdx(n)) = optRibbons(r, n);
            runRasterMap(r, startIdx(n):endIdx(n)) = runRibbons(r, n);
        end
    end
    
    xform2 = [1, 0, 0; -strkSet(i).shear, 1, 0; 0, 0, 1];
    tform2 = maketform('affine', xform2);
    
    xdata = strkSet(i).padding(4) + [1 subWidth];
    if strcmp(strkSet(i).strkDir, 'UP')
        xdata = xdata + halfPeriod;
    end
    
    ydata = strkSet(i).padding(1) + [1 subHeight];
    
    optData = imtransform(optRasterMap, tform2, 'nearest', 'FillValues', nan, ...
        'XData', xdata, 'YData', ydata);
    runData = imtransform(runRasterMap, tform2, 'nearest', 'FillValues', nan, ...
        'XData', xdata, 'YData', ydata);
    
    supDepth = supDepth + strkSet(i).data;
    optSupDepth = optSupDepth + optData;
    runSupDepth = runSupDepth + runData;
    
end
