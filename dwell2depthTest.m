%% transform dwell time into etch depth
%
close all; clear all; clc

n = 1; % microsteps per full step, number of periods per ionBeamWidth
nTiers = 2*n;
ionBeamWidth = 36; % mm @ionBeamWidth ------------------------------------@
ionBeamWidth = round(ionBeamWidth/nTiers) * nTiers;
disp(['Ion beam width is revised as ' num2str(ionBeamWidth)])

leafWidth = 60; % mm
subHeight = 400; % mm
subWidth = 400; % mm
vStroke = leafWidth + subHeight; % vertical stroke

halfPeriod = ionBeamWidth / nTiers;

% provided that mean etch rate is 1 nm/s
meanEtchRate = 0.75;

%%
load('data/DWELL_TIME_236_735924.806031.mat', 'strokeTime', 'strkSet');


load('data/OPT_RESULT_735924.878627.mat', 'rmses');

S3 = load('data/ETCH_DEPTH.mat', 'depth');
etchDepth = S3.depth(end:-1:1, end:-1:1);

load('data/ERROR_MODEL_735924.945739.mat', 'dwellTime', 'optDwellTime', 'runDwellTime');


nTiers = numel(strkSet);
supDepth = zeros(subHeight, subWidth);
optSupDepth = supDepth;
runSupDepth = supDepth;

nanPad = nan(size(strkSet(1).data));

for i = 1:nTiers
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

% meanDepth = mean(etchDepth(:));
rmse1 = rmse(etchDepth, supDepth);
rmse2 = rmse(etchDepth, optSupDepth);
rmse3 = rmse(etchDepth, runSupDepth);

figure

subplot(2, 2, 1)
imshow( mat2gray(etchDepth) )
colormap jet
set(gca, 'XDir', 'reverse', 'YDir', 'normal', 'Visible', 'off', 'YAxisLocation', 'right')

subplot(2, 2, 2)
imshow( mat2gray(supDepth) )
colormap jet
set(gca, 'XDir', 'reverse', 'YDir', 'normal', 'Visible', 'off', 'YAxisLocation', 'right')
title(['sup RMSE: ' num2str(rmse1)])

subplot(2, 2, 3)
imshow( mat2gray(optSupDepth) )
colormap jet
set(gca, 'XDir', 'reverse', 'YDir', 'normal',  'YAxisLocation', 'right')
title(['opt RMSE: ' num2str(rmse2)])

subplot(2, 2, 4)
imshow( mat2gray(runSupDepth) )
colormap jet
set(gca, 'XDir', 'reverse', 'YDir', 'normal',  'YAxisLocation', 'right')
title(['run RMSE: ' num2str(rmse3)])

dwellTimeWithoutNan = dwellTime(~isnan(dwellTime));
% normalizedRmses = [rmses{:}] / (max(dwellTimeWithoutNan) - min(dwellTimeWithoutNan));

disp('RMSD: (s)')
disp(num2str([rmses{1:2:end}], 3))
disp(num2str([rmses{2:2:end}], 3))


runDiffDepth = runSupDepth - etchDepth;
figure,
hf = histfit(runDiffDepth(:), 51, 'normal');


% h = findobj(gca,'Type','patch');
set(hf(1),'FaceColor','k','EdgeColor','w')
xlim([-0.4 0.4])
xlabel('Depth difference (nm)')
ylabel('Frequency')

[muhat,sigmahat,muci,sigmaci] = normfit(runDiffDepth(:),alpha);

% [D, PD] = allfitdist(runDiffDepth(:),'PDF');

% save('data/ERROR_MODEL.mat', 'supDepth', 'optSupDepth', 'runSupDepth', '-append')
