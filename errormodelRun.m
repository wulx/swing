%% run error model

% close all; clear all; clc

if useDefaultSettings
    % load caculated data
    load('data/ETCH_DEPTH.mat', 'depth');
    load(selatest('data', 'DWELL_TIME_*.mat'), 'strokeTime',  'strkSet', 'dwellTime');
    load(selatest('data', 'OPT_RESULT_*.mat'), 'C', 'optDwellTime', 'rmsds');
    load(selatest('data', 'DATA_TABLE_*.mat'), 'DataTable');
    
    leafWidth = 60; % 60 mm
    stepAngleDeg = 1.8 / 8;
    timeStep = 0.001; % 1 ms
    
    ionBeamWidth = 36; % mm
    etchRate = 0.75;
end

outDwellTime = data2dwell(strokeTime, strkSet, C, DataTable, leafWidth, stepAngleDeg, timeStep);

nLayers = numel(strkSet);
ionBeamWidth = round(ionBeamWidth/nLayers) * nLayers;
disp(['Ion beam width: ' num2str(ionBeamWidth) ', Number of layers: ' num2str(nLayers)])

halfPeriod = ionBeamWidth / nLayers;

subSize = size(strkSet(1).data);
rasEtchDepth = zeros(subSize);

for i = 1:nLayers
    rasEtchDepth = rasEtchDepth + strkSet(i).data;
end

% outEtchDepth = dwell2depth(strkSet, inDwellTime, strokeTime, etchRate, ionBeamWidth);
optEtchDepth = dwell2depth(strkSet, optDwellTime, strokeTime, etchRate, ionBeamWidth);
outEtchDepth = dwell2depth(strkSet, outDwellTime, strokeTime, etchRate, ionBeamWidth);

%%
etchDepth = depth(end:-1:1, end:-1:1);

% meanDepth = mean(etchDepth(:));
rmsd1 = rmsd(etchDepth, rasEtchDepth);
rmsd2 = rmsd(etchDepth, optEtchDepth);
rmsd3 = rmsd(etchDepth, outEtchDepth);

figure

subplot(2, 2, 1)
imshow( mat2gray(etchDepth) )
colormap jet
set(gca, 'XDir', 'reverse', 'YDir', 'normal', 'Visible', 'off', 'YAxisLocation', 'right')

subplot(2, 2, 2)
imshow( mat2gray(rasEtchDepth) )
colormap jet
set(gca, 'XDir', 'reverse', 'YDir', 'normal', 'Visible', 'off', 'YAxisLocation', 'right')
title(['Rastered etch depth, RMSD: ' num2str(rmsd1)])

subplot(2, 2, 3)
imshow( mat2gray(optEtchDepth) )
colormap jet
set(gca, 'XDir', 'reverse', 'YDir', 'normal',  'YAxisLocation', 'right')
title(['Optimized etch depth, RMSD: ' num2str(rmsd2)])

subplot(2, 2, 4)
imshow( mat2gray(outEtchDepth) )
colormap jet
set(gca, 'XDir', 'reverse', 'YDir', 'normal',  'YAxisLocation', 'right')
title(['Output etch depth, RMSD: ' num2str(rmsd3)])

%%
save(['data/ERROR_MODEL_' num2str(now, 12) '.mat'], ...
    'dwellTime', 'optDwellTime', 'outDwellTime', ...
    'rasEtchDepth', 'optEtchDepth', 'outEtchDepth');

%%
disp('RMSD: (s)')
disp(num2str([rmsds{1:2:end}], 3))
disp(num2str([rmsds{2:2:end}], 3))

%%
outDiffDepth = outEtchDepth - etchDepth;

figure,
hf = histfit(outDiffDepth(:), 51, 'normal');

% h = findobj(gca,'Type','patch');
set(hf(1), 'FaceColor', 'k', 'EdgeColor', 'w')
xlim([-0.4 0.4])
xlabel('Depth difference (nm)')
ylabel('Frequency')

alpha = 0.05; % 100(1-alpha)%
[muhat, sigmahat, muci, sigmaci] = normfit(outDiffDepth(:), alpha);

