%% Run raster scanning 

if useDefaultSettings
    load data/ETCH_DEPTH.mat % load etch-depth map
    
    nTiers = 2; % the number of tiers should be even, such as 2, 4, 6 ...
    ionBeamWidth = 36; % mm
    leafWidth = 60; % mm
    etchRate = 0.75; % provided that mean etch rate is 0.75 nm/s
    showImage = true;
end

x1 = 0.5:399.5; % mm, 400 points in total
y1 = 0.5:399.5; % mm
z1 = depth; % in nm

figure('Name','etch depth map (mesh)')
mesh(x1, y1, z1)
axis ij
view([-32.5000, 75.0000])

[r, dwellTime, strokeTime, maxTunableRatio, strkSet] = rascan(depth, nTiers, ionBeamWidth, leafWidth, etchRate, showImage);

fileName = ['data/DWELL_TIME_' num2str(nTiers) '_' num2str(ionBeamWidth) '_' num2str(now, 12) '.mat'];
save(fileName, 'dwellTime', 'strokeTime', 'maxTunableRatio', 'strkSet')
disp(['... saved as ' fileName])