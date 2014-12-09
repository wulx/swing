%% parameters optimization of the raster scanning algorithm

if useDefaultSettings
    load data/ETCH_DEPTH.mat % load the etch-depth map
    
    leafWidth = 60; % mm
    etchRate = 0.75; % provided that mean etch rate is 0.75 nm/s
    
    N = [2 4 6];
    W = [12 36 60];
end

nN = numel(N);
nW = numel(W);

R = nan(nN, nW);
for i = 1:nN
    for j = 1:nW
        R(i, j) = rascan(depth, N(i), W(j), leafWidth, etchRate);
    end
end

%% show figures
% # RMSD vs number of tiers
wStrs = cell(nW, 1);
wMarkers = 'osd'; % Circle, Square and Diamond

figure, hold on;
for m = 1:nW
    plot(N, R(:,m), 'LineStyle', '-', 'Marker', wMarkers(m), 'color', 'black', 'MarkerFaceColor', 'black')
    wStrs{m} = ['W_{beam} = ' num2str(W(m))];
end

legend(wStrs)

set(gca, 'XTick', N, 'XLim', [3*N(1)/4 N(end)+N(1)/4])
xlabel('Number of layers, 2n')
ylabel('RMSD[nm]')

% # RMSD vs ion beam width
nStrs = cell(nW, 1);

figure, hold on;
for n = 1:nN
    plot(W, R(n,:), 'LineStyle', '-', 'Marker', wMarkers(n), 'color', 'black', 'MarkerFaceColor', 'black')
    nStrs{n} = ['2n = ' num2str(N(n))];
end

legend(nStrs, 'Location', 'NorthWest')

set(gca, 'XTick', W, 'XLim', W(1)/2 + [0 W(end)])
xlabel('Ion beam width, W_{beam}[mm]')
ylabel('RMSD[nm]')
