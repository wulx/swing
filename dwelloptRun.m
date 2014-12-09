%% Run optimiztion of dwell times

% (c) wulx, 2014/7/12, 2014/12/1

if useDefaultSettings
    stepAngleDeg = 1.8 / 8;
    leafWidth = 60; % 60 mm
    timeStep = 0.001; % 1ms
    
    % load the latest edition of DWELL_TIME_*.mat
    D = load(selatest('data/DWELL_TIME_*.mat'));
end


% Open pool of MATLAB sessions for parallelcomputation
isOpen = matlabpool('size') > 0;
if ~isOpen
    % fix bug, convert str to double firstly
    nProcessers = str2double(getenv('NUMBER_OF_PROCESSORS'));
    if nProcessers > 3
        matlabpool('local', nProcessers-1);
    else
        matlabpool('local', nProcessers);
    end
end

dwellTime = D.dwellTime;
strokeTime = D.strokeTime;
strkSet = D.strkSet;
% C = D.C;

% Note that
% leafWidth/vStroke = 2*maxEtchDepthContrast/(1+maxEtchDepthContrast)
dwellTimeUpperLimit = 2 * strokeTime * D.maxEtchDepthContrast / (1 + D.maxEtchDepthContrast);

% % divide splines into S-curves
% opts = {'MINPEAKHEIGHT', -0.004, 'MINPEAKDISTANCE', 15, 'NPEAKS', 5};
% C = spdiv(dwellTime, dwellTimeUpperLimit, strkSet, stepAngleDeg, opts);
% spdivplot(C) % plot spdiv results

nLayers = numel(strkSet);
nStrks = sum([strkSet.nStrks]);

params = cell(nStrks, 1);
rmsds = cell(nStrks, 1);
optDwellTime = nan( size(dwellTime) );

psopts = psoptimset('Cache', 'on', 'Vectorized','off', 'MaxIter', 20, ...
    'UseParallel', 'always', 'CompletePoll', 'on', 'TolFun', 0.01, ...
    'TolX', 0.01, 'PollMethod', 'GPSPositiveBasis2N', ...
    'MaxMeshSize', 0.1, 'InitialMeshSize', 0.1, 'MeshAccelerator', 'on', ...
    'Display', 'iter'); %, ...
%     'PlotFcns', {@psplotbestf,@psplotmeshsize,@psplotfuncount,@psplotbestx}, ...
%     'PlotInterval', 1);

for j = 1:nLayers
    for k = strkSet(j).indices
        
        xs = C(k).scaleDivs;
        
        figure, hold on;
        hSp = plot(xs, C(k).spline, 'k-');
        
        % mark all critical points at the spline
        for ni = 1:numel(C(k).locs)
            li = C(k).locs(ni);
            ti = C(k).types(ni);
            
            segmarker('make', [xs(li), C(k).spline(li)], ti, hSp);
        end
        
        
        nw = numel(C(k).steps) - 1;
        
        % weighted factors to be optimized
        w = [1/3*ones(1,2*nw), 0.4*ones(1,nw), 0.3, 0.3];
        
        % [a, b, c] = dwopt(params);
        lb = [zeros(1,3*nw), 0.1, 0.1];  % Lower bounds
        ub = [0.99*ones(1,3*nw), 0.9, 0.9];  % Upper bounds
        
        % params = [wa, wd, wf, w1, w2]
        dwopt = @(w) dwellopt(w(1:nw), w(nw+(1:nw)), w(2*nw+(1:nw)), w(end-1), w(end), ...
            C(k), strokeTime, timeStep, stepAngleDeg, leafWidth);
        
        % pattern search parallely
        [params{k}, rmsds{k}] = patternsearch(dwopt, w, [], [], [], [], lb, ub, psopts);
        
%         [params{k}, rmsds{k}] = particleswarm(dwopt, 3*nw+2, lb, ub);
%         [params{k}, rmsds{k}] = simulannealbnd(dwopt, w, lb, ub);
        
        % plot results
        [~, dwell1] = dwopt(params{k});
        
        xx = C(k).scaleDivs(2:end-1);

        figure, hold on;
        plot(xx, C(k).ogee, 'k-');
        plot(xx, dwell1, 'b-');
        title(['RMSD: ' num2str(rmsds{k})]);
        
        if strcmp(strkSet(j).strkDir, 'DOWN')
            optDwellTime(:, k) = dwell1(end:-1:1);
        else
            optDwellTime(:, k) = dwell1;
        end
    end
end

if isOpen
    matlabpool close
end
