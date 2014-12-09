%% generate data table

% by wulx, 2014/7/12, 2014/7/13, 2014/12/5


if useDefaultSettings
    % load stroke time
    latestDwellFile = selatest('data', 'DWELL_TIME_*.mat');
    load(latestDwellFile, 'strokeTime');
    
    % load optimization results
    latestOptFile = selatest('data', 'OPT_RESULT_*.mat');
    load(latestOptFile, 'C', 'params');
end

%% datatable

% use struct array to save objective data
nStrks = numel(params);
emptyCell = cell(1, nStrks);
DataTable = struct('pulseFreqs', emptyCell, ...
    'pulseNums', emptyCell, ...
    'nRamps', emptyCell, ...
    'crtSteps', emptyCell);

for  n = 1:nStrks

    % calculate frequencies and time sequencies
    crt = C(n);
    nRamps = numel(C(n).steps) - 1;
    prm = params{n};
    wa = prm(1:nRamps);
    wd = prm(nRamps+(1:nRamps));
    wf = prm(2*nRamps+(1:nRamps));
    w1 = prm(end-1);
    w2 = prm(end);

    [pulseFreqs, pulseNums, nRamps, crtSteps] = dataopt(wa, wd, wf, w1, w2, crt, strokeTime);
    
    DataTable(n).pulseFreqs = pulseFreqs;
    DataTable(n).pulseNums = pulseNums;
    DataTable(n).nRamps = nRamps;
    DataTable(n).crtSteps = crtSteps;

end

%% datatable2 (reverse datatable)

DataTable2 = DataTable(end:-1:1);

for m = 1:nStrks
    strk = DataTable2(m);
    
    % reverse critical steps
    DataTable2(m).crtSteps = strk.crtSteps(end:-1:1);

    for i = 1:strk.nRamps
        % reverse pulse numbers
        pulseNum = strk.pulseNums{end+1-i};
        DataTable2(m).pulseNums{i} = -pulseNum(end:-1:1);
        
        % reverse pulse frequencies
        pulseFreq = strk.pulseFreqs{end+1-i};
        DataTable2(m).pulseFreqs{i} = pulseFreq(end:-1:1);
    end

end

