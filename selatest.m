function fileLatest = selatest(fileDir, fileFormat)
%GETLATEST select the latest edition from the files with similar names in a
% specified directory.

F = dir(fullfile(fileDir, fileFormat));
fileNames = { F.name };

[~, ix] = sort([F.datenum]);
fileLatest = fullfile(fileDir, fileNames{ix(end)});
