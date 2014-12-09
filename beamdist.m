function [exErrMap, errMap, mus] = beamdist(nMajor, nMinor, nMajorRands, nMinorSamps, lambda, nsig)
%BEAMDIST ion beam current density distribution

% (c) wulx, 2014/11/30 

if nargin<3, nsig = 0.6; end

% uncerntainty: +/-5% (lambda)
% ref.: Zhou et al._2007_Amending the uniformity of ion beam current density profile
% ion beam current density distribution along major axis

mus = 1.0 + lambda * (rand(nMajorRands, 1) - 0.5);

coffs = 0.01 * (rand(nMajorRands, 1) - 0.5);
nSigmas = nsig + 0.01 * (rand(nMajorRands, 1) - 0.5);

errMap = nan(nMajorRands, nMinorSamps);

% ion beam current density distribution along minor axis
% standard normal distribution: mu = 0, sigma = 1
stdNormDist = makedist('Normal');

for i = 1:nMajorRands
    mu = mus(i);
    nSigma = nSigmas(i);
    coff = coffs(i); % center offset, range (-1, 1)
    
    % sampling points
    sampPts = linspace((coff-1)*nSigma, (coff+1)*nSigma, nMinorSamps);
    
    samps = pdf(stdNormDist, sampPts);
    
    errMap(i, :) = mu/mean(samps) * samps;
end

[X, Y] = meshgrid(1:nMinorSamps, 1:nMajorRands);
[XI, YI] = meshgrid(linspace(1, nMinorSamps, nMinor), linspace(1, nMajorRands, nMajor));

exErrMap = interp2(X, Y, errMap, XI, YI, 'spline');
