%% Demo on generating the ion beam current density distribution

% (c) wulx, 2014/11/30, 2014/12/1

if useDefaultSettings
    vStroke = 460; % vertical stroke (mm)
    ionBeamWidth = 36; % mm 
    nMajorRands = 32;
    nMinorSamps = 8;
    lambda = 10/100; % +/-5% deviations
    nSig = 0.6 * ionBeamWidth/60;
end

[exErrMap, errMap, mus] = beamdist(vStroke, ionBeamWidth, nMajorRands, nMinorSamps, lambda, nSig);
xx = linspace(1, nMajorRands, vStroke);
yy = spline(1:nMajorRands, mus, xx);

figure, hold on;
xlim([1 nMajorRands])
ylim([0 1+2*lambda])
plot(xx, yy, 'k-')
stem(mus, 'k:')

plot(xlim, [1 1]-0.5*lambda, 'r:')
plot(xlim, [1 1]+0.5*lambda, 'r:')

% figure, hold on;
% surf(errMap)
% axis tight

figure, hold on;
meshz(exErrMap)
shading flat
axis equal

