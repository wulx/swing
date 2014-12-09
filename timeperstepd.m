function [F, T, timeDiff] = timeperstepd(nAcc, nConst, nDec, targTime, wi, s)
%TIMEPERSTEP time per step algorithm for DWELL2

% by wulx, 2014/7/11, 2014/7/13

% narginchk(5, 6)
% climb only one stair every one microstep (default)
if nargin < 6
    s = 1;
end

if any([nAcc, nConst, nDec] < [2 0 2])
    error('timeperstep:nAcc_nConst_nDec', ...
        'The number of steps is a non-negative integer. For acc-/dec-eleration ramps, that should be no less than 2.');
end

if wi<0 || wi>=1
    error('timeperstep:wi', 'The weighted factor of the initial frequency is out of the range of [0, 1)!');
end

nTotal = nAcc + nConst + nDec;

% initial frequency ------------------------------------------------------%
% range: [2, initFreqUpperLimit)

initFreqUpperLimit = nTotal * s / targTime; % greater than the initial frequency

if initFreqUpperLimit < 3
    error('timeperstep:initFreqUpperLimit', 'The upper limit of the initial frequency should be no less than 3.')
else
    initFreq = fix(2 + wi*(initFreqUpperLimit-2)); % wi is the weighted factor in the range of [0, 1)
end

% maximum frequency (iteration!) -----------------------------------------%

% range of total time: (T1, T1+T2)
% T1 = 2*s*(nAcc+nDec)/(initFreq+maxFreq) + nConst*s/maxFreq;
% T2 = 2*s*(1/initFreq + 1/maxFreq - 4/(initFreq+maxFreq));

% initializing the iteration variables
A = targTime;
B = targTime*initFreq - (nAcc+nDec+nTotal)*s;
C = -nConst*initFreq*s;
% the lower limit of the maximum frequency
iMaxFreq = fix(2*A \ (sqrt(B^2 - 4*A*C) - B));

% the total time will decrease with the increment of the maximum frequency
iTotalTime = inf;

%#TODO If the upper limit of the maximum frequency was introduced, a more 
% efficient iteration algorithm would be developed.
% !IMPORTANT we are to make sure that the total time is no greater than the
% target time.
while iTotalTime > targTime
    % +1 per iteration
    iMaxFreq = iMaxFreq + 1;
    
    % acceleration and deceleration
    A1 = 0.125*s^2*(1/initFreq^2 + 1/iMaxFreq^2);
    B1 = (nAcc - 1)*s;
    B2 = (nDec - 1)*s;
    C1 = 0.5*(initFreq^2 - iMaxFreq^2);
    
    iAcc = 2*A1 \ (sqrt(B1^2 - 4*A1*C1) - B1);
    iDec = -2*A1 \ (sqrt(B2^2 - 4*A1*C1) - B2);
    
    % refresh total time
    iTotalTime = (iMaxFreq - initFreq)*(1/iAcc - 1/iDec) + s*(1/initFreq + (nConst+1)/iMaxFreq);
end

totalTime = iTotalTime;
maxFreq = iMaxFreq;
acc1 = iAcc;
dec1 = iDec;

% time difference is less than or equals to 0
timeDiff = totalTime - targTime;

% frequencies and time periods (DATA TABLE) ------------------------------%

% preallocation for frequencies
nF = nAcc+nDec-1;
F = zeros(nF, 1);
F([1, nAcc, end]) = [initFreq, maxFreq, initFreq];

% acceleration zone
accV = zeros(nAcc+1, 1);

accV([1 end]) = [initFreq - 0.5*acc1*s/initFreq, maxFreq + 0.5*acc1*s/maxFreq];
accV(2:end-1) = sqrt(accV(end)^2 - 2*acc1*s*((nAcc-1):-1:1)');

F(2:nAcc-1) = 0.5 * (accV(2:end-2) + accV(3:end-1));

% deceleration zone
decV = zeros(nDec+1, 1);

decV([1 end]) = [maxFreq - 0.5*dec1*s/maxFreq, initFreq + 0.5*dec1*s/initFreq];
decV(2:end-1) = sqrt(decV(1)^2 + 2*dec1*s*(1:nDec-1)');

F(nAcc+1 : end-1) = 0.5 * (decV(2:end-2) + decV(3:end-1));

F = round(F);

% time periods
T = F .\ s;
T(nAcc) = (nConst+2)*s / maxFreq;

