function [stepNums, timeSeqs] = stepsamp(wa, wd, wf, w1, w2, crt, strokeTime)
%STEPSAMP time per STEP SAMPling

% by wulx, 2014/7/11, 2014/12/5


% variable shortcuts
scaleDivs = crt.scaleDivs;
locs = crt.locs;
types = crt.types;
steps = crt.steps;

nRamps = numel(crt.steps) - 1;

% steps difference should not be less than 4
stepDiffs = steps(2:end) - steps(1:end-1);
narrows = stepDiffs<4 | -stepDiffs<4;
broadenMask = [narrows(1); narrows(1:end)] | [narrows(1:end); narrows(end)];

% broaden differences of neighbors
amp1 = 10 * w1 * (types>0); % pull up peaks
amp2 = -10 * w2 * (types<0); % pull down valeys

amp3 = 2 * broadenMask .* types; % broaden the narrow gaps

% % provided that the number of critical points is no less than 3
% if nc > 3
%     if broadenMask(1)
%         amp3(2) = amp3(2) + 2*types(2);
%     end
%     if broadenMask(end)
%         amp3(end-1) = amp3(end-1) + 2*types(end-1);
%     end
% else % n == 3
%     if broadenMask(1)
%         amp3(2) = amp3(2) + 2*types(2);
%     end
% end

% whos broadenMask types amp1 amp2 amp3 steps

steps = round(steps + amp1 + amp2 + amp3);

% pre-allocation for steps and time sequencies
stepNumsInCell = cell(nRamps, 1);
timeSeqsInCell = cell(nRamps, 1);
% timeDiffs = zeros(nc-1, 1);

timeDiff = 0;
for j = 1:nRamps
    nTotal = steps(j+1) - steps(j);
    
    % go backward (-1) or forward (1)
    dirSign = 1;
    if nTotal < 0
        dirSign = -1;
        nTotal = dirSign * nTotal;
    end
    
    % distance
    targDistance = scaleDivs(locs(j+1)) - scaleDivs(locs(j));
    strkLength = scaleDivs(end) - scaleDivs(1);
    targTime = strokeTime*targDistance/strkLength - timeDiff;
    
    if nTotal > 3
        % numbers of steps, [nAcc nConst nDec]
        nAcc = round(2 + wa(j)*(nTotal-4));
        nDec = round(2 + wd(j)*(nTotal-nAcc-2));
        nConst = nTotal - nAcc - nDec;
        % nTotal = nAcc + nDec + nConst;
        
        % time per step algorithm
        [F, T, timeDiff] = timeperstepd(nAcc, nConst, nDec, targTime, wf(j));
        
        % time periods and stepper numbers
        timeSeqsInCell{j} = [T(1:nAcc-1); ones(nConst+2,1)/F(nAcc); T(nAcc+1:end)];
        % stepNumsInCell{j} = (1:nTotal)';
        stepNumsInCell{j} = steps(j) + dirSign*(1:nTotal)';
    elseif nTotal > 0 % nTotal = 1, 2, 3
        disp('total number: 1, 2, 3')

        F = round(nTotal/targTime);
%         while F < 2
%             steps(j+1) = steps(j+1) + dirSign;
%             
%             nTotal = nTotal + 1;
%             F = round(nTotal/targTime);
%         end
                
        T = nTotal/F;
        timeDiff = T - targTime;
        
        timeSeqsInCell{j} = 1/F * ones(nTotal,1);
        stepNumsInCell{j} = steps(j) + dirSign*(1:nTotal)';
    end
    

end

% time compensation ==> sum(timeSeqs) == strokeTime
timeSeqs = [cell2mat(timeSeqsInCell); timeDiff]; 
stepNums = cell2mat(stepNumsInCell);
stepNums = [stepNums; stepNums(end)]; % add the tail

