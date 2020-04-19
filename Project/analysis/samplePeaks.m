% Purpose: Randomly samples 1000 new traces for memorability and identity
% information using the means and standard deviations given as the output
% of corrByTime and imageIDbyTime. For each trace, save the time at which
% the peak occurs and then run a t-test to determine whether the
% distribution of peak times for memorability and identity information is
% significantly different. 

% Written 4.17.2020 by CMH

%% Load data, set variables, etc.
load('defaultVals_error.mat'); % Load data points
numSimulations = 20;

%% Loop through 
idPeaks = []; memPeaks = [];
for a = 1:numSimulations % Simulate several times
    a
    
    for b = 1:86 % Look at first 300 ms 
        % NOTE: Chose to look at the first 300 ms instead of all the info I
        % had data for, because memorability trace has sustained response
        % and that turns up some very late 
        
        memVals(b) = normrnd(memInfo(b), memError(b));
        idVals(b) = normrnd(idInfo(b), idError(b));
        
    end
    
    % Find peak for mem and id and save to separate variables
    idPeaks(a) = timePts2(idVals == max(idVals));
    memPeaks(a) = timePts1(memVals == max(memVals));
    
end

%% Stats -- Mann-Whitney U test
% Non-parametric test to see if distributions are different
combinedPeaks = [idPeaks memPeaks];
[sortedPeaks, idx] = sort(combinedPeaks);
distNum = [ones(length(idPeaks), 1); ones(length(memPeaks), 1) + 1]; % 1 = identity, 2 = memorability
numID = length(idPeaks);
numMem = length(memPeaks);

% Assign ranks where tied values are grouped
possibleTimes = unique(combinedPeaks);

ranks = zeros(1, length(combinedPeaks));
for a = 1:length(possibleTimes)
    
    currIdx = find(sortedPeaks == possibleTimes(a));
    assignedRank = mean(currIdx);
    ranks(currIdx) = assignedRank;
    
end

sumID = sum(ranks(distNum == 1));
U1 = sumID - ((numID * (numID + 1)) / 2);
sumMem = sum(ranks(distNum == 2));
U2 = sumMem - ((numMem * (numMem + 1)) / 2);
mU = (numID * numMem) / 2;
% sigmaU = sqrt((numID * numMem * (numID + numMem + 1))/12);
U = min([U1 U2]); % Use smaller value of U for significance tests

possibleRanks = unique(ranks);
for a = 1:length(possibleRanks) % Correct for having lots of ties
    
    numRank = possibleRanks(a);
    currNum = length(find(ranks == numRank));
    vals(a) = (currNum^3 - currNum)/((numID + numMem) * ((numID + numMem) - 1));
    
end

sumCorrection = sum(vals);

sigmaU = sqrt((numID * numMem /12) * ((numID + numMem + 1) - sumCorrection));
z = (U - mU)/sigmaU;
p = normcdf(z);
f = U/(numID * numMem);
fprintf('z = %.3f, p = %.7f, f = %.3f\n', z, p, f);

%% Plot
[dummy, centerBoth] = hist([idPeaks; memPeaks], 80);
[idCounts, idCenters] = hist(idPeaks, centerBoth);
idProb = transpose(idCounts/sum(idCounts));
[memCounts, memCenters] = hist(memPeaks, centerBoth);
memProb = transpose(memCounts/sum(memCounts));

figure(1)
hold on

% bar(idCenters, idProb);
% bar(memCenters, memProb);
b = bar(centerBoth, [idProb, memProb], 5);
b(1).FaceColor = [0 0 0];
b(2).FaceColor = [1 0 0];
set(gca, 'LineWidth', 2, 'FontSize', 18);
xlim([75 300]);
xlabel('Time (ms)');
ylabel('Probability');
legend('Identity Peak', 'Memorability Peak');
box off

hold off

