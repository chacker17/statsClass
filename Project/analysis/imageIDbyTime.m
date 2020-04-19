function [timePts, idAccByTime, errorByTime] = imageIDbyTime(window_width, step_size, total_time)

% Purpose: Plots the accuracy of a 107-way classifier using a sliding
% window analysis.

% Inputs: (all should be given in microseconds)
    % window_width: Width of the sliding window 
    % step_size: Size of the step for sliding window analysis
    % total_time: Total window over which to run the analysis
    % Defaults:
        % window_width: 50 ms
        % step_size: 5 ms
        % total_time: [-100 400]

% Outputs:
    % timePts: Row vecotr of pts with the center of each window
    % idAccByTime: Row vector of pts with value of accuracy of 107-way
        % classifier 

% Written 4.13.2020 by CMH

%% Load data, set variables, etc.
if nargin ~= 3
    total_time = [-100 400] .* 10^3;
    window_width = 50 * 10^3;
    step_size = 5 * 10^3;
end

load('BOTH_procData_180_260.mat');

spikeTimes.N = procData.novelST; % Pull spike info
spikeTimes.F = procData.familiarST;
mems = procData.memorability; % Image memorabilities
mems = mean(mems, 1);
[units, images] = size(spikeTimes.N);

window = [total_time(1) - (window_width/2), total_time(1) + (window_width/2)]; % Generate window info
startPt = total_time(1) - (window_width/2);
endPt = total_time(2) - (window_width/2);
timePts = startPt:step_size:endPt; % time points 
numPoints = length(timePts);

%% Classifier performance over time
for a = 1:numPoints % Slide window over time
    a
    
    curr_window = window;
    
    for b = 1:images
        
        for c = 1:units
            
            currSTs.N = spikeTimes.N{c, b};
            currSTs.F = spikeTimes.F{c, b};
            
            spikes.N = sum(currSTs.N >= curr_window(1) & currSTs.N < curr_window(2));
            spikes.F = sum(currSTs.F >= curr_window(1) & currSTs.F < curr_window(2));
            
            SC.N(c, b) = spikes.N;
            SC.F(c, b) = spikes.F;
            
        end
    end
    
    [perf1, error1] = classifierPerformance(SC.N, SC.F); % Train fam, test nov
    [perf2, error2] = classifierPerformance(SC.F, SC.N); % Train nov, test fam
    
    idAccByTime(a) = (perf1 + perf2) / 2; % Average training with fam or nov
    errorByTime(a) = (error1 + error2) / 2;
    
    window = window + step_size;
    
end

timePts = timePts * 10^-3; % Put into ms for output

end

%% Helper functions
function [perf, err] = classifierPerformance(spikes_train, spikes_test)

[units, images] = size(spikes_train);
angleDist = nan(images, images);
numReplicates = 30;

for a = 1:images % Training set
    
    for b = 1:images % With all images in test set
        
        trainVec = spikes_train(:, a);
        testVec = spikes_test(:, b);
        
        angleDist(a, b, 1) = abs(corr(trainVec(:), testVec(:), 'Type', 'Pearson'));
        SEM = sqrt((1-angleDist(a, b)^2)/(images-2)); % Get standard error of the correlation
        
        for c = 2:numReplicates % Randomly choose a correlation from the distribution and save 100 times (to be used to generate error terms
            
            angleDist(a, b, c) = normrnd(angleDist(a, b, 1), SEM);
            
        end
        
    end
    
end

perf_temp = [];
for d = 1:numReplicates

    max_dist = max(angleDist(:, :, d), [], 2);
    testMax = angleDist(:, :, d) == max_dist;
    actMax = eye(images, images);

    correct_ind = testMax .* actMax;
    numCorrect = sum(correct_ind, 2);

    perf_temp(d) = sum(numCorrect) / length(numCorrect); % Actual performance
end

perf = mean(perf_temp);
err = std(perf_temp) / sqrt(numReplicates);

end