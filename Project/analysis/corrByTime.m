function [timePts, corrOverTime, errorOverTime] = corrByTime(window_width, step_size, total_time)
% Purpose: Computes the correlation between response magnitude and image
% memorability over time using a sliding window analysis.

% Inputs: (all should be given in microseconds)
    % total_time: Total window over which to run the analysis
    % window_width: Width of the sliding window 
    % step_size: Size of the step for sliding window analysis
    % Defaults:
        % window_width: 50 ms
        % step_size: 5 ms
        % total_time: [-100 400]

% Outputs:
    % timePts: Row vecotr of pts with the center of each window
    % corrOverTime: Row vector of pts with value of correlation in each
        % window centered with the times given in timePts

% Written 4.13.2020 by CMH
% Updated 4.16.2020 to add error bars by CMH

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

%% Compute correlation over time
corrOverTime = [];
for a = 1:numPoints % Slide the window across all points
    
    curr_window = window;
    RMperIm = [];
    
    for b = 1:images % Calculate the correlation for all images in that window
        
        currSTs.N = spikeTimes.N(:, b);
        currSTs.F = spikeTimes.F(:, b);
        
        currImg = [];
        for c = 1:units % Count across the units pooled together in the pseudopopulation
            currSession = [];
            ST_sess.N = currSTs.N{c};
            ST_sess.F = currSTs.F{c};
            
            spikes.N = sum(ST_sess.N >= curr_window(1) & ST_sess.N < curr_window(2));
            spikes.F = sum(ST_sess.F >= curr_window(1) & ST_sess.F < curr_window(2));
            currSession = (spikes.N + spikes.F)/2;
            
            currImg = [currImg; currSession];
        end
        
        currRM = sqrt(sum(currImg .^2, 1));
        RMperIm = [RMperIm, currRM];
        
    end
    
    currCorr = corrcoef(mems, RMperIm);
    corrOverTime(a) = currCorr(1, 2); % Correlation
    errorOverTime(a) = sqrt((1-corrOverTime(a)^2)/(images-2)); % SEM of the correlation
    
    window = window + step_size;
    
end

timePts = timePts * 10^-3; % Put back into ms to output