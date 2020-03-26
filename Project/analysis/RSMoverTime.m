% function RSMoverTime
% Purpose: Computes a RSM using a sliding window analysis and plots the
% average correlation in that RSM over time.

% Written 3.26.2020 by CMH

%% Load data, set variables, etc.
total_time = [-100 400];
window_width = 50;
window = [total_time(1) - (window_width/2), total_time(1) + (window_width/2)];
step_size = 10;
startPt = total_time(1) - (window_width/2);
endPt = total_time(2) - (window_width/2);
tps = startPt:step_size:endPt; % time points 
numPoints = length(tps);

%% Compute RSM in each window -- save in 3d matrix and compute avg corr of RSM in another array
for a = 1:numPoints
    curr_window = window;
    fprintf('Computing RSM for %.0f to %.0f ms...\n', curr_window(1), curr_window(2));
    
    RSM = computeRSM('BOTH', curr_window, 0); % Compute RSM but don't plot
	RSM(eye(length(RSM)) == 1) = NaN; % Don't include identity line with correlations of 1
    
    RSMbyTime(:, :, a) = RSM;
    avgCorrs(a) = nanmean(RSM(:));
    
    corrsUnique = tril(RSM);
    corrsUnique(isnan(NaN)) = 0;
    numUnique = length(find(corrsUnique ~= 0));
    semCorrs(a) = nanstd(RSM(:)) / sqrt(numUnique);
    
    window = window + step_size;
end

%% Plot average corr from RSM over time
le = avgCorrs - semCorrs;
ue = avgCorrs + semCorrs;

yP=[le,fliplr(ue)];
xP=[tps,fliplr(tps)];
figure(1)
hold on

plot(tps, avgCorrs, 'k.-', 'LineWidth', 2);
H.patch=patch(xP,yP,1); % adapted from shadedErrorBar for this plot
set(H.patch,'facecolor','k', ...
    'edgecolor','none', ...
    'facealpha',0.2, ...
    'HandleVisibility', 'off', ...
    'Tag', 'shadedErrorBar_patch')
set(gca, 'LineWidth', 2, 'FontSize', 18);
xlabel('Time (ms)');
ylabel('Average correlation in RSM');
ylim([-0.015 0]);

hold off

% end