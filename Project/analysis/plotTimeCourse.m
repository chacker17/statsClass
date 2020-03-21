% Purpose: Plots the timecourse of firing rate averaged across
% images. Also plots this timecourse separately for highest, middle and
% lowest third of image memorabilities.

% Uses a sliding window approach by counting spikes in each window and then
% averaging number of spikes and calculating a firing rate at that
% timepoint.

% Window width: 50 ms
% Step size: 10 ms

% Written 3.19.2020 by CMH

clear all
load('BOTH_procData_180_260.mat');

spikeTimes.N = procData.novelST;
spikeTimes.F = procData.familiarST;
mems = procData.memorability;
mems = mean(mems, 1);
[units, images] = size(spikeTimes.N);

total_time = [-100 400] .* 10^3;
window_width = 50 * 10^3;
window = [total_time(1) - (window_width/2), total_time(1) + (window_width/2)];
step_size = 10 * 10^3;
startPt = total_time(1) - (window_width/2);
endPt = total_time(2) - (window_width/2);
tps = startPt:step_size:endPt; % time points 
numPoints = length(tps);

overallTimecourse = [];
for a = 1:numPoints
    
    curr_window = window;
    
    for b = 1:images %  Calculate the timecourse of FR for each image
        
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
        
        overallTimecourse(b, a) = (sum(currImg)/length(currImg)) / ((curr_window(2) - curr_window(1)) * 10^-6); % FR for response magnitude

    end
    
    window = window + step_size;
    
end

%% Plot average timecourse
avgTimeCourse_all = mean(overallTimecourse, 1);
SEM_timeCourse_all = std(overallTimecourse, 0, 1) ./ sqrt(images);
t = tps ./ 1000;

le = avgTimeCourse_all - SEM_timeCourse_all;
ue = avgTimeCourse_all + SEM_timeCourse_all;

yP=[le,fliplr(ue)];
xP=[t,fliplr(t)];

% Plot
figure(1)
hold on

H.patch=patch(xP,yP,1); % adapted from shadedErrorBar for this plot
set(H.patch,'facecolor','k', ...
    'edgecolor','none', ...
    'facealpha',0.2, ...
    'HandleVisibility', 'off', ...
    'Tag', 'shadedErrorBar_patch')
plot(t, avgTimeCourse_all, 'LineWidth', 2, 'color', 'k');
% plot(t, ue, 'LineWidth', 2, 'color', 'k', 'LineStyle', '--');
% plot(t, le, 'LineWidth', 2, 'color', 'k', 'LineStyle', '--');
xlabel('Time (ms)');
ylabel('Firing rate (spks/s)');
set(gca, 'FontSize', 18, 'LineWidth', 2);

hold off

%% Plot timcourse for low, medium and high memorability images
% Sort into high, medium and low memorability and plot timecourse 
clear le; clear ue
clear xP; clear yP;
[mems_sort, sortIdx] = sort(mems);
mems_low = mems_sort(2:36);
firing_low = mean(overallTimecourse(sortIdx(2:36), :), 1);
SEM_low = std(overallTimecourse(sortIdx(2:36), :), 0, 1) ./ sqrt(length(mems_low));
le.low = firing_low - SEM_low;
ue.low = firing_low + SEM_low;
yP.low = [le.low,fliplr(ue.low)];
xP.low = [t,fliplr(t)];

mems_med = mems_sort(37:71);
firing_med = mean(overallTimecourse(sortIdx(37:71), :), 1);
SEM_med = std(overallTimecourse(sortIdx(37:71), :), 0, 1) ./ sqrt(length(mems_med));
le.med = firing_med - SEM_med;
ue.med = firing_med + SEM_med;
yP.med = [le.med,fliplr(ue.med)];
xP.med = [t,fliplr(t)];

mems_high = mems_sort(72:106);
firing_high = mean(overallTimecourse(sortIdx(72:106), :), 1);
SEM_high = std(overallTimecourse(sortIdx(72:106), :), 0, 1) ./ sqrt(length(mems_high));
le.high = firing_high - SEM_high;
ue.high = firing_high + SEM_high;
yP.high = [le.high,fliplr(ue.high)];
xP.high = [t,fliplr(t)];

% Plot
figure(2)
hold on

plot(t, firing_low, 'LineWidth', 2, 'color', 'k');
plot(t, firing_med, 'LineWidth', 2, 'color', 'b');
plot(t, firing_high, 'LineWidth', 2, 'color', 'r');
H.patch=patch(xP.low,yP.low,1); % adapted from shadedErrorBar for this plot
set(H.patch,'facecolor','k', ...
    'edgecolor','none', ...
    'facealpha',0.2, ...
    'HandleVisibility', 'off', ...
    'Tag', 'shadedErrorBar_patch')
H.patch=patch(xP.med,yP.med,1); % adapted from shadedErrorBar for this plot
set(H.patch,'facecolor','b', ...
    'edgecolor','none', ...
    'facealpha',0.2, ...
    'HandleVisibility', 'off', ...
    'Tag', 'shadedErrorBar_patch')
H.patch=patch(xP.high,yP.high,1); % adapted from shadedErrorBar for this plot
set(H.patch,'facecolor','r', ...
    'edgecolor','none', ...
    'facealpha',0.2, ...
    'HandleVisibility', 'off', ...
    'Tag', 'shadedErrorBar_patch')
xlabel('Time (ms)');
ylabel('Firing rate (spks/s)');
set(gca, 'FontSize', 18, 'LineWidth', 2);
legend({'low: 0.44', 'med: 0.67', 'high: 0.84'}, 'FontSize', 18);

hold off
    