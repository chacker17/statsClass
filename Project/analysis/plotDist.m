% Purpose: Plots the distributions of spike counts across units for each
% image.

% Written 3.20.2020 by CMH

load('BOTH_procData_180_260.mat');

novelSpikes = procData.novelSC;
famSpikes = procData.familiarSC;
spikeCount = (novelSpikes + famSpikes) ./ 2;
mems = procData.memorability;
mems = mean(mems, 1);

%% Average response across all images
avgRespPerIm = mean(spikeCount, 2);

% Histogram of average spike response for 707 units to an image
figure(1)
hold on

histogram(avgRespPerIm);
set(gca, 'FontSize', 18, 'LineWidth', 2);
xlabel('Spike Count');
ylabel('Number of neurons');
box off

hold off

%% High v. med v. low memorability images
[mems_sort, sortIdx] = sort(mems);
mems_low = mems_sort(2:36);
avgResp_low = mean(spikeCount(:, sortIdx(2:36)));

mems_med = mems_sort(37:71);
avgResp_med = mean(spikeCount(:, sortIdx(37:71)));

mems_high = mems_sort(72:106);
avgResp_high = mean(spikeCount(:, sortIdx(72:106)));

figure(2)
hold on

histogram(avgResp_low, 7);
histogram(avgResp_med, 7);
histogram(avgResp_high, 7);

set(gca, 'FontSize', 18, 'LineWidth', 2);
xlabel('Spike Count');
ylabel('Number of neurons');
box off

legend({'low: 0.44', 'med: 0.67', 'high: 0.84'}, 'FontSize', 18);

hold off
