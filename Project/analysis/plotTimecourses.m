% Purpose: Computes and plots timecourse of memorability and identity
% information in IT. Memorability is measured as the magnitude of the
% correlation between population response magnitude and image memorability
% and identity is measured as the performance of a 107-way classifier.

% Written 4.13.2020 by CMH

%% Compute vectors
 total_time = [-100 400] .* 10^3;
window_width = 50 * 10^3;
step_size = 5 * 10^3;

[timePts1, memInfo] = corrByTime(window_width, step_size, total_time); % Custom functions
[timePts2, idInfo] = imageIDbyTime(window_width, step_size, total_time);

%% Plot
figure('defaultAxesColorOrder', [[0 0 0]; [1 0 0]]);

hold on

yyaxis left % Identity -- black
plot(timePts2 * 10^-3, idInfo, 'k-', 'LineWidth', 2);
ylabel('Classifier Performance'); 

yyaxis right % Memorability -- red
plot(timePts1, memInfo, 'r-', 'LineWidth', 2);
ylabel('Pearson Corr'); % Should I try this with Spearman later?

hold off

xlabel('Time (ms)');
set(gca, 'FontSize', 18, 'LineWidth', 2);
box off

