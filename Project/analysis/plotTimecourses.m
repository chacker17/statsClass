function plotTimecourses(defaultVals)
% Purpose: Computes and plots timecourse of memorability and identity
% information in IT. Memorability is measured as the magnitude of the
% correlation between population response magnitude and image memorability
% and identity is measured as the performance of a 107-way classifier.

% Inputs
    % defaultVals: String with the name of the variables that contain the
        % values to be plotted (saves a considerable amount of 
        % computational time)

% Written 4.13.2020 by CMH

%% Compute vectors
total_time = [-100 400] .* 10^3;
window_width = 50 * 10^3;
step_size = 5 * 10^3;

if nargin < 1
    [timePts1, memInfo, memError] = corrByTime(window_width, step_size, total_time); % Custom functions
    [timePts2, idInfo, idError] = imageIDbyTime(window_width, step_size, total_time);
else
    load(defaultVals);
end

%% Plot
le_mem = memInfo - memError; % Compute error bar values for memorability info
ue_mem = memInfo + memError;
yP_mem = [le_mem, fliplr(ue_mem)];
xP_mem =[timePts1, fliplr(timePts1)];

le_id = idInfo - idError; % Compute error bar values for identity info
ue_id = idInfo + idError;
yP_id = [le_id, fliplr(ue_id)];
xP_id = [timePts2, fliplr(timePts2)];

figure('defaultAxesColorOrder', [[0 0 0]; [1 0 0]]);

hold on

yyaxis left % Identity -- black
H.patch=patch(xP_id,yP_id,1); % adapted from shadedErrorBar for this plot
set(H.patch,'facecolor','k', ...
    'edgecolor','k', ...
    'facealpha',0.2, ...
    'HandleVisibility', 'off', ...
    'Tag', 'shadedErrorBar_patch')
plot(timePts2, idInfo, 'k-', 'LineWidth', 2);
ylabel('Classifier Performance'); 

yyaxis right % Memorability -- red
H.patch=patch(xP_mem,yP_mem,1); % adapted from shadedErrorBar for this plot
set(H.patch,'facecolor','r', ...
    'edgecolor','r', ...
    'facealpha',0.2, ...
    'HandleVisibility', 'off', ...
    'Tag', 'shadedErrorBar_patch')
plot(timePts1, memInfo, 'r-', 'LineWidth', 2);
ylabel('Pearson Corr'); % Should I try this with Spearman later?

hold off

xlabel('Time (ms)');
set(gca, 'FontSize', 18, 'LineWidth', 2);
box off

end
