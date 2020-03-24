% Purpose: Plots the distribution of memorabilities for the 107 images in
% the pseudopopulation. 

% Memorability for an image is the average of the memorabilities of the
% images that were pooled together to make the pseudopopulation. In this
% case, 27 images.

% Written 3.23.2020 by CMH

%% Load data, set variables, etc.
load('BOTH_procData_180_260.mat');

mems = procData.memorability;
mems = mean(mems, 1);

%% Plot
figure(1)
hold on

histogram(mems);
set(gca, 'FontSize', 18, 'LineWidth', 2);
xlabel('Memorability');
ylabel('Number of images');

hold off