% Purpose: Correlates firing rate and memorability separately for
% novel and familiar presentations of the images and plots.

% Data: 2019 eLife, Figure 2
% Adapated from VM181113_respVsMemPlot3_forPrint2

% Written 11.15.19 by CMH
% Adapted 3.23.20 by CMH to produce simpler plot for data visualization presentation

%% Set up: Set variables, load data, etc.
clear all
monkeyInitials = 'BOTH'; % Use data from both monkeys
monkeyInitials = upper(monkeyInitials);

timeBin_ms = [180 260];
widthWindow = (timeBin_ms(2) - timeBin_ms(1))/1000; % in seconds

fileName2load = sprintf('%s_procData_%0.0f_%0.0f.mat', monkeyInitials, timeBin_ms(1), timeBin_ms(2));

load(fileName2load)

memorability = mean(procData.memorability, 1); 

NMag = mean(procData.novelSC/widthWindow, 1); % Firing rate
FMag = mean(procData.familiarSC/widthWindow,1);

avgMag = (NMag + FMag) ./ 2; % Average response across novel v. familiar presentation

%% Linear regression
linModel = struct;
x = memorability;
y = avgMag;

fitRes = fitlm(x, y);
CoefficientNames = fitRes.CoefficientNames;
interceptInd = strcmp(CoefficientNames, '(Intercept)');
linModel.intercept.estimate =  fitRes.Coefficients.Estimate(interceptInd);
linModel.slope.estimate = fitRes.Coefficients.Estimate(~interceptInd);
linModel.intercept.pValue =  fitRes.Coefficients.pValue(interceptInd);
linModel.slope.pValue = fitRes.Coefficients.pValue(~interceptInd);
linModel.rSq = fitRes.Rsquared.Ordinary;

m = linModel.slope.estimate;
d = linModel.intercept.estimate;

%% Correlation
rPearson = corr(x(:), y(:), 'type', 'pearson');

%% Plot 
hold on

Pts = plot(memorability, y, 'marker', '.', 'MarkerSize', 15, 'linestyle', 'none', 'color', 'k'); % Ind. images
x = [0.2, 1];
line.N = plot(x, m*x+d, 'marker', 'none', 'linestyle', '-', 'linewidth', 1.5, 'color', 'k'); % Linear regression

% Clean up plot
xlim([0.1, 1.1]);
ylim([10 18]);
set(gca, 'FontSize', 18, 'LineWidth', 2);
xlabel('Image memorability');
ylabel('Firing rate (spks/s)');

hold off