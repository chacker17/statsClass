function corrMatrix = computeRSM(monkeyInitials, count_window, plotOrNot)
% Purpose: Produces an RSA matrix for images sorted by memorability.

% Inputs:
    % monkeyInitials: String with initials of monkey to pull procData from
    % count_window: 1 x 2 vector with start and end val of window for spike counts
    % plotOrNot: 1 (default) = plots, 0 = doesn't plot
    
% Outputs:
    % corrMatrix contains the n x n correlation matrix with correlations
        % between spike counts for all images and all images
        
% Note: Currrently uses average spike count between novel and familiar
% presentation as spike count for each image.

% Data: 2019 eLife pseudopopulation

% Written 3.20.2020 by CMH

%% Load data, set variables, etc.
if nargin == 2 % Default is to plot
    plotOrNot = 1;
end

if ~exist(['data/' monkeyInitials '_procData_' num2str(count_window(1)) '_' num2str(count_window(2)) '.mat']) % Count spikes manually if data file doesn't already exist
    load(['data/' monkeyInitials '_procData_180_260.mat']);
    spikeTimes.N = procData.novelST;
    spikeTimes.F = procData.familiarST;
    [nUnits, nImages] = size(spikeTimes.N);
    
    for a = 1:nImages % Compute a spike count in that window for each image
        currSTs.N = spikeTimes.N(:, a);
        currSTs.F = spikeTimes.F(:, a);
        
        currImg = [];
        for c = 1:nUnits % Count across the units pooled together in the pseudopopulation
            ST_sess.N = currSTs.N{c};
            ST_sess.F = currSTs.F{c};
            
            spikes.N = sum(ST_sess.N >= (count_window(1) * 10^3) & ST_sess.N < (count_window(2) * 10^3));
            spikes.F = sum(ST_sess.F >= (count_window(1) * 10^3) & ST_sess.F < (count_window(2) * 10^3));
            currSession = (spikes.N + spikes.F)/2;
            
            currImg = [currImg; currSession];
        end
        spikeCounts(:,a) = currImg;
    end
    
else
    load(['data/' monkeyInitials '_procData_' num2str(count_window(1)) '_' num2str(count_window(2)) '.mat']);
    [nUnits, nImages] = size(spikeTimes.N);
    novelSC = procData.novelSC;
    famSC = procData.familiarSC;
    spikeCounts = (novelSC + famSC) / 2; % Averaging spike count across nov/fam presentation
end

mems = procData.memorability; % Memorability for each image as avg of imgs pooled in pseudopopulation
mems = mean(mems, 1);

%% Z-Score by neuron
meanPerUnit = mean(spikeCounts, 2); % Subract for each unit from the response to each unit
meanSub = spikeCounts - meanPerUnit;
stdPerUnit = std(spikeCounts, 0, 2); % Divide by the standard deviation for that unit
zCounts = meanSub ./ stdPerUnit;

%% Arrange by memorability
[mems_sort, sortIdx] = sort(mems); % Arrange by memorability (low to high)
zSorted = zCounts(:, sortIdx);

%% Compute correlations
corrMatrix = [];
for a = 1:nImages
    im1 = zSorted(:, a);
    for b = 1:nImages
        im2 = zSorted(:, b);
        corrMatrix(a, b) = corr(im1, im2);
    end
end

%% Plot
if plotOrNot == 1
    figure(1)
    hold on

    imagesc(corrMatrix)
    colormap(jet)
    axis('square');
    colorbar
    ylabel('Image number');
    xlabel('Image number');
    xlim([0 107]);
    ylim([0 107]);
    set(gca, 'FontSize', 18, 'LineWidth', 2);

    hold off
end
end
