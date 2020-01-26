% Purpose: Code to answer the exercises in quantal release example for
% binomial distribution reading for quant neuro course.

% Link to exercise page: https://canvas.upenn.edu/courses/1358934/discussion_topics/5002068

% Written 1.24.20 by CMH

%% Exercise 1
n = 10; % 10 quanta
p = 0.2; % Each with release probability of 0.2

theoreticalProbs = binopdf(0:n, n, p);
hold on
plot(0:10,theoreticalProbs,'ro-', 'LineWidth', 2, 'MarkerSize', 10);
title(sprintf('Predicted number of %d quanta released, each with PoR of %.2f', n, p));
xlabel('Number of quanta released');
ylabel('Probability');
hold off

%% Exercise 2
n = 14; % 14 quanta available for release
p = 0.1; % Release probability of 0.1 for each
theoreticalProbs = binopdf(0:n, n, p);

clf;
figure(1)
hold on
plot(0:14,theoreticalProbs,'ro-', 'LineWidth', 2, 'MarkerSize', 10);
title(sprintf('Predicted number of %d quanta released, each with PoR of %.2f', n, p));
xlabel('Number of quanta released');
ylabel('Probability');
hold off

% Probability of getting a value of 8 with given parameters
fprintf('Probability of getting a value of 8 with release probability of 0.1 is %.5f\n', theoreticalProbs(9));

% Probability of getting a value of 8 with increasing release probability
p = 0;
probsByProb = zeros(1);
probsCounted = 1;
for a = 0:0.1:1
    tempProbs = binopdf(0:n, n, a);
    fprintf('Probability of getting a value of 8 with release probability of %.2f is %.3f percent\n', a, tempProbs(9) * 100);
    probsByProb(probsCounted) = tempProbs(9) * 100;
    probsCounted = probsCounted + 1;
end

% Plot probability of release of 8 quanta for each release probability
figure(2)
hold on
bar(0:0.1:1, probsByProb);
title('Probability of releasing 8/14 quanta across release probabilities');
xlabel('Probability of release for each quanta');
ylabel('Probability that 8/14 quanta are released');
hold off

%% Exercise 3
n = 14; % 14 quanta available for release
p = 0.1; % Release probability of 0.1 for each
theoreticalProbs = binopdf(0:n, n, p);
prob8 = theoreticalProbs(9);
prob5 = theoreticalProbs(6);
totalLik = prob8 * prob5;

fprintf('Total likelihood of these events is %.8f\n', totalLik);

logLik = log(prob8) + log(prob5);

fprintf('Log-likelihood of these events is %.0f\n', logLik);

% Repeat for all deciles of PoR
totLiks = zeros(1);
logLiks = zeros(1);
releaseProbs = 0:0.05:1;
cnt = 1;
for a = releaseProbs
    tempProb = binopdf(0:n, n, a);
    totLik = prob8 + tempProb(6);
    totLiks(cnt) = totLik;
    logLik = log(prob8) + log(tempProb(6));
    logLiks(cnt) = logLik;
    
    fprintf('Total likelihood with %.2f PoR is %.5f\n', a, totLik);
    fprintf('Log likelihood with %.2f PoR is %.5f\n', a, logLik);
    
    cnt = cnt + 1;
end

clf;
figure(1)
hold on

plot(releaseProbs, totLiks, 'ro-', 'LineWidth', 2, 'MarkerSize', 10);
xlabel('Release probability');
ylabel('Total likelihood of 8 then 5 released');
title('Total likelihood of 8 then 5 quanta released for each release prob');

hold off

figure(2)
hold on

plot(releaseProbs, logLiks, 'ro-', 'LineWidth', 2, 'MarkerSize', 10);
xlabel('Release probability');
ylabel('Log likelihood of 8 then 5 released');
title('Log likelihood of 8 then 5 quanta released for each release prob');

hold off

% See what happens as you increase the sample size
% Animation that updates the figure in a loop as you 

totNumSamples = 50;
n = 14;
p = 0.3;
for a = 1:3:51

    quantaReleased = binornd(n, p, a, 1);
    
    % Find likelihood of these quanta being released at each PoR
    releaseProbs = 0:0.05:1;
    totLiks = zeros(1, length(releaseProbs));
    logLiks = zeros(1, length(releaseProbs));
    allProbs = zeros(1, a);
    cnt = 1;
    for c = releaseProbs  
        tempProb = binopdf(0:n, n, c);
        
        % Loop through all samples and find probabilities
        for d = 1:length(quantaReleased)
            currQuanta = quantaReleased(d);
            allProbs(d) = tempProb(currQuanta + 1); % Offset index by 1 bc indexing in matlab starts at 1
        end

        totLiks(cnt) = sum(allProbs);
        logLiks(cnt) = sum(log(allProbs));
        cnt = cnt + 1;
        
    end
    
    
    figure(3);
    
    subplot(1, 2, 1) % Total likelihood 
    cla reset;
    hold on
    
    plot(releaseProbs, totLiks, 'ro-', 'LineWidth', 2, 'MarkerSize', 10);
    xlabel('Release probability');
    ylabel('Total likelihood');
    title(sprintf('Total likelihood of random releases with each PoR for %d samples', a));
    
    subplot(1, 2, 2)
    cla reset;
    hold on

    plot(releaseProbs, logLiks, 'ro-', 'LineWidth', 2, 'MarkerSize', 10);
    xlabel('Release probability');
    ylabel('Log likelihood');
    title(sprintf('Log likelihood of random releases with each PoR for %d samples', a));
    
    
    pause(0.4);
end


%% Exercise 4


%% Exercise 5

