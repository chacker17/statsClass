% Purpose: Coded answers to exercises in the linear regression discusion on
% Canvas. 

% Canvas Discussion: https://canvas.upenn.edu/courses/1358934/discussion_topics/5116781

% Josh Answers: 
    % Ex. 1: https://github.com/PennNGG/Statistics/blob/master/Simple%20Linear%20Regression/LinearRegressionExample.m
    % Ex. 2: https://github.com/PennNGG/Statistics/blob/master/Simple%20Linear%20Regression/TestofLinearity.m
    
% Written 2.18.2020 by CMH

%%%%%%%%% Exercise 1 %%%%%%%%%%%%%%%

age = [3 4 5 6 7 8 9 11 12 14 15 16 17];
age = [3 4 5 6 8 9 10 11 12 14 15 16 17]; % Updated to match age in Josh's code so I can check answers
wingLength = [1.4 1.5 2.2 2.4 3.1 3.2 3.2 3.9 4.1 4.7 4.5 5.2 5];
n = length(age);

%% Q1 Plot 
figure(1) 
hold on

plot(age, wingLength, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
set(gca, 'FontSize', 18, 'LineWidth', 2);
xlabel('Age (years)');
ylabel('Wing Length(cm)');
box off

hold off

%% Q2 Calcualte the regression line 
b = (sum(age .* wingLength) - ((sum(age) * sum(wingLength))/n)) / (sum(age.^2) - ((sum(age) ^ 2)/n));
xbar = mean(age);
ybar = mean(wingLength);
a = ybar - (b*xbar);

figure(1) % Add to the plot
hold on

plot([1 19], b*[1 20] + a, 'k-', 'LineWidth', 2);

hold off

% y = 0.27x + 0.71

%% Q3 Can you reject the null b = 0?
df = n-2;
num = (sum(age .* wingLength) - (sum(age)*sum(wingLength)/n))^2 / ((sum(age.^2) - (sum(age)^2/n)));
denom = ((sum(wingLength.^2) - (sum(wingLength)^2/n)) - num)/df;
F = num/denom;
xunder = 1./max(0,F); % Adapted from Josh's code
xunder(isnan(F)) = NaN;
p = fcdf(xunder,df,1);

% Yes, with a p value << 0.001 you can reject the null that b = 0

%% Q4 What is r?
r_table = corrcoef(age, wingLength);
r = r_table(1, 2);

% r is 0.99

%% Q5 Add some noise and see what happens
% Completed by changing the values up top and running the code several
% times