% Purpose: Code written in response to Exercise 1 of the parametric
% correlation tests discussion on Canvas.

% Link: https://canvas.upenn.edu/courses/1358934/discussion_topics/5600885

% Josh's solution: 

% Written 2.16.2020 by CMH

wingLength = [10.4 10.8 11.1 10.2 10.3 10.2 10.7 10.5 10.8 11.2 10.6 11.4]; % X - cm
mean_WL = mean(wingLength);
tailLength = [7.4 7.6 7.9 7.2 7.4 7.1 7.4 7.2 7.8 7.7 7.8 8.3]; % Y - cm
mean_TL = mean(tailLength);
n = length(wingLength);
alpha = 0.05;

%% Plot against each other
figure(1) 
hold on

plot(wingLength, tailLength, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
set(gca, 'FontSize', 18, 'LineWidth', 2);
xlabel('Wing Length (cm)');
ylabel('Tail Length(cm)');
xlim([10 11.6]);
box off

hold off

% 1. Yes, these two variables look like they are related.

%% Calculate r(x,y) and r(y,x) and then with corrcoef
num = sum((wingLength - mean_WL) .* (tailLength - mean_TL));
x_denom = sqrt(sum((wingLength - mean_WL) .^ 2));
y_denom = sqrt(sum((tailLength - mean_TL) .^ 2));

r_xy_hand = num/(x_denom * y_denom);
r_yx_hand = num/(y_denom * x_denom);
r_matlab = corrcoef(wingLength, tailLength);
r_xy_matlab = r_matlab(1, 2);
r_yx_matlab = r_matlab(2, 1);
 
r = r_xy_hand;

% 2. Yes, the value for r I calcualted by hand matches the one in Matlab.

%% Find standard error and 95% CI
SE = sqrt((1 - r^2)/(n - 2));

FZ = 0.5 * log((1 + r)/(1 - r));
sz = sqrt(1/(n-3));

z_lower = FZ + norminv(alpha/2) * sz;
z_upper = FZ - norminv(alpha/2) * sz;

r_upper = (exp(2*z_upper) - 1)/(exp(2*z_upper) + 1);
r_lower = (exp(2*z_lower) - 1)/(exp(2*z_lower) + 1);

% 3. Standard error is 0.1557 and CI is [0.5923 0.9632];

%% Find significance
t = r/SE;
p = 1-tcdf(t, n - 2);

% 4. Yes, p = 0.00015

%% Is r = 0.75 different?
zm = FZ;
rh = 0.75;
zh = 0.5 * log((1 + rh)/(1 - rh));

z = (zm - zh)/(sqrt(1/(n - 3)));

t_new = 1-normcdf(z/2);

% 5. No, in this case p = 0.2938 so we do not reject the null

%% Calculate power and sample size needed to reject r = 0 when r > 0.5
% Pieces of code below adapted from Josh's answer key
df = n-2;
tcrit=tinv(1-alpha/2,df);
rcrit=sqrt(tcrit^2/(tcrit^2+(df)));
zr=0.5*log((1+rcrit)/(1-rcrit));

Zb=(zm-zr)*sqrt(n-3);
power=normcdf(Zb);

% 6. The power is 0.98

percReject = 0.01;
rho = 0.5;
Zb = tinv(1 - percReject,inf);
Za = tinv(1 - alpha/2, inf);
zeta= 0.5 * log((1 + rho)/(1 - rho));
SampleSize = round(((Zb + Za)/zeta)^2+3);

% 6. We need a sample size of 64