% laterTutorial_modelFits
%
% Script for fitting the LATER model to real data
%
% Recall that the idea behind the LATER (Linear Approach to Threshold
%  with Ergodic Rate) model is that a saccadic decision is made when a
%  "ballistic" decision variable rises from a starting value to an ending
%  value at a constant ("ergodic") rate that is fixed for a particular
%  saccadic decision but varies as a Gaussian across decisions reaches a
%  fixed threshold.
%
% Because the rate of rise is distributed as a Gaussian, and the distance
%  traveled by the decision variable is fixed, we know that the amount of
%  time, representing the reaction time (RT), is just:
%
%  RT = distance / rate
%
% And thus RT is distributed as 1/rate, or an inverse Gaussian -- which is
%  the primary empirical observation that the LATER model was developed to
%  explain.
%
% Copyright 2019 by Joshua I. Gold, University of Pennsylvania

%% Set up fitting
%
% We will use the same fitting procedure on a bunch of data sets, so we set
%  the options here. We will be using GlobalSearch. The general
%  advantage of this approach is to avoid local minima; for details, see:
%  https://www.mathworks.com/help/gads/how-globalsearch-and-multistart-work.html
%  
% These options seem to work well, but I don't have a stronger
%  rationale for using them. See the Matlab documentation if you really
%  want to dive in and understand them, and let me know if you find
%  better settings!
opts = optimoptions(@fmincon,    ... % "function minimization with constraints"
   'Algorithm',   'active-set',  ...
   'MaxIter',     3000,          ...
   'MaxFunEvals', 3000);

%% Example session and fits
%
%  We use the same example from the scriptLATER_rawData script,
%  corresponding to Fig. 2, Kim et al, J Neuroscience

% Use a particular subject tag
SUBJECT_TAG = 'JT';

% Load the data from that subject, given the base directory
load(fullfile('LATERdata', 'data_mgl', 'F', ...
   [SUBJECT_TAG '_RT.mat']));

% Select the appropriate data (see scriptLATER_rawData for details)
%  C_L,0:  Left choices, change-point trials
%  C_L,1+: Left choices, non-change-point trials
%  C_R,0:  Right choices, change-point trials
%  C_R,1+: Right choices, non-change-point trials
labels        = {'C_L_,_0', 'C_L_,_1_+', 'C_R_,_0', 'C_R_,_1_+'};
colors        = {'b' 'r' 'y' 'm'};
Ltrials       = percorrSum == 1 & tRxnSum < 1.2;
LtrialSubsets = [ ...
   Ltrials & numdirSum == -1 & labelSum == 1; ...
   Ltrials & numdirSum == -1 & labelSum ~= 1; ...
   Ltrials & numdirSum ==  1 & labelSum == 1; ...
   Ltrials & numdirSum ==  1 & labelSum ~= 1];

% Pick a cutoff for express saccades
expressCutoff = 0.2;

% Save the fits
fits = nan.*zeros(length(labels), 2);

% Open a figure
figure
subplot(2,1,1); hold on;

% Loop through each data set
for ii = 1:length(labels)

   % Get the data, convert to column vector
   RTs = tRxnSum(LtrialSubsets(ii,:))';
   
   % Convert to reciprocal RT and take only non-express saccades for fitting
   rRTs = 1./RTs(RTs>expressCutoff);
   
   % Pick initial values, using empirical mean/std of reciprocal RTs
   deltaS0 = 1/std(rRTs);
   muR0 = mean(rRTs)*deltaS0;
   
   % Define the objective function. This could be set up as Matlab
   %  function, as it is here:
   %  https://github.com/PennNGG/Statistics/blob/master/LATER%20model/reciprobit_err.m
   %
   %  But here we use an Anonymous Function, so you can see it
   %
   %  Arguments:
   %  fits is a 1x2 vector of the two free parameters:  [muR deltaS]
   %     As desribed in scriptLATER_modelParameters, the model assumes that:
   %        1/RT ~N(muR/deltaS, 1/deltaS)
   %  Also uses rRTs, the nx1 vector of per-trial reciprocal rts (in seconds)
   %        defined above
   %
   %  Returns the negative summed log-likelihood:
   %  a. Compute the likelihood as:
   %        p(data | model parameters) =
   %        p( 1/RT | muR, deltaS) =
   %        normpdf( 1/RT, muR/deltaS, 1/deltaS) =
   %        normpdf(rrts, fits(1)/fits(2), 1/fits(2))
   %  b. Take the logarithm of each likelihood
   %  c. Sum them all together
   %  d. Take the negative, because fmincon finds the minimum (thus
   %        corresponding to the maximum likelihood)
   laterErrFcn = @(fits) -sum(log(normpdf(rRTs, fits(1)/fits(2), 1/fits(2))));
   
   % Set up the optimization problem
   problem = createOptimProblem('fmincon',    ...
      'objective',   laterErrFcn,      ... % Use the objective function
      'x0',          [deltaS0 muR0],   ... % Initial conditions
      'lb',          [0.001 0.001],    ... % Parameter lower bounds
      'ub',          [1000 1000],      ... % Parameter upper bounds
      'options',     opts);                % Options defined above

   % Create a GlobalSearch object
   gs = GlobalSearch;
   
   % Run it, returning the best-fitting parameter values and the negative-
   % log-likelihood returned by the objective function
   [fits(ii,:), nllk] = run(gs,problem);
   
   % Plot using our utility function, which expects RT in msec
   reciprobit_plot(RTs, fits(ii,:), gca, [], colors{ii})
end

% Plot the fits
for xx = 1:2
   
   subplot(2,2,2+xx); cla reset; hold on;
   plot(fits(:,xx), 'ko-');
end
subplot(2,2,3);
axis([1 4 2 8]);
subplot(2,2,4);
axis([1 4 0.5 2.5]);

function reciprobit_plot(rts, fits, ax, xlm, color)
% function reciprobit_plot(rts, fits, ax, xlm, color)
%
% fits are muR, deltaS
% rts are in SECONDS
%
% Reciprobit is reciprocal latency vs cumulative z-score (probit)
%   see Carpenter & Williams 1995
%       Reddi & Carpenter, 2000, 2003
%       http://www.cudos.ac.uk/later.html
%

% Created by jig 09-30-10

if nargin < 1 || isempty(rts)
   return
end

if nargin < 3 || isempty(ax)
   %figure
   ax = gca;
else
   axes(ax);
end

if nargin < 4 || isempty(xlm)
   xlm = [0.1 1];
end

if nargin < 5 || isempty(color)
   color = 'k';
end

% get transformed data
[xs,ys] = reciprobit_xys(rts);

% plot raw data, make axes real purdy
plot(-xs, ys, '.', 'Color', color);

% Label axes with RT (x) and probability (y) values that make sense, but
% corresond to actual -1/RT (x) and probit (y) values that are plotted
XTICK  = 10.^(linspace(log10(xlm(1)), log10(xlm(2)), 4));
YTICK  = [.1 1 5 10 50 90 95 99 99.9];
YTICKi = norminv(YTICK./100,0,1);

set(gca, ...
   'FontSize',     14,                 ...
   'YTick',        YTICKi,             ...
   'YTickLabel',   YTICK,              ...
   'YLim',         YTICKi([1 end]),    ...
   'XTick',        -1./XTICK,          ...
   'XTickLabel',   XTICK,              ...
   'XLim',         -1./XTICK([1 end]));
xlabel('RT (sec)');
ylabel('Probability');

% possibly show fits
if nargin > 1 && ~isempty(fits)
   fxs = -1./xlm;
   plot(fxs, (fxs+fits(1)/fits(2)).*fits(2), '-', 'Color', color);
end
end

function [xvals_, yvals_] = reciprobit_xys(rts)
% function [xvals_, yvals_] = reciprobit_xys(rts)
%
% Compute x (-1./rt) and y (probit) values
%   from RT distribution

if nargin < 1 || isempty(rts)
    xvals_ = [];
    yvals_ = [];
    return
end

% get empirical cumulative distribution
rtss  = sort(rts);
n     = length(rtss);
rtsc  = ((1:n)./n)';
Lgood = diff([-999; rtss])~=0;

xvals_ = 1./rtss(Lgood);
yvals_ = norminv(rtsc(Lgood),0,1);
end