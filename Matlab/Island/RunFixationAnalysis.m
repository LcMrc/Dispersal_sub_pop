% ========================================================================
% RunFixationAnalysis.m
% Author: Loïc Marrec
%
% Description:
%   Runs Gillespie simulations of deme fixation dynamics in the island
%   model. Computes:
%       - Fixation probability (pfix)
%       - Mean number of local fixations (nfixmean)
%       - Mean fixation time (tfixmean)
%   and generates summary plots for analysis.
%
% Dependencies:
%   - GillespieDemeFixation.m
%   - GillespieDemeFixationConditioned.m
%
% ========================================================================

%% Initialization

clear; close all; clc;

% Simulation parameters
Nit = 1e2;                 % Number of stochastic replicates
D = 5;                     % Number of demes
K = 1e2;                   % Carrying capacity 

% Wild-type parameters
bW = 1;                    % Birth rate
dW = 0.1;                  % Death rate
mW = 1e-6;                 % Dispersal rate

% Mutant parameters
bM = 0.95 : 0.01 : 1.1;     % Range of mutant birth rates
dM = 0.1;                   % Death rate
mM = 2e-6;                  % Dispersal rate
alpha = mM / mW;            % Dispersal rate ratio (mutant/wild-type)

% Derived quantities
N = K * (1 - dW / bW);      % Wild-type deme equilibrium size 

% Preallocate result arrays
DMlist = NaN(length(bM), Nit);
nfixlist = NaN(length(bM), Nit);
tfixlist = NaN(length(bM), Nit);
pfix = NaN(length(bM), 1);

%% Main simulation loop

for i = 1 : length(bM)

    % General simulation: fixation OR loss
    [DMlist(i, :), ~, ~] = GillespieDemeFixation(Nit, D, K, bW, dW, bM(i), dM, mW, mM);

    % Compute fixation probability
    pfix(i) = sum(DMlist(i, :) ~= 0) / Nit;

    % Conditioned simulation: fixation only
    [~, nfixlist(i, :), tfixlist(i, :)] = GillespieDemeFixationConditioned(Nit, D, K, bW, dW, bM(i), dM, mW, mM);

    % Progress display
    fprintf('Completed %d / %d\n', i, length(bM));

end

%% Statistical analysis

% Mean and confidence intervals for number of local fixations
nfixmean = nanmean(nfixlist, 2);
nfixstd = nanstd(nfixlist, 0, 2);
nfixci = 1.96 .* nfixstd ./ sqrt(Nit);

% Mean and confidence intervals for fixation time
tfixmean = nanmean(tfixlist, 2);
tfixstd = nanstd(tfixlist, 0, 2);
tfixci = 1.96 .* tfixstd ./ sqrt(Nit);

% Selection coefficients
sb = bM ./ bW - 1;
sd = dW ./ dM - 1;
s = sb + sd;

%% Plot settings

plotStyle = {'LineStyle', 'none', 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k'};
axisFontSize = 16;

%% Figure 1: Fixation Probability 

figure(1); clf; hold on;
plot(s, pfix, plotStyle{:});
plot([0 0], [0 1], '--k', 'LineWidth', 1.5);
plot([min(s) max(s)], [1/D 1/D], '--k', 'LineWidth', 1.5);
hold off;

xlabel('Selection coefficient, s', 'FontSize', axisFontSize);
ylabel('Single-deme fixation probability, U', 'FontSize', axisFontSize);
set(gca, 'FontName', 'Arial', 'FontSize', axisFontSize, ...
    'Box', 'off', 'TickDir', 'out', 'LineWidth', 1, ...
    'XMinorTick', 'on', 'YMinorTick', 'on');
xlim([min(s) max(s)]);
ylim([0 1]);
axis tight;

%% Figure 2: Number of Local Fixations 

figure(2); clf; hold on;
plot(s, nfixmean, plotStyle{:});
plot([0 0], [min(nfixmean) max(nfixmean)], '--k', 'LineWidth', 1.5);
hold off;

xlabel('Selection coefficient, s', 'FontSize', axisFontSize);
ylabel('Number of local fixations, n_{fix}', 'FontSize', axisFontSize);
set(gca, 'FontName', 'Arial', 'FontSize', axisFontSize, ...
    'Box', 'off', 'TickDir', 'out', 'LineWidth', 1, ...
    'XMinorTick', 'on', 'YMinorTick', 'on');
xlim([min(s) max(s)]);
ylim([min(nfixmean) max(nfixmean)]);
axis tight;

%% Figure 3: Fixation Time 

figure(3); clf; hold on;
plot(s, tfixmean, plotStyle{:});
plot([0 0], [min(tfixmean) max(tfixmean)], '--k', 'LineWidth', 1.5);
hold off;

xlabel('Selection coefficient, s', 'FontSize', axisFontSize);
ylabel('Global fixation time, t_{fix}', 'FontSize', axisFontSize);
set(gca, 'FontName', 'Arial', 'FontSize', axisFontSize, ...
    'Box', 'off', 'TickDir', 'out', 'LineWidth', 1, ...
    'XMinorTick', 'on', 'YMinorTick', 'on');
xlim([min(s) max(s)]);
ylim([min(tfixmean) max(tfixmean)]);
axis tight;

