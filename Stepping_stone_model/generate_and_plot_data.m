% This code was written by Loïc Marrec 

clear all;
close all;

Nit = 1e3;                                  % Number of stochastic replicates
D = 5;                                      % Number of demes
K = 1e2;                                    % Carrying capacity
fA = 1;                                     % Wild-type birth rate
gA = .1;                                    % Wild-type death rate
fB = .95 : .01 : 1.1;                       % Mutant birth rate
gB = .1;                                    % Mutant death rate
tmig = 9000;                                % Total time between two dispersal events
N = K*(1-gA/fA);                            % Wild-type equilibrium size
alpha = 1/5;                                % Ratio of mutant and wild-type dispersal rates
mA = alpha/(tmig*(1+alpha)*N*D);            % Wild-type dispersal rate            
mB = 1/(tmig*(1+alpha)*N*D);                % Mutant dispersal rate
XBlist = NaN(length(fB), Nit);              % List of final number of mutants 
nfixlist = NaN(length(fB), Nit);            % List of number of local fixations
tfixlist = NaN(length(fB), Nit);            % List of fixation times
pfix = NaN(length(fB), 1);                  % List of fixation probabilities

for i = 1 : length(fB)
    
    [XBlist(i, :), ~, ~] = Gillespie_fct(Nit, D, K, fA, gA, fB(i), gB, mA, mB);
    
    pfix(i) = length(XBlist(i, XBlist(i, :) ~= 0))/Nit;     % Compute the fixation probability
    
    [~, nfixlist(i, :), tfixlist(i, :)] = Gillespie_fix_fct(Nit, D, K, fA, gA, fB(i), gB, mA, mB);
    
    Say = sprintf('%d/%d', i, length(fB));
    disp(Say)
    
end

nfixmean = nanmean(nfixlist, 2);      % Compute the mean number of local fixations
nfixstd = nanstd(nfixlist, 0, 2);     % Compute the standard deviation of the number of local fixations
nfixci = 1.96.*nfixstd./Nit;          % Compute the 95% interval confidence of the number of local fixations

tfixmean = nanmean(tfixlist, 2);      % Compute the mean global fixation time
tfixstd = nanstd(tfixlist, 0, 2);     % Compute the standard deviation of the global fixation time
tfixci = 1.96.*tfixstd./Nit;          % Compute the 95% interval confidence of the global fixation time

figure(1)
hold on  

    plot(fB-fA, pfix, 'LineStyle', 'None', 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k');
    plot([log(alpha)/K log(alpha)/K], [0 1], 'LineStyle', ':', 'LineWidth', 1.5, 'Color', 'k');
    plot([0 0], [0 1], 'LineStyle', '--', 'LineWidth', 1.5, 'Color', 'k');
    plot([min(fB)-fA max(fB)-fA], [1/D 1/D], 'LineStyle', '--', 'LineWidth', 1.5, 'Color', 'k');

hold off
hXLabel = xlabel('Selection coefficient, s', 'Color', 'k');
hYLabel = ylabel('Fixation probability, U', 'Color', 'k');
set( gca                       , ...
    'FontName'   , 'Arial'   , 'FontSize'   , 16);
set([hXLabel, hYLabel], ...
    'FontName'   , 'Arial'   , 'FontSize'   , 16);
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'LineWidth'   , 1         );
ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
axis tight
xlim([min(fB)-fA max(fB)-fA])
ylim([0 1])

figure(2)
hold on  

    plot(fB-fA, nfixmean, 'LineStyle', 'None', 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k');
    plot([log(alpha)/K log(alpha)/K], [min(nfixmean) max(nfixmean)], 'LineStyle', ':', 'LineWidth', 1.5, 'Color', 'k');
    plot([0 0], [min(nfixmean) max(nfixmean)], 'LineStyle', '--', 'LineWidth', 1.5, 'Color', 'k');

hold off
hXLabel = xlabel('Selection coefficient, s', 'Color', 'k');
hYLabel = ylabel('Number of local fixations, n_{fix}', 'Color', 'k');
set( gca                       , ...
    'FontName'   , 'Arial'   , 'FontSize'   , 16);
set([hXLabel, hYLabel], ...
    'FontName'   , 'Arial'   , 'FontSize'   , 16);
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'LineWidth'   , 1         );
ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
axis tight
xlim([min(fB)-fA max(fB)-fA])
ylim([min(nfixmean) max(nfixmean)])

figure(3)
hold on  

    plot(fB-fA, tfixmean, 'LineStyle', 'None', 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k');
    plot([log(alpha)/K log(alpha)/K], [min(tfixmean) max(tfixmean)], 'LineStyle', ':', 'LineWidth', 1.5, 'Color', 'k');
    plot([0 0], [min(tfixmean) max(tfixmean)], 'LineStyle', '--', 'LineWidth', 1.5, 'Color', 'k');

hold off
hXLabel = xlabel('Selection coefficient, s', 'Color', 'k');
hYLabel = ylabel('Fixation time, t_{fix}', 'Color', 'k');
set( gca                       , ...
    'FontName'   , 'Arial'   , 'FontSize'   , 16);
set([hXLabel, hYLabel], ...
    'FontName'   , 'Arial'   , 'FontSize'   , 16);
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'LineWidth'   , 1         );
ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
axis tight
xlim([min(fB)-fA max(fB)-fA])
ylim([min(tfixmean) max(tfixmean)])