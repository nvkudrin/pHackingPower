% ------------------------------------------------------------------------
% Figure 13: p-Curves from covariate selection with two-sided tests
% ------------------------------------------------------------------------
clear all
mkdir('Figures/Figures_analytical_examples/two_sided')
rng(12345)

h_grid = linspace(0, 6, 10000); % grid of values of h
cv2 = @(p) norminv(1 - p/2);
p = linspace(0, 1, 203); % x-axis
bound = zeros(1, length(p)); %Upper bound for 2-sided t-tests
for j = 1:length(p)
    bound(j) = max((exp(h_grid * cv2(p(j)) - h_grid.^2/2) + ...
                    exp(-h_grid * cv2(p(j)) - h_grid.^2/2)) / 2);
end

alpha = 0.05; %significance level
N = 1e6; % single MC sample size

gamma = [0.1, 0.5, 0.9]; %values of gamma
Rho  = 1 - gamma.^2;
H = [0, 1, 2]; %values of h
P0_s = struct(); % Null
Pr_s = struct(); %p-hacked thresholding
Pr_min_s = struct(); %p-hacked minimum
M = 100; % improve accuracy by averaging 100 Monte Carlos

for k = 1:length(gamma)
    rho = Rho(k);
    for j = 1:3
        for m = 0:M
            T1 = randn(N, 1);
            T2 = rho * T1 + sqrt(1 - rho^2) * randn(N, 1);
            h = H(j);
            T1h = abs(T1 + h);
            T2h = abs(T2 + h);
            p1 = 2 * normcdf(-T1h);
            p2 = 2 * normcdf(-T2h);

            P0 = p1;
            Pr_min = (p2 .* (p2 <= p1) + p1 .* (p1 < p2));
            Pr = p1 .* (p1 <= alpha) + (p1 > alpha) .* Pr_min;

            [centers0, N0] = Edges_N(P0);
            [centers1, N1] = Edges_N(Pr);
            [centers1min, N1min] = Edges_N(Pr_min);

            % Initialize storage arrays for first iteration
            if m == 0
                K0 = 400;
                P0_s.(['h' num2str(H(j)) 'g' num2str(gamma(k)*10)]) = zeros(K0, 2);
                Pr_s.(['h' num2str(H(j)) 'g' num2str(gamma(k)*10)]) = zeros(K0, 2);
                Pr_min_s.(['h' num2str(H(j)) 'g' num2str(gamma(k)*10)]) = zeros(K0, 2);
            end

            % Average over MC runs
            if m > 0
                P0_s.(['h' num2str(H(j)) 'g' num2str(gamma(k)*10)]) = ...
                    P0_s.(['h' num2str(H(j)) 'g' num2str(gamma(k)*10)]) + [centers0', N0']/M;
                Pr_s.(['h' num2str(H(j)) 'g' num2str(gamma(k)*10)]) = ...
                    Pr_s.(['h' num2str(H(j)) 'g' num2str(gamma(k)*10)]) + [centers1', N1']/M;
                Pr_min_s.(['h' num2str(H(j)) 'g' num2str(gamma(k)*10)]) = ...
                    Pr_min_s.(['h' num2str(H(j)) 'g' num2str(gamma(k)*10)]) + [centers1min', N1min']/M;
            end
        end
    end
end

close all

% ------------------------------------------------------------------------
% FIGURE 13 (a) [Online Appendix]
% ------------------------------------------------------------------------
figure('Visible', 'off')
hold on
plot(Pr_s.h0g5(:,1), Pr_s.h0g5(:,2), 'LineStyle','--', 'LineWidth',2, 'Color', 'black')
plot(Pr_s.h1g5(:,1), Pr_s.h1g5(:,2), 'LineStyle',':',  'LineWidth',2, 'Color', 'black')
plot(Pr_s.h2g5(:,1), Pr_s.h2g5(:,2), 'LineStyle','-.', 'LineWidth',2, 'Color', 'black')
plot(p, bound,              'LineStyle','-',  'LineWidth',2, 'Color', 'black')
ylim([0 15])
set(gca, 'FontSize', 18)
xlabel('$$p$$', 'FontSize', 25, 'interpreter', 'latex')
ylabel('$$g_1^t(p)$$', 'FontSize', 25, 'interpreter', 'latex')
title('(a) Covariate selection: threshold, $\gamma=0.5$', ...
    'fontweight', 'bold', 'FontSize', 20, 'interpreter', 'latex')
legend('$h=0$', '$h=1$', '$h=2$', 'Bound', ...
    'Orientation','horizontal', 'interpreter', 'latex')
axes('position', [0.35, 0.34, 0.45, 0.45])
box on
indexOfInterest = (centers0 < 0.2) & (centers0 > 0.0001);
plot(Pr_s.h0g5(indexOfInterest,1), Pr_s.h0g5(indexOfInterest,2), 'LineStyle','--',  'LineWidth',2, 'Color', 'black')
hold on
plot(Pr_s.h1g5(indexOfInterest,1), Pr_s.h1g5(indexOfInterest,2), 'LineStyle',':',   'LineWidth',2, 'Color', 'black')
plot(Pr_s.h2g5(indexOfInterest,1), Pr_s.h2g5(indexOfInterest,2), 'LineStyle','-.',  'LineWidth',2, 'Color', 'black')
plot(p(indexOfInterest), bound(indexOfInterest), 'LineStyle','-', 'LineWidth',2, 'Color', 'black')
xlim([0 0.1])
ylim([0 12])
saveas(gcf, 'Figures/Figures_analytical_examples/two_sided/Fig13a', 'epsc')
close all

% ------------------------------------------------------------------------
% FIGURE 13 (b) [Online Appendix]
% ------------------------------------------------------------------------
close all
figure('Visible', 'off')
hold on
plot(Pr_s.h1g1(:,1), Pr_s.h1g1(:,2), 'LineStyle','--', 'LineWidth',2, 'Color', 'black')
plot(Pr_s.h1g5(:,1), Pr_s.h1g5(:,2), 'LineStyle',':',  'LineWidth',2, 'Color', 'black')
plot(Pr_s.h1g9(:,1), Pr_s.h1g9(:,2), 'LineStyle','-.', 'LineWidth',2, 'Color', 'black')
plot(p, bound,               'LineStyle','-', 'LineWidth',2, 'Color', 'black')
ylim([0 15])
set(gca, 'FontSize', 18)
xlabel('$$p$$', 'FontSize', 25, 'interpreter', 'latex')
ylabel('$$g_1^t(p)$$', 'FontSize', 25, 'interpreter', 'latex')
title('(b) Covariate selection: threshold, $h=1$', ...
    'fontweight', 'bold', 'FontSize', 20, 'interpreter', 'latex')
legend('$\gamma=0.1$', '$\gamma=0.5$', '$\gamma=0.9$', 'Bound', ...
    'Orientation','horizontal', 'interpreter', 'latex')
axes('position', [0.35, 0.34, 0.45, 0.45])
box on
indexOfInterest = (centers0 < 0.2) & (centers0 > 0.0001);
plot(Pr_s.h1g1(indexOfInterest,1), Pr_s.h1g1(indexOfInterest,2), 'LineStyle','--',  'LineWidth',2, 'Color', 'black')
hold on
plot(Pr_s.h1g5(indexOfInterest,1), Pr_s.h1g5(indexOfInterest,2), 'LineStyle',':',   'LineWidth',2, 'Color', 'black')
plot(Pr_s.h1g9(indexOfInterest,1), Pr_s.h1g9(indexOfInterest,2), 'LineStyle','-.',  'LineWidth',2, 'Color', 'black')
plot(p(indexOfInterest), bound(indexOfInterest), 'LineStyle','-', 'LineWidth',2, 'Color', 'black')
xlim([0 0.1])
ylim([0 15])
saveas(gcf, 'Figures/Figures_analytical_examples/two_sided/Fig13b', 'epsc')
close all

% ------------------------------------------------------------------------
% FIGURE 13 (c) [Online Appendix]
% ------------------------------------------------------------------------
close all
figure('Visible', 'off')
hold on
plot(Pr_s.h1g5(:,1),    Pr_s.h1g5(:,2),    'LineStyle','--', 'LineWidth',2, 'Color', 'black')
plot(Pr_min_s.h1g5(:,1), Pr_min_s.h1g5(:,2), 'LineStyle',':', 'LineWidth',2, 'Color', 'black')
plot(P0_s.h1g5(:,1),     P0_s.h1g5(:,2),    'LineStyle','-.', 'LineWidth',2, 'Color', 'black')
plot(p, bound,           'LineStyle','-',   'LineWidth',2, 'Color', 'black')
ylim([0 15])
set(gca, 'FontSize', 18)
xlabel('$$p$$', 'FontSize', 25, 'interpreter', 'latex')
ylabel('$$g_1(p)$$', 'FontSize', 25, 'interpreter', 'latex')
title('(c) Covariate selection: threshold vs. minimum, $h = 1$', ...
    'fontweight', 'bold', 'FontSize', 20, 'interpreter', 'latex')
legend('Threshold', 'Minimum','No $p$-hacking','Bound', ...
    'Orientation','horizontal', 'interpreter', 'latex', 'NumColumns',2)
axes('position', [0.35, 0.34, 0.45, 0.45])
box on
indexOfInterest = (centers0 < 0.2) & (centers0 > 0.0001);
plot(Pr_s.h1g5(indexOfInterest,1),    Pr_s.h1g5(indexOfInterest,2),    'LineStyle','--',  'LineWidth',2, 'Color', 'black')
hold on
plot(Pr_min_s.h1g5(indexOfInterest,1), Pr_min_s.h1g5(indexOfInterest,2), 'LineStyle',':',  'LineWidth',2, 'Color', 'black')
plot(P0_s.h1g5(indexOfInterest,1),     P0_s.h1g5(indexOfInterest,2),    'LineStyle','-.',  'LineWidth',2, 'Color', 'black')
plot(p(indexOfInterest), bound(indexOfInterest), 'LineStyle','-', 'LineWidth',2, 'Color', 'black')
xlim([0 0.1])
ylim([0 15])
saveas(gcf, 'Figures/Figures_analytical_examples/two_sided/Fig13c', 'epsc')
close all

    % --- Local function for computing histogram as PDF ---
    function [centers, N] = Edges_N(P)
    %EDGESN  Compute histogram for Monte Carlo p-curves
        %   [centers, N] = Edges_N(P)
        %   Computes a normalized histogram (PDF) of p-values using 400 bins on [0,1].
        %   Returns the bin centers and bin counts (N).
    
        edges = linspace(0, 1, 401);  % 400 bins between 0 and 1
        [N, edges] = histcounts(P, edges, 'Normalization', 'pdf');
    
        % Convert bin edges to bin centers
        centers = edges(2:end) - (edges(2) - edges(1))/2;
    end