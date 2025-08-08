% ------------------------------------------------------------------------
% Figure 14: p-Curves from IV selection with two-sided tests
% ------------------------------------------------------------------------
clear all
mkdir('Figures/Figures_analytical_examples/two_sided')
rng(12345)

alpha = 0.05; %significance level
N = 1e6; % single MC sample size
H = [0, 1, 2]; % values of h

h_grid = linspace(0, 6, 10000); % grid of h's
cv2 = @(p) norminv(1 - p/2);
p = linspace(0, 1, 203); % x-axis
bound = zeros(1, length(p)); % Upper bound for 2-sided tests
for j = 1:length(p)
    bound(j) = max((exp(h_grid * cv2(p(j)) - h_grid.^2/2) + ...
                    exp(-h_grid * cv2(p(j)) - h_grid.^2/2)) / 2);
end

M = 100; % improve accuracy by averaging 100 Monte Carlos
for j = 1:3
    for m = 0:M
        W1 = randn(N, 1);
        W2 = randn(N, 1);
        W12 = (W1 + W2) / sqrt(2);

        h = H(j);
        T1h  = abs(W1  + h);
        T2h  = abs(W2  + h);
        T12h = abs(W12 + sqrt(2) * h);

        p1  = 2 * normcdf(-T1h);
        p2  = 2 * normcdf(-T2h);
        p12 = 2 * normcdf(-T12h);

        P0     = p12;
        Pr_min = min([p1, p2, p12]')';
        Pr     = P0 .* (P0 <= alpha) + (P0 > alpha) .* Pr_min;

        [centers0, N0]     = Edges_N(P0);
        [centers1, N1]     = Edges_N(Pr);
        [centers1min, N1min] = Edges_N(Pr_min);

        % Initialize storage for first iteration
        if m == 0
            K0 = 400;
            P0_s.(['h' num2str(H(j)) 'iv'])     = zeros(K0, 2);
            Pr_s.(['h' num2str(H(j)) 'iv'])     = zeros(K0, 2);
            Pr_min_s.(['h' num2str(H(j)) 'iv']) = zeros(K0, 2);
        end
        % Average over MC draws
        if m > 0
            P0_s.(['h' num2str(H(j)) 'iv'])     = P0_s.(['h' num2str(H(j)) 'iv'])     + [centers0', N0']     / M;
            Pr_s.(['h' num2str(H(j)) 'iv'])     = Pr_s.(['h' num2str(H(j)) 'iv'])     + [centers1', N1']     / M;
            Pr_min_s.(['h' num2str(H(j)) 'iv']) = Pr_min_s.(['h' num2str(H(j)) 'iv']) + [centers1min', N1min']/ M;
        end
    end
end

% ------------------------------------------------------------------------
% FIGURE 14 (a) [Online Appendix]: IV selection, threshold
% ------------------------------------------------------------------------
close all
figure('Visible', 'off')
hold on
plot(Pr_s.h0iv(:,1), Pr_s.h0iv(:,2), 'LineStyle','--', 'LineWidth',2, 'Color', 'black')
plot(Pr_s.h1iv(:,1), Pr_s.h1iv(:,2), 'LineStyle',':',  'LineWidth',2, 'Color', 'black')
plot(Pr_s.h2iv(:,1), Pr_s.h2iv(:,2), 'LineStyle','-.', 'LineWidth',2, 'Color', 'black')
plot(p, bound,               'LineStyle','-', 'LineWidth',2, 'Color', 'black')
ylim([0 15])
set(gca, 'FontSize', 18)
xlabel('$$p$$',   'FontSize', 25, 'interpreter', 'latex')
ylabel('$$g_2^t(p)$$', 'FontSize', 25, 'interpreter', 'latex')
title('(a) IV selection: threshold', ...
    'fontweight', 'bold', 'FontSize', 20, 'interpreter', 'latex')
legend('$h=0$', '$h=1$', '$h=2$', 'Bound', ...
    'Orientation', 'horizontal', 'interpreter', 'latex')
axes('position', [0.35, 0.34, 0.45, 0.45])
box on
indexOfInterest = (centers0 < 0.2) & (centers0 > 0.0001);
plot(Pr_s.h0iv(indexOfInterest,1), Pr_s.h0iv(indexOfInterest,2), 'LineStyle','--',  'LineWidth',2, 'Color', 'black')
hold on
plot(Pr_s.h1iv(indexOfInterest,1), Pr_s.h1iv(indexOfInterest,2), 'LineStyle',':',   'LineWidth',2, 'Color', 'black')
plot(Pr_s.h2iv(indexOfInterest,1), Pr_s.h2iv(indexOfInterest,2), 'LineStyle','-.',  'LineWidth',2, 'Color', 'black')
plot(p(indexOfInterest), bound(indexOfInterest), 'LineStyle','-', 'LineWidth',2, 'Color', 'black')
xlim([0 0.1])
ylim([0, 10])
saveas(gcf, 'Figures/Figures_analytical_examples/two_sided/Fig14a', 'epsc')
close all

% ------------------------------------------------------------------------
% FIGURE 14 (b) [Online Appendix]: IV selection, minimum
% ------------------------------------------------------------------------
close all
figure('Visible', 'off')
hold on
plot(Pr_min_s.h0iv(:,1), Pr_min_s.h0iv(:,2), 'LineStyle','--', 'LineWidth',2, 'Color', 'black')
plot(Pr_min_s.h1iv(:,1), Pr_min_s.h1iv(:,2), 'LineStyle',':',  'LineWidth',2, 'Color', 'black')
plot(Pr_min_s.h2iv(:,1), Pr_min_s.h2iv(:,2), 'LineStyle','-.', 'LineWidth',2, 'Color', 'black')
plot(p, bound,                 'LineStyle','-', 'LineWidth',2, 'Color', 'black')
ylim([0 15])
set(gca, 'FontSize', 18)
xlabel('$$p$$',   'FontSize', 25, 'interpreter', 'latex')
ylabel('$$g_2^m(p)$$', 'FontSize', 25, 'interpreter', 'latex')
title('(b) IV selection: minimum', ...
    'fontweight', 'bold', 'FontSize', 20, 'interpreter', 'latex')
legend('$h=0$', '$h=1$', '$h=2$', 'Bound', ...
    'Orientation', 'horizontal', 'interpreter', 'latex')
axes('position', [0.35, 0.34, 0.45, 0.45])
box on
indexOfInterest = (centers0 < 0.2) & (centers0 > 0.0001);
plot(Pr_min_s.h0iv(indexOfInterest,1), Pr_min_s.h0iv(indexOfInterest,2), 'LineStyle','--',  'LineWidth',2, 'Color', 'black')
hold on
plot(Pr_min_s.h1iv(indexOfInterest,1), Pr_min_s.h1iv(indexOfInterest,2), 'LineStyle',':',   'LineWidth',2, 'Color', 'black')
plot(Pr_min_s.h2iv(indexOfInterest,1), Pr_min_s.h2iv(indexOfInterest,2), 'LineStyle','-.',  'LineWidth',2, 'Color', 'black')
plot(p(indexOfInterest), bound(indexOfInterest), 'LineStyle','-', 'LineWidth',2, 'Color', 'black')
xlim([0 0.1])
ylim([0, 10])
saveas(gcf, 'Figures/Figures_analytical_examples/two_sided/Fig14b', 'epsc')
close all

% ------------------------------------------------------------------------
% FIGURE 14 (c) [Online Appendix]: IV selection, threshold vs. minimum
% ------------------------------------------------------------------------
close all
figure('Visible', 'off')
hold on
plot(Pr_s.h1iv(:,1),      Pr_s.h1iv(:,2),      'LineStyle','--', 'LineWidth',2, 'Color', 'black')
plot(Pr_min_s.h1iv(:,1),  Pr_min_s.h1iv(:,2),  'LineStyle',':',  'LineWidth',2, 'Color', 'black')
plot(P0_s.h1iv(:,1),      P0_s.h1iv(:,2),      'LineStyle','-.', 'LineWidth',2, 'Color', 'black')
plot(p, bound,            'LineStyle','-',     'LineWidth',2, 'Color', 'black')
ylim([0 15])
set(gca, 'FontSize', 18)
xlabel('$$p$$',   'FontSize', 25, 'interpreter', 'latex')
ylabel('$$g_2(p)$$', 'FontSize', 25, 'interpreter', 'latex')
title('(c) IV selection: threshold vs. minimum, $h = 1$', ...
    'fontweight', 'bold', 'FontSize', 20, 'interpreter', 'latex')
legend('Threshold', 'Minimum','No $p$-hacking','Bound', ...
    'Orientation','horizontal', 'interpreter', 'latex', 'NumColumns',2)
axes('position', [0.35, 0.34, 0.45, 0.45])
box on
indexOfInterest = (centers0 < 0.2) & (centers0 > 0.0001);
plot(Pr_s.h1iv(indexOfInterest,1),      Pr_s.h1iv(indexOfInterest,2),      'LineStyle','--',  'LineWidth',2, 'Color', 'black')
hold on
plot(Pr_min_s.h1iv(indexOfInterest,1),  Pr_min_s.h1iv(indexOfInterest,2),  'LineStyle',':',   'LineWidth',2, 'Color', 'black')
plot(P0_s.h1iv(indexOfInterest,1),      P0_s.h1iv(indexOfInterest,2),      'LineStyle','-.',  'LineWidth',2, 'Color', 'black')
plot(p(indexOfInterest), bound(indexOfInterest), 'LineStyle','-', 'LineWidth',2, 'Color', 'black')
xlim([0 0.1])
ylim([0, 15])
saveas(gcf, 'Figures/Figures_analytical_examples/two_sided/Fig14c', 'epsc')
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