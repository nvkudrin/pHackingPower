% ------------------------------------------------------------------------
% Figure 17: p-Curves from lag length selection with two-sided tests
% ------------------------------------------------------------------------
clear all
mkdir('Figures/Figures_analytical_examples/two_sided')
rng(12345)

alpha = 0.05; % significance level
Mc = 1e7;    % Number of Monte Carlo simulations
N = 200;     % Sample size
H = [0, 1, 2]; % h-values

Pr_min = zeros(Mc, 3); % Stores minimum p-values for each h
Pr     = zeros(Mc, 3); % Stores thresholded p-values for each h

% Main Monte Carlo loop
for m = 1:Mc
    U = randn(N,1);
    for j = 1:3
        h = H(j);
        T0 = abs(h + sqrt(N)*mean(U)); % t-stat for h
        rho_hat = mean((U(2:end)-mean(U)).*(U(1:(end-1))-mean(U)));
        omega2_hat = (1 + rho_hat); % plug-in estimator for variance inflation

        % Handle possible negative omega2_hat
        if omega2_hat <= 0
            T1 = T0;
        else
            T1 = T0 / sqrt(omega2_hat);
        end

        p0 = 2 * normcdf(-T0);
        p1 = 2 * normcdf(-T1);

        Pr_min(m,j) = min(p0, p1);
        Pr(m,j) = p0 * (p0 < alpha) + (p0 > alpha) * min(p0, p1);
    end
end

% Histogram centers and counts for p-value curves
[centers10, N10] = Edges_N(Pr(:,1));
[centers11, N11] = Edges_N(Pr(:,2));
[centers12, N12] = Edges_N(Pr(:,3));
[centers1min0, N1min0] = Edges_N(Pr_min(:,1));
[centers1min1, N1min1] = Edges_N(Pr_min(:,2));
[centers1min2, N1min2] = Edges_N(Pr_min(:,3));

centers0 = centers10; % use for inset/indexing

% Compute analytical bound for two-sided case
h_grid = linspace(0, 6, 10000);
cv2 = @(p) norminv(1-p/2);
p = centers10;
bound = zeros(size(p));
for j = 1:length(p)
    bound(j) = max((exp(h_grid*cv2(p(j)) - h_grid.^2/2) + ...
                    exp(-h_grid*cv2(p(j)) - h_grid.^2/2))/2);
end

% ------------------------------------------------------------------------
% FIGURE 17 (a) [Online Appendix]: Lag length selection, threshold
% ------------------------------------------------------------------------
close all
figure('Visible', 'off')
hold on
plot(centers10, N10, 'LineStyle','--', 'LineWidth',2, 'Color', 'black')
plot(centers11, N11, 'LineStyle',':',  'LineWidth',2, 'Color', 'black')
plot(centers12, N12, 'LineStyle','-.', 'LineWidth',2, 'Color', 'black')
plot(p, bound,    'LineStyle','-',  'LineWidth',2, 'Color', 'black')
ylim([0 15])
set(gca, 'FontSize', 18)
xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('$$g_3^t(p)$$', 'FontSize',25, 'interpreter', 'latex')
title('(a) Lag length selection: threshold', ...
    'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
legend('$h=0$', '$h=1$', '$h=2$', 'Bound', ...
    'Orientation','horizontal', 'interpreter', 'latex')

axes('position',[.35 .34 .45 .45])
box on
indexOfInterest = (centers0 < 0.2) & (centers0 > 0.0001);
plot(centers10(indexOfInterest), N10(indexOfInterest), 'LineStyle','--', 'LineWidth',2, 'Color', 'black')
hold on
plot(centers11(indexOfInterest), N11(indexOfInterest), 'LineStyle',':',  'LineWidth',2, 'Color', 'black')
plot(centers12(indexOfInterest), N12(indexOfInterest), 'LineStyle','-.', 'LineWidth',2, 'Color', 'black')
plot(p(indexOfInterest), bound(indexOfInterest), 'LineStyle','-', 'LineWidth',2, 'Color', 'black')
ylim([0 10])
xlim([0 0.1])
saveas(gcf, 'Figures/Figures_analytical_examples/two_sided/Fig17a', 'epsc')
close all

% ------------------------------------------------------------------------
% FIGURE 17 (b) [Online Appendix]: Lag length selection, minimum
% ------------------------------------------------------------------------
close all
figure('Visible', 'off')
hold on
plot(centers1min0, N1min0, 'LineStyle','--', 'LineWidth',2, 'Color', 'black')
plot(centers1min1, N1min1, 'LineStyle',':',  'LineWidth',2, 'Color', 'black')
plot(centers1min2, N1min2, 'LineStyle','-.', 'LineWidth',2, 'Color', 'black')
plot(p, bound,         'LineStyle','-',  'LineWidth',2, 'Color', 'black')
ylim([0 15])
set(gca, 'FontSize', 18)
xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('$$g_3^m(p)$$', 'FontSize',25, 'interpreter', 'latex')
title('(b) Lag length selection: minimum', ...
    'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
legend('$h=0$', '$h=1$', '$h=2$', 'Bound', ...
    'Orientation','horizontal', 'interpreter', 'latex')

axes('position',[.35 .34 .45 .45])
box on
indexOfInterest = (centers0 < 0.6) & (centers0 > 0.0001);
plot(centers1min0(indexOfInterest), N1min0(indexOfInterest), 'LineStyle','--', 'LineWidth',2, 'Color', 'black')
hold on
plot(centers1min1(indexOfInterest), N1min1(indexOfInterest), 'LineStyle',':',  'LineWidth',2, 'Color', 'black')
plot(centers1min2(indexOfInterest), N1min2(indexOfInterest), 'LineStyle','-.', 'LineWidth',2, 'Color', 'black')
plot(p(indexOfInterest), bound(indexOfInterest), 'LineStyle','-', 'LineWidth',2, 'Color', 'black')
ylim([0 2])
xlim([0 0.6])
saveas(gcf, 'Figures/Figures_analytical_examples/two_sided/Fig17b', 'epsc')
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