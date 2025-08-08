% ------------------------------------------------------------------------
% Figure 1: Covariate Selection p-curves, one-sided tests
% ------------------------------------------------------------------------
clear all
mkdir('Figures/Figures_analytical_examples/one_sided')
rng(12345)

alpha = 0.05; % significance level
z = @(h, p) norminv(1 - p) - h;
gamma = [0.5, 0.5, 0.5, 0.1, 0.5, 0.9]; % first 3 for (a), last 3 for (b)
H     = [0, 1, 2, 1, 1, 1]; % first 3 for (a), last 3 for (b)
P = linspace(0.0001, 0.9999, 1000); % x-axis
pcurves = zeros(length(P), 6); % store p-curves here
Bound = exp(z(0, P).^2 / 2) .* (P <= 0.5) + (P > 0.5); % Upper bound

% Calculate 6 p-curves for combinations of gamma and h: first 3 for (a), last 3 for (b)
for d = 1:6
    rho = 1 - gamma(d)^2;
    h   = H(d);
    for j = 1:length(P)
        p = P(j);
        Upsilon_covsel_t = (1 + normcdf((z(h, alpha) - rho * z(h, p)) / sqrt(1 - rho^2))) * (p <= alpha) ...
            + 2 * normcdf(z(h, p) * sqrt((1 - rho) / (1 + rho))) * (p > alpha);
        pcurves(j, d) = exp(h * z(0, p) - h^2 / 2) * Upsilon_covsel_t;
    end
end

% === Figure 1 (a) [Main text] ===========================================
figure('Visible','off')
plot(P, pcurves(:, 1), '--',  'LineWidth', 2, 'Color', 'black'); hold on
plot(P, pcurves(:, 2), ':',   'LineWidth', 2, 'Color', 'black');
plot(P, pcurves(:, 3), '-.',  'LineWidth', 2, 'Color', 'black');
plot(P, Bound,         '-',   'LineWidth', 2, 'Color', 'black');
ylim([0 15])
set(gca, 'FontSize', 18)
xlabel('$$p$$', 'FontSize', 25, 'Interpreter', 'latex')
ylabel('$$g_1^t(p)$$', 'FontSize', 25, 'Interpreter', 'latex')
title('(a) Covariate selection: threshold, $\gamma=0.5$', ...
    'FontWeight', 'bold', 'FontSize', 20, 'Interpreter', 'latex')
legend({'$h=0$', '$h=1$', '$h=2$', 'Bound'}, ...
    'Orientation', 'horizontal', 'Interpreter', 'latex')
% Inset (zoom-in) for small p-values
axes('Position', [.35 .34 .45 .45])
box on
indexOfInterest = (P < 0.2) & (P > 0.0001);
plot(P(indexOfInterest), pcurves(indexOfInterest, 1), '--',  'LineWidth', 2, 'Color', 'black'); hold on
plot(P(indexOfInterest), pcurves(indexOfInterest, 2), ':',   'LineWidth', 2, 'Color', 'black');
plot(P(indexOfInterest), pcurves(indexOfInterest, 3), '-.',  'LineWidth', 2, 'Color', 'black');
plot(P(indexOfInterest), Bound(indexOfInterest), '-', 'LineWidth', 2, 'Color', 'black');
ylim([0 12])
xlim([0 0.1])
saveas(gcf, 'Figures/Figures_analytical_examples/one_sided/Fig1a', 'epsc')
close all

% === Figure 1 (b) [Main text] ===========================================
figure('Visible','off')
plot(P, pcurves(:, 4), '--',  'LineWidth', 2, 'Color', 'black'); hold on
plot(P, pcurves(:, 5), ':',   'LineWidth', 2, 'Color', 'black');
plot(P, pcurves(:, 6), '-.',  'LineWidth', 2, 'Color', 'black');
plot(P, Bound,         '-',   'LineWidth', 2, 'Color', 'black');
ylim([0 15])
set(gca, 'FontSize', 18)
xlabel('$$p$$', 'FontSize', 25, 'Interpreter', 'latex')
ylabel('$$g_1^t(p)$$', 'FontSize', 25, 'Interpreter', 'latex')
title('(b) Covariate selection: threshold, $h=1$', ...
    'FontWeight', 'bold', 'FontSize', 20, 'Interpreter', 'latex')
legend({'$\gamma=0.1$', '$\gamma=0.5$', '$\gamma=0.9$', 'Bound'}, ...
    'Orientation', 'horizontal', 'Interpreter', 'latex')
% Inset (zoom-in) for small p-values
axes('Position', [.35 .34 .45 .45])
box on
indexOfInterest = (P < 0.6) & (P > 0.0001);
plot(P(indexOfInterest), pcurves(indexOfInterest, 4), '--',  'LineWidth', 2, 'Color', 'black'); hold on
plot(P(indexOfInterest), pcurves(indexOfInterest, 5), ':',   'LineWidth', 2, 'Color', 'black');
plot(P(indexOfInterest), pcurves(indexOfInterest, 6), '-.',  'LineWidth', 2, 'Color', 'black');
plot(P(indexOfInterest), Bound(indexOfInterest), '-', 'LineWidth', 2, 'Color', 'black');
ylim([0 15])
xlim([0 0.1])
saveas(gcf, 'Figures/Figures_analytical_examples/one_sided/Fig1b', 'epsc')
close all

% (c) Covariate selection: threshold vs minimum vs no p-hacking 
h = 1;
gamma = 0.5;
rho = 1 - gamma^2;
pcurves = zeros(length(P), 3);
Bound = exp(z(0, P).^2 / 2) .* (P <= 0.5) + (P > 0.5);

for j = 1:length(P)
    p = P(j);
    Upsilon_covsel_t = (1 + normcdf((z(h, alpha) - rho * z(h, p)) / sqrt(1 - rho^2))) * (p <= alpha) ...
        + 2 * normcdf(z(h, p) * sqrt((1 - rho) / (1 + rho))) * (p > alpha);
    Upsilon_covsel_m = 2 * normcdf(z(h, p) * sqrt((1 - rho) / (1 + rho)));
    pcurves(j, 1) = exp(h * z(0, p) - h^2 / 2) * Upsilon_covsel_t;
    pcurves(j, 2) = exp(h * z(0, p) - h^2 / 2) * Upsilon_covsel_m;
    pcurves(j, 3) = exp(h * z(0, p) - h^2 / 2); % no p-hacking
end

% === Figure 1 (c) [Main text] ===========================================
figure('Visible','off')
plot(P, pcurves(:, 1), '--',  'LineWidth', 2, 'Color', 'black'); hold on
plot(P, pcurves(:, 2), ':',   'LineWidth', 2, 'Color', 'black');
plot(P, pcurves(:, 3), '-.',  'LineWidth', 2, 'Color', 'black');
plot(P, Bound,         '-',   'LineWidth', 2, 'Color', 'black');
ylim([0 15])
xlim([0 1])
set(gca, 'FontSize', 18)
xlabel('$$p$$', 'FontSize', 25, 'Interpreter', 'latex')
ylabel('$$g_1(p)$$', 'FontSize', 25, 'Interpreter', 'latex')
title('(c) Covariate selection: threshold vs. minimum, $h = 1$', ...
    'FontWeight', 'bold', 'FontSize', 20, 'Interpreter', 'latex')
legend({'Threshold', 'Minimum', 'No $p$-hacking', 'Bound'}, ...
    'Orientation', 'horizontal', 'Interpreter', 'latex', 'NumColumns', 2)
% Inset (zoom-in) for small p-values
axes('Position', [.35 .34 .45 .45])
box on
indexOfInterest = (P < 0.2) & (P > 0.0001);
plot(P(indexOfInterest), pcurves(indexOfInterest, 1), '--',  'LineWidth', 2, 'Color', 'black'); hold on
plot(P(indexOfInterest), pcurves(indexOfInterest, 2), ':',   'LineWidth', 2, 'Color', 'black');
plot(P(indexOfInterest), pcurves(indexOfInterest, 3), '-.',  'LineWidth', 2, 'Color', 'black');
plot(P(indexOfInterest), Bound(indexOfInterest), '-', 'LineWidth', 2, 'Color', 'black');
ylim([0 15])
xlim([0 0.1])
saveas(gcf, 'Figures/Figures_analytical_examples/one_sided/Fig1c', 'epsc')
close all