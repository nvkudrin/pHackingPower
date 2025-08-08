% ------------------------------------------------------------------------
% Figure 3: Lag Length Selection p-curves, one-sided tests
% ------------------------------------------------------------------------
clear all
mkdir('Figures/Figures_analytical_examples/one_sided')
rng(12345)

% Simulation setup for lag length selection
n = 200; %sample size per paper
M = 1e6; % number of MC draws for empirical distribution of sample autocorrelation
rho = zeros(M,1);
% calculate rho for each draw
for m = 1:M
    X = randn(n,1);
    rho(m) = sum((X(2:n) - mean(X)) .* (X(1:(n-1)) - mean(X))) / (n - 1);
end

alpha = 0.05;
kappa = 1/2;
L = -1/(2*kappa);
z = @(h, p) norminv(1 - p) - h;
omega = @(r) sqrt(max(1 + 2*kappa*r, 0));
H0 = mean(rho <= 0);
HL = mean(rho <= L);

P = linspace(0.0001, 0.9999, 1000); %x-axis
pcurve_t = zeros(length(P), 3); %thresholding p-curve
pcurve_m = zeros(length(P), 3); %minimum p-curve
Bound = exp(z(0, P).^2 / 2) .* (P <= 0.5) + (P > 0.5); %Upper Bound

for h = 0:2
    for j = 1:length(P)
        p = P(j);
        lp = L * (1 - (z(0, alpha) / z(0, p))^2);
        g = omega(rho) .* normpdf(z(0, p) * omega(rho) - h);

        % Thresholding approach
        if (p <= alpha)
            I = mean(g .* (rho <= lp) .* (rho > L));
            Upsilon3_t = 1 + I / normpdf(z(h, p));
        elseif (p > alpha) && (p <= 0.5)
            I = mean(g .* (rho <= 0) .* (rho > L));
            Upsilon3_t = 1 - H0 + HL + I / normpdf(z(h, p));
        else % (p > 0.5)
            I = mean(g .* (rho > 0));
            Upsilon3_t = H0 + I / normpdf(z(h, p));
            Upsilon3_m = Upsilon3_t;
        end

        % Minimum approach
        if (p <= 0.5)
            I = mean(g .* (rho <= 0) .* (rho > L));
            Upsilon3_m = 1 - H0 + HL + I / normpdf(z(h, p));
        end

        pcurve_t(j, h+1) = exp(h * z(0, p) - h^2 / 2) * Upsilon3_t;
        pcurve_m(j, h+1) = exp(h * z(0, p) - h^2 / 2) * Upsilon3_m;
    end
end

% ------------------------------------------------------------------------
% LAG LENGTH SELECTION: Figure 3 (a) [Main text] — Thresholding
% ------------------------------------------------------------------------
figure('Visible','off')
plot(P, pcurve_t(:,1), '--',  'LineWidth', 2, 'color', 'black'); hold on
plot(P, pcurve_t(:,2), ':',   'LineWidth', 2, 'color', 'black')
plot(P, pcurve_t(:,3), '-.',  'LineWidth', 2, 'color', 'black')
plot(P, Bound,         '-',   'LineWidth', 2, 'color', 'black')
ylim([0 15])
set(gca, 'FontSize', 18)
xlabel('$$p$$', 'FontSize', 25, 'interpreter', 'latex')
ylabel('$$g_3^t(p)$$', 'FontSize', 25, 'interpreter', 'latex')
title('(a) Lag length selection: threshold', ...
    'fontweight', 'bold', 'FontSize', 20, 'interpreter', 'latex')
legend('$h=0$', '$h=1$', '$h=2$', 'Bound', ...
    'Orientation', 'horizontal', 'interpreter', 'latex')
axes('position', [0.35, 0.34, 0.45, 0.45])
box on
indexOfInterest = (P < 0.2) & (P > 0.0001);
plot(P(indexOfInterest), pcurve_t(indexOfInterest,1), '--',  'LineWidth', 2, 'color', 'black'); hold on
plot(P(indexOfInterest), pcurve_t(indexOfInterest,2), ':',   'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest), pcurve_t(indexOfInterest,3), '-.',  'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest), Bound(indexOfInterest), '-', 'LineWidth', 2, 'color', 'black')
ylim([0 10])
xlim([0 0.1])
saveas(gcf, 'Figures/Figures_analytical_examples/one_sided/Fig3a', 'epsc')
close all

% ------------------------------------------------------------------------
% LAG LENGTH SELECTION: Figure 3 (b) [Main text] — Minimum
% ------------------------------------------------------------------------
figure('Visible','off')
plot(P, pcurve_m(:,1), '--',  'LineWidth', 2, 'color', 'black'); hold on
plot(P, pcurve_m(:,2), ':',   'LineWidth', 2, 'color', 'black')
plot(P, pcurve_m(:,3), '-.',  'LineWidth', 2, 'color', 'black')
plot(P, Bound,         '-',   'LineWidth', 2, 'color', 'black')
ylim([0 15])
set(gca, 'FontSize', 18)
xlabel('$$p$$', 'FontSize', 25, 'interpreter', 'latex')
ylabel('$$g_3^m(p)$$', 'FontSize', 25, 'interpreter', 'latex')
title('(b) Lag length selection: minimum', ...
    'fontweight', 'bold', 'FontSize', 20, 'interpreter', 'latex')
legend('$h=0$', '$h=1$', '$h=2$', 'Bound', ...
    'Orientation', 'horizontal', 'interpreter', 'latex')
axes('position', [0.35, 0.34, 0.45, 0.45])
box on
indexOfInterest = (P < 0.6) & (P > 0.0001);
plot(P(indexOfInterest), pcurve_m(indexOfInterest,1), '--',  'LineWidth', 2, 'color', 'black'); hold on
plot(P(indexOfInterest), pcurve_m(indexOfInterest,2), ':',   'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest), pcurve_m(indexOfInterest,3), '-.',  'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest), Bound(indexOfInterest), '-', 'LineWidth', 2, 'color', 'black')
ylim([0 2])
xlim([0 0.6])
saveas(gcf, 'Figures/Figures_analytical_examples/one_sided/Fig3b', 'epsc')
close all