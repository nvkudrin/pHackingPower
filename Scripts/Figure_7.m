% ------------------------------------------------------------------------
% Figure 7: Data set selection p-curve, thresholding, one-sided tests
% ------------------------------------------------------------------------
clear all
mkdir('Figures/Figures_analytical_examples/one_sided')
rng(12345)

alpha = 0.05; % significance level
z = @(h, p) norminv(1 - p) - h;
H = [0, 1, 2]; %values of h

P = linspace(0.0001, 0.9999, 1000); %x-axis
pcurves4_t = zeros(length(P), 3); %for p-curves
Bound = exp(z(0, P).^2 / 2) .* (P <= 0.5) + (P > 0.5); %Upper Bound

% Calculate p-curves for different h
for d = 1:3
    h = H(d);
    for j = 1:length(P)
        p = P(j);
        Upsilon4_t = (1 + normcdf(z(h, alpha))) * (p <= alpha) + 2 * normcdf(z(h, p)) * (p > alpha);
        pcurves4_t(j, d) = exp(h * z(0, p) - h^2 / 2) * Upsilon4_t;
    end
end

% Figure 7 [Online Appendix]
figure('Visible','off')
plot(P, pcurves4_t(:, 1), '--',  'LineWidth', 2, 'color', 'black'); hold on
plot(P, pcurves4_t(:, 2), ':',   'LineWidth', 2, 'color', 'black');
plot(P, pcurves4_t(:, 3), '-.',  'LineWidth', 2, 'color', 'black');
plot(P, Bound,             '-',  'LineWidth', 2, 'color', 'black');
ylim([0 15])
set(gca, 'FontSize', 18)
xlabel('$$p$$', 'FontSize', 25, 'interpreter', 'latex')
ylabel('$$g_4^t(p)$$', 'FontSize', 25, 'interpreter', 'latex')
title('Dataset selection: threshold', ...
    'fontweight', 'bold', 'FontSize', 20, 'interpreter', 'latex')
legend('$h=0$', '$h=1$', '$h=2$', 'Bound', ...
    'Orientation', 'horizontal', 'interpreter', 'latex')
axes('position', [0.35, 0.34, 0.45, 0.45])
box on
indexOfInterest = (P < 0.2) & (P > 0.0001);
plot(P(indexOfInterest), pcurves4_t(indexOfInterest, 1), '--',  'LineWidth', 2, 'color', 'black'); hold on
plot(P(indexOfInterest), pcurves4_t(indexOfInterest, 2), ':',   'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest), pcurves4_t(indexOfInterest, 3), '-.',  'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest), Bound(indexOfInterest), '-', 'LineWidth', 2, 'color', 'black')
ylim([0 10])
xlim([0 0.1])
saveas(gcf, 'Figures/Figures_analytical_examples/one_sided/Fig7', 'epsc')
close all
