% ------------------------------------------------------------------------
% Figure 8: Data set selection p-curve, thresholding, one-sided tests
% ------------------------------------------------------------------------

clear all
mkdir('Figures/Figures_analytical_examples/one_sided')
rng(12345)

z = @(h, p) norminv(1 - p) - h;
K = [2, 5, 20, 2, 5, 20]; % number of datasets
H = [0, 0, 0, 1, 1, 1]; % values of h

P = linspace(0.0001, 0.9999, 1000); %x-axis
pcurves4_m = zeros(length(P), 6); %for p-curves
Bound = exp(z(0, P).^2 / 2) .* (P <= 0.5) + (P > 0.5); %Upper Bound

%Calculate p-curves
for d = 1:6
    h = H(d);
    k = K(d);
    for j = 1:length(P)
        p = P(j);
        Upsilon4_m = k * normcdf(z(h, p))^(k - 1);
        pcurves4_m(j, d) = exp(h * z(0, p) - h^2 / 2) * Upsilon4_m;
    end
end

% Figure 8 (a) [Online Appendix]
figure('Visible','off')
plot(P, pcurves4_m(:, 1), '--',  'LineWidth', 2, 'color', 'black'); hold on
plot(P, pcurves4_m(:, 2), ':',   'LineWidth', 2, 'color', 'black')
plot(P, pcurves4_m(:, 3), '-.',  'LineWidth', 2, 'color', 'black')
plot(P, Bound,            '-',   'LineWidth', 2, 'color', 'black')
ylim([0 15])
set(gca, 'FontSize', 18)
xlabel('$$p$$', 'FontSize', 25, 'interpreter', 'latex')
ylabel('$$g_4^m(p;K)$$', 'FontSize', 25, 'interpreter', 'latex')
title('(a) Dataset selection: minimum, $h=0$', ...
    'fontweight', 'bold', 'FontSize', 20, 'interpreter', 'latex')
legend('$K=2$', '$K=5$', '$K=20$', 'Bound', ...
    'Orientation', 'horizontal', 'interpreter', 'latex')
axes('position', [0.35, 0.34, 0.45, 0.45])
box on
indexOfInterest = (P < 0.4) & (P > 0.0001);
plot(P(indexOfInterest), pcurves4_m(indexOfInterest, 1), '--',  'LineWidth', 2, 'color', 'black'); hold on
plot(P(indexOfInterest), pcurves4_m(indexOfInterest, 2), ':',   'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest), pcurves4_m(indexOfInterest, 3), '-.',  'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest), Bound(indexOfInterest), '-', 'LineWidth', 2, 'color', 'black')
ylim([0 20])
xlim([0 0.25])
saveas(gcf, 'Figures/Figures_analytical_examples/one_sided/Fig8a', 'epsc')
close all

% Figure 8 (b) [Online Appendix]
figure('Visible','off')
plot(P, pcurves4_m(:, 4), '--',  'LineWidth', 2, 'color', 'black'); hold on
plot(P, pcurves4_m(:, 5), ':',   'LineWidth', 2, 'color', 'black')
plot(P, pcurves4_m(:, 6), '-.',  'LineWidth', 2, 'color', 'black')
plot(P, Bound,            '-',   'LineWidth', 2, 'color', 'black')
ylim([0 15])
set(gca, 'FontSize', 18)
xlabel('$$p$$', 'FontSize', 25, 'interpreter', 'latex')
ylabel('$$g_4^m(p;K)$$', 'FontSize', 25, 'interpreter', 'latex')
title('(b) Dataset selection: minimum, $h=1$', ...
    'fontweight', 'bold', 'FontSize', 20, 'interpreter', 'latex')
legend('$K=2$', '$K=5$', '$K=20$', 'Bound', ...
    'Orientation', 'horizontal', 'interpreter', 'latex')
axes('position', [0.35, 0.34, 0.45, 0.45])
box on
indexOfInterest = (P < 0.4) & (P > 0.0001);
plot(P(indexOfInterest), pcurves4_m(indexOfInterest, 4), '--',  'LineWidth', 2, 'color', 'black'); hold on
plot(P(indexOfInterest), pcurves4_m(indexOfInterest, 5), ':',   'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest), pcurves4_m(indexOfInterest, 6), '-.',  'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest), Bound(indexOfInterest), '-', 'LineWidth', 2, 'color', 'black')
ylim([0 20])
xlim([0 0.25])
saveas(gcf, 'Figures/Figures_analytical_examples/one_sided/Fig8b', 'epsc')
close all