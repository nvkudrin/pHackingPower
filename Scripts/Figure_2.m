% ------------------------------------------------------------------------
% Figure 2: IV Selection p-curves, one-sided tests
% ------------------------------------------------------------------------
clear all
mkdir('Figures/Figures_analytical_examples/one_sided')
rng(12345)

alpha = 0.05; %significance level
z = @(h, p) norminv(1 - p) - h;
zeta = @(p) 1 - 2 * normcdf((1 - sqrt(2)) * z(0, p));
D = @(h, p) sqrt(2) * z(0, p) - 2 * h;
H = [0, 1, 2]; %values of h for (a) and (b)

P = linspace(0.0001, 0.9999, 1000); %x-axis
pcurves2_t = zeros(length(P), 3); %thresholding p-curves
pcurves2_m = zeros(length(P), 3); %minimum p-curves
Bound = exp(z(0, P).^2 / 2) .* (P <= 0.5) + (P > 0.5); % Upper Bound
No_phacking = exp(1 * z(0, P) - 1^2 / 2); % for h = 1

for d = 1:3
    h = H(d);
    for j = 1:length(P)
        p = P(j);

        if (p <= alpha)
            Upsilon_2_t = normpdf(z(sqrt(2) * h, p)) / normpdf(z(h, p)) + 2 * normcdf(D(h, alpha) - z(h, p));
            Upsilon_2_m = zeta(p) * normpdf(z(sqrt(2) * h, p)) / normpdf(z(h, p)) + 2 * normcdf(D(h, p) - z(h, p));
        elseif (p > alpha) && (p <= 0.5)
            Upsilon_2_t = zeta(p) * normpdf(z(sqrt(2) * h, p)) / normpdf(z(h, p)) + 2 * normcdf(D(h, p) - z(h, p));
            Upsilon_2_m = Upsilon_2_t;
        else % p > 0.5
            Upsilon_2_t = 2 * normcdf(z(h, p));
            Upsilon_2_m = Upsilon_2_t;
        end

        pcurves2_t(j, d) = exp(h * z(0, p) - h^2 / 2) * Upsilon_2_t;
        pcurves2_m(j, d) = exp(h * z(0, p) - h^2 / 2) * Upsilon_2_m;
    end
end

% === Figure 2 (a) [Main text]: Thresholding IV Selection ===============
figure('Visible','off')
plot(P, pcurves2_t(:, 1), '--',  'LineWidth', 2, 'Color', 'black'); hold on
plot(P, pcurves2_t(:, 2), ':',   'LineWidth', 2, 'Color', 'black');
plot(P, pcurves2_t(:, 3), '-.',  'LineWidth', 2, 'Color', 'black');
plot(P, Bound,            '-',   'LineWidth', 2, 'Color', 'black');
ylim([0 15])
set(gca, 'FontSize', 18)
xlabel('$$p$$', 'FontSize', 25, 'Interpreter', 'latex')
ylabel('$$g_2^t(p)$$', 'FontSize', 25, 'Interpreter', 'latex')
title('(a) IV selection: threshold', ...
    'FontWeight', 'bold', 'FontSize', 20, 'Interpreter', 'latex')
legend({'$h=0$', '$h=1$', '$h=2$', 'Bound'}, ...
    'Orientation', 'horizontal', 'Interpreter', 'latex')
axes('position', [.35 .34 .45 .45])
box on
indexOfInterest = (P < 0.2) & (P > 0.0001);
plot(P(indexOfInterest), pcurves2_t(indexOfInterest, 1), '--',  'LineWidth', 2, 'Color', 'black'); hold on
plot(P(indexOfInterest), pcurves2_t(indexOfInterest, 2), ':',   'LineWidth', 2, 'Color', 'black');
plot(P(indexOfInterest), pcurves2_t(indexOfInterest, 3), '-.',  'LineWidth', 2, 'Color', 'black');
plot(P(indexOfInterest), Bound(indexOfInterest), '-', 'LineWidth', 2, 'Color', 'black');
ylim([0 10])
xlim([0 0.1])
saveas(gcf, 'Figures/Figures_analytical_examples/one_sided/Fig2a', 'epsc')
close all

% === Figure 2 (b) [Main text]: Minimum IV Selection ===================
figure('Visible','off')
plot(P, pcurves2_m(:, 1), '--',  'LineWidth', 2, 'Color', 'black'); hold on
plot(P, pcurves2_m(:, 2), ':',   'LineWidth', 2, 'Color', 'black');
plot(P, pcurves2_m(:, 3), '-.',  'LineWidth', 2, 'Color', 'black');
plot(P, Bound,            '-',   'LineWidth', 2, 'Color', 'black');
ylim([0 15])
set(gca, 'FontSize', 18)
xlabel('$$p$$', 'FontSize', 25, 'Interpreter', 'latex')
ylabel('$$g_2^m(p)$$', 'FontSize', 25, 'Interpreter', 'latex')
title('(b) IV selection: minimum', ...
    'FontWeight', 'bold', 'FontSize', 20, 'Interpreter', 'latex')
legend({'$h=0$', '$h=1$', '$h=2$', 'Bound'}, ...
    'Orientation', 'horizontal', 'Interpreter', 'latex')
axes('position', [.35 .34 .45 .45])
box on
indexOfInterest = (P < 0.6) & (P > 0.0001);
plot(P(indexOfInterest), pcurves2_m(indexOfInterest, 1), '--',  'LineWidth', 2, 'Color', 'black'); hold on
plot(P(indexOfInterest), pcurves2_m(indexOfInterest, 2), ':',   'LineWidth', 2, 'Color', 'black');
plot(P(indexOfInterest), pcurves2_m(indexOfInterest, 3), '-.',  'LineWidth', 2, 'Color', 'black');
plot(P(indexOfInterest), Bound(indexOfInterest), '-', 'LineWidth', 2, 'Color', 'black');
ylim([0 10])
xlim([0 0.1])
saveas(gcf, 'Figures/Figures_analytical_examples/one_sided/Fig2b', 'epsc')
close all

% === Figure 2 (c) [Main text]: Threshold vs. Minimum ===================
figure('Visible','off')
plot(P, pcurves2_t(:, 2), '--',  'LineWidth', 2, 'Color', 'black'); hold on
plot(P, pcurves2_m(:, 2), ':',   'LineWidth', 2, 'Color', 'black');
plot(P, No_phacking,      '-.',  'LineWidth', 2, 'Color', 'black');
plot(P, Bound,            '-',   'LineWidth', 2, 'Color', 'black');
ylim([0 15])
xlim([0 1])
set(gca, 'FontSize', 18)
xlabel('$$p$$', 'FontSize', 25, 'Interpreter', 'latex')
ylabel('$$g_2(p)$$', 'FontSize', 25, 'Interpreter', 'latex')
title('(c) IV selection: threshold vs. minimum, $h = 1$', ...
    'FontWeight', 'bold', 'FontSize', 20, 'Interpreter', 'latex')
legend({'Threshold', 'Minimum', 'No $p$-hacking', 'Bound'}, ...
    'Orientation', 'horizontal', 'Interpreter', 'latex', 'NumColumns', 2)
axes('position', [.35 .34 .45 .45])
box on
indexOfInterest = (P < 0.6) & (P > 0.0001);
plot(P(indexOfInterest), pcurves2_t(indexOfInterest, 2), '--',  'LineWidth', 2, 'Color', 'black'); hold on
plot(P(indexOfInterest), pcurves2_m(indexOfInterest, 2), ':',   'LineWidth', 2, 'Color', 'black');
plot(P(indexOfInterest), No_phacking(indexOfInterest), '-.',  'LineWidth', 2, 'Color', 'black');
plot(P(indexOfInterest), Bound(indexOfInterest), '-', 'LineWidth', 2, 'Color', 'black');
ylim([0 15])
xlim([0 0.1])
saveas(gcf, 'Figures/Figures_analytical_examples/one_sided/Fig2c', 'epsc')
close all