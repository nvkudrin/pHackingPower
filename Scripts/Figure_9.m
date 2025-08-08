% ------------------------------------------------------------------------
% Figure 9: Size distortions under p-hacking, covariate selection
% ------------------------------------------------------------------------

clear all
mkdir('Figures/Figures_analytical_examples/one_sided')
rng(12345)

% === Figure 9 (a) [Online Appendix] ====================================

gamma = linspace(0.0001, 0.9999, 1000); %values of gamma
rho   = 1 - gamma.^2;
alpha = 0.05; %significance level
z = @(h, p) norminv(1 - p) - h;
RR = zeros(length(rho), 1); % rejection rates

% Calculate rejection rates
for j = 1:length(rho)
    RR(j) = 1 - mvncdf([z(0, alpha), z(0, alpha)], [0, 0], [1, rho(j); rho(j), 1]);
end

figure('Visible','off')
plot(gamma, RR, '-', 'LineWidth', 2, 'Color', 'black')
ylim([0.05 0.1])
xlim([0 1])
set(gca, 'FontSize', 18)
xlabel('$$\gamma$$', 'FontSize', 25, 'Interpreter', 'latex')
ylabel('Rejection rate', 'FontSize', 25, 'Interpreter', 'latex')
title('(a) Covariate selection: size distortion', ...
    'FontWeight', 'bold', 'FontSize', 20, 'Interpreter', 'latex')
saveas(gcf, 'Figures/Figures_analytical_examples/one_sided/Fig9a', 'epsc')
close all

% === Figure 9 (b) [Online Appendix] ====================================

alpha = linspace(0.0001, 0.15, 1000); %significance levels
gamma = [0.1, 0.5, 0.9]; % values of gamma
rho = 1 - gamma.^2;
RR = zeros(length(alpha), 3); % Rejection Rates

% Calculate rejection rates
for j = 1:length(alpha)
    for k = 1:3
        RR(j, k) = 1 - mvncdf([z(0, alpha(j)), z(0, alpha(j))], [0, 0], [1, rho(k); rho(k), 1]);
    end
end

figure('Visible','off')
plot(alpha, RR(:, 1), '--',  'LineWidth', 2, 'Color', 'black'); hold on
plot(alpha, RR(:, 2), ':',   'LineWidth', 2, 'Color', 'black');
plot(alpha, RR(:, 3), '-.',  'LineWidth', 2, 'Color', 'black');
plot(alpha, alpha,    '-',   'LineWidth', 2, 'Color', 'black');
ylim([0 0.3])
xlim([0 0.15])
set(gca, 'FontSize', 18)
xlabel('$$\alpha$$', 'FontSize', 25, 'Interpreter', 'latex')
ylabel('Rejection rate', 'FontSize', 25, 'Interpreter', 'latex')
title('(b) Covariate selection: size distortion', ...
    'FontWeight', 'bold', 'FontSize', 20, 'Interpreter', 'latex')
legend({'$\gamma=0.1$', '$\gamma=0.5$', '$\gamma=0.9$', 'Nominal size'}, ...
    'Orientation', 'vertical', 'Interpreter', 'latex', 'Location', 'northwest')
saveas(gcf, 'Figures/Figures_analytical_examples/one_sided/Fig9b', 'epsc')
close all