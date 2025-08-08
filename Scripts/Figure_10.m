% ------------------------------------------------------------------------
% Figure 10: Bias from covariate selection
% ------------------------------------------------------------------------
clear all
mkdir('Figures/Figures_analytical_examples/one_sided')
rng(12345)

% Figure 10 (a) [Online Appendix] - Bias

z = @(h, p) norminv(1 - p) - h;
alpha = 0.05; %significance level
h = linspace(0, 4, 1000); %values of h

% rho's implied by gamma = 0.1, 0.5, 0.9
rho1 = 1 - 0.1^2;
rho2 = 1 - 0.5^2;
rho3 = 1 - 0.9^2;

%Bias from minimum
Bias_min = [ ...
    normpdf(0) * sqrt(2) * 0.1 / sqrt(rho1) * ones(length(h), 1), ...
    normpdf(0) * sqrt(2) * 0.5 / sqrt(rho2) * ones(length(h), 1), ...
    normpdf(0) * sqrt(2) * 0.9 / sqrt(rho3) * ones(length(h), 1) ...
];

% Bias from thresholding
Bias_select = zeros(length(h), 3);
params = [rho1, rho2, rho3; 0.1, 0.5, 0.9];

for k = 1:3
    rho_k = params(1, k);
    gamma_k = params(2, k);
    Bias_select(:, k) = (1 / sqrt(rho_k)) * ( ...
        sqrt((1 - rho_k) / pi) * normcdf(sqrt(2 / (1 + rho_k)) * z(h, alpha)) + ...
        (1 - rho_k) * normpdf(z(h, alpha)) .* normcdf(-z(h, alpha) * sqrt((1 - rho_k) / (1 + rho_k))) ...
    );
end

% plot

figure('Visible','off')
plot(h, Bias_min(:, 1), '-',  'LineWidth', 2, 'Color', 'black'); hold on
plot(h, Bias_min(:, 2), '--', 'LineWidth', 2, 'Color', 'black');
plot(h, Bias_min(:, 3), ':',  'LineWidth', 2, 'Color', 'black');
plot(h, Bias_select(:, 1), '-o',  'MarkerIndices', 1:50:1000, 'LineWidth', 2, 'Color', 'black');
plot(h, Bias_select(:, 2), '--o', 'MarkerIndices', 1:50:1000, 'LineWidth', 2, 'Color', 'black');
plot(h, Bias_select(:, 3), ':o',  'MarkerIndices', 1:50:1000, 'LineWidth', 2, 'Color', 'black');
xlim([0 4])
set(gca, 'FontSize', 18)
xlabel('$$h$$', 'FontSize', 25, 'Interpreter', 'latex')
ylabel('$\sqrt{N}\times$Bias', 'FontSize', 25, 'Interpreter', 'latex')
title('(a) Covariate selection: bias', ...
    'FontWeight', 'bold', 'FontSize', 20, 'Interpreter', 'latex')
legend({'Minimum, $\gamma=0.1$', 'Minimum, $\gamma=0.5$', 'Minimum, $\gamma=0.9$', ...
        'Threshold, $\gamma=0.1$', 'Threshold, $\gamma=0.5$', 'Threshold, $\gamma=0.9$'}, ...
        'Orientation', 'vertical', 'Interpreter', 'latex', 'Location', 'west')
saveas(gcf, 'Figures/Figures_analytical_examples/one_sided/Fig10a', 'epsc')
close all

% Figure 10 (b) [Online Appendix] - Relative Bias

figure('Visible','off')
plot(h, sqrt(rho1)*Bias_min(:, 1)./h', '-',  'LineWidth', 2, 'Color', 'black'); hold on
plot(h, sqrt(rho2)*Bias_min(:, 2)./h', '--', 'LineWidth', 2, 'Color', 'black');
plot(h, sqrt(rho3)*Bias_min(:, 3)./h', ':',  'LineWidth', 2, 'Color', 'black');
plot(h, sqrt(rho1)*Bias_select(:, 1)./h', '-o',  'MarkerIndices', 1:50:1000, 'LineWidth', 2, 'Color', 'black');
plot(h, sqrt(rho2)*Bias_select(:, 2)./h', '--o', 'MarkerIndices', 1:50:1000, 'LineWidth', 2, 'Color', 'black');
plot(h, sqrt(rho3)*Bias_select(:, 3)./h', ':o',  'MarkerIndices', 1:50:1000, 'LineWidth', 2, 'Color', 'black');
xlim([0 4])
ylim([0 1])
yticks([0, 0.25, 0.5, 0.75, 1]);
yticklabels({'0%', '25%', '50%', '75%', '100%'});
set(gca, 'FontSize', 18)
xlabel('$$h$$', 'FontSize', 25, 'Interpreter', 'latex')
ylabel('Relative Bias', 'FontSize', 25, 'Interpreter', 'latex')
title('(b) Covariate selection: relative bias', ...
    'FontWeight', 'bold', 'FontSize', 20, 'Interpreter', 'latex')
legend({'Minimum, $\gamma=0.1$', 'Minimum, $\gamma=0.5$', 'Minimum, $\gamma=0.9$', ...
        'Threshold, $\gamma=0.1$', 'Threshold, $\gamma=0.5$', 'Threshold, $\gamma=0.9$'}, ...
        'Orientation', 'vertical', 'Interpreter', 'latex', 'Location', 'northeast')
saveas(gcf, 'Figures/Figures_analytical_examples/one_sided/Fig10b', 'epsc')
close all
