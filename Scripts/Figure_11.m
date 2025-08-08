% ------------------------------------------------------------------------
% Figure 11: Size distortion and relative bias under p-hacking, IV
% selection
% ------------------------------------------------------------------------
clear all
mkdir('Figures/Figures_analytical_examples/one_sided')
rng(12345)

z = @(h, p) norminv(1 - p) - h;
alpha = linspace(0.0001, 0.15, 1000); %significance levels

RR = zeros(length(alpha), 1); %Rejection Rates
%Calculate rejection rates
for j = 1:length(alpha)
    I = integral(@(x) normpdf(x) .* normcdf(sqrt(2) * z(0, alpha(j)) - x), ...
                 (sqrt(2) - 1) * z(0, alpha(j)), z(0, alpha(j)));
    RR(j) = 1 - normcdf(z(0, alpha(j))) * normcdf((sqrt(2) - 1) * z(0, alpha(j))) - I;
end


% === Figure 11 (a) [Online Appendix] ===================================
figure('Visible','off')
plot(alpha, RR, '--',  'LineWidth', 2, 'Color', 'black'); hold on
plot(alpha, alpha, '-', 'LineWidth', 2, 'Color', 'black');
ylim([0 0.3])
xlim([0 0.15])
set(gca, 'FontSize', 18)
xlabel('$$\alpha$$', 'FontSize', 25, 'Interpreter', 'latex')
ylabel('Rejection rate', 'FontSize', 25, 'Interpreter', 'latex')
title('(a) IV selection: size distortion', ...
    'FontWeight', 'bold', 'FontSize', 20, 'Interpreter', 'latex')
legend({'Under $p$-hacking', 'Nominal size'}, ...
    'Orientation', 'vertical', 'Interpreter', 'latex', 'location', 'northwest')
saveas(gcf, 'Figures/Figures_analytical_examples/one_sided/Fig11a', 'epsc')
close all

% === Figure 11 (b) [Online Appendix] ===================================
h = linspace(0, 2.5, 1000); %values of h
alpha = 0.05; %significance level
%Calculate biases
B2m = (1 / sqrt(2 - sqrt(2))) * normpdf(h * sqrt((sqrt(2) - 1) / sqrt(2))) .* (normcdf(h / sqrt(2 - sqrt(2)))) ...
    + sqrt(2) * normpdf(0) * (1 - normcdf(sqrt(2) * h));
B2t = B2m - (1 / sqrt(2 - sqrt(2))) * normpdf(h * sqrt((sqrt(2) - 1) / sqrt(2))) ...
    .* (normcdf(-sqrt(4 - 2 * sqrt(2)) * z(0, alpha) + h / sqrt(2 - sqrt(2))));

figure('Visible','off')
plot(h, B2t./h, '--', 'LineWidth', 2, 'Color', 'black'); hold on
plot(h, B2m./h, ':',  'LineWidth', 2, 'Color', 'black');
%ylim([0 0.6])
ylim([0 1])
yticks([0, 0.25, 0.5, 0.75, 1]);
yticklabels({'0%', '25%', '50%', '75%', '100%'});
xlim([0 2.5])
set(gca, 'FontSize', 18)
xlabel('$$h$$', 'FontSize', 25, 'Interpreter', 'latex')
ylabel('Relative Bias', 'FontSize', 25, 'Interpreter', 'latex')
title('(b) IV selection: relative bias', ...
    'FontWeight', 'bold', 'FontSize', 20, 'Interpreter', 'latex')
legend({'Threshold', 'Minimum'}, ...
    'Orientation', 'vertical', 'Interpreter', 'latex')
saveas(gcf, 'Figures/Figures_analytical_examples/one_sided/Fig11b', 'epsc')
close all