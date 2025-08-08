% ------------------------------------------------------------------------
% Figure 12: Size distortion under p-hacking from lag length selection
% ------------------------------------------------------------------------
clear all
mkdir('Figures/Figures_analytical_examples/one_sided')
rng(12345)

% Simulation setup for lag length selection
n = 200; % sample size per paper
M = 1e6; % number of MC draws for empirical distribution of sample autocorrelation
rho = zeros(M,1); 
% calculate autocorr per draw
for m = 1:M
    X = randn(n,1);
    rho(m) = sum((X(2:n) - mean(X)) .* (X(1:(n-1)) - mean(X))) / (n - 1);
end

kappa = 1/2;
L = -1/(2*kappa);
z = @(h, p) norminv(1 - p) - h;
omega = @(r) sqrt(max(1 + 2*kappa*r, 0));
H0 = mean(rho <= 0);
HL = mean(rho <= L);

%LAG LENGTH SELECTION: Rejection rate (Figure 12, Online Appendix)
alpha = linspace(0.0001, 0.15, 1000); % significance levels
RR = zeros(length(alpha), 1); %Rejection Rates
% Calculate rejection rates
for j = 1:length(alpha)
    g = normcdf(z(0, alpha(j)) * omega(rho));
    RR(j) = alpha(j) + (1 - alpha(j)) * (H0 - HL) - mean(g .* (rho <= 0) .* (rho > L));
end

figure('Visible','off')
plot(alpha, RR, '--',  'LineWidth', 2, 'color', 'black'); hold on
plot(alpha, alpha, '-', 'LineWidth', 2, 'color', 'black')
ylim([0 0.15])
xlim([0 0.15])
set(gca, 'FontSize', 18)
xlabel('$$\alpha$$', 'FontSize', 25, 'interpreter', 'latex')
ylabel('Rejection rate', 'FontSize', 25, 'interpreter', 'latex')
title('Lag length selection: size distortion', ...
    'fontweight', 'bold', 'FontSize', 20, 'interpreter', 'latex')
legend('Under $p$-hacking', 'Nominal size', ...
    'Orientation', 'vertical', 'interpreter', 'latex', 'location', 'northwest')
saveas(gcf, 'Figures/Figures_analytical_examples/one_sided/Fig12', 'epsc')
close all