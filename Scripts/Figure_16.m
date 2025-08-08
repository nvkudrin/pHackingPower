% ------------------------------------------------------------------------
% Figure 16: p-Curves from dataset selection based on the minimum approach with two-sided tests
% ------------------------------------------------------------------------
clear all
mkdir('Figures/Figures_analytical_examples/two_sided')
rng(12345)
% ---------------------- Minimum: Calculate Curves ------------------------
z = @(h,p) norminv(1 - p) - h;
K = [2, 5, 20, 2, 5, 20]; % number of datasets
H = [0, 0, 0, 1, 1, 1]; % values of h
P = linspace(0.0001, 0.9999, 1000); %x-axis
pcurves4_m = zeros(length(P), 6); % for p-curves

% Recompute Bound for each P value using two-sided formula
Bound = zeros(size(P));
h_grid = linspace(0, 6, 10000); % grid of values of h
cv2 = @(p) norminv(1 - p/2);
for j = 1:length(P)
    Bound(j) = max((exp(h_grid*cv2(P(j)) - h_grid.^2/2) + ...
                    exp(-h_grid*cv2(P(j)) - h_grid.^2/2))/2);
end

% Calculate p-curves
for d = 1:6
    h = H(d);
    k = K(d);
    for j = 1:length(P)
        p = P(j);
        Upsilon_4_m = k * (normcdf(z(-h,p/2)) - normcdf(-z(h,p/2)))^(k-1);
        pcurves4_m(j,d) = 0.5 * (exp(h*z(0,p/2) - h^2/2) + exp(-h*z(0,p/2) - h^2/2)) * Upsilon_4_m;
    end
end

% -------------------------- Figure 16 (a) --------------------------------
close all
figure('Visible', 'off')
plot(P, pcurves4_m(:,1), '--' , 'LineWidth', 2, 'color', 'black')
hold on
plot(P, pcurves4_m(:,2), ':' , 'LineWidth', 2, 'color', 'black')
plot(P, pcurves4_m(:,3), '-.' , 'LineWidth', 2, 'color', 'black')
plot(P, Bound, '-' , 'LineWidth', 2, 'color', 'black')
ylim([0 15])
set(gca,'FontSize',18)
xlabel('$$p$$',   'FontSize',25, 'interpreter', 'latex')
ylabel('$$g_4^m(p;K)$$', 'FontSize',25, 'interpreter', 'latex')
title('(a) Dataset selection: minimum, $h=0$', 'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
legend('$K=2$', '$K=5$', '$K=20$', 'Bound', 'Orientation','horizontal', 'interpreter', 'latex')

axes('position',[.35 .34 .45 .45])
box on 
indexOfInterest = (P < 0.4) & (P > 0.0001);
plot(P(indexOfInterest), pcurves4_m(indexOfInterest,1), '--' , 'LineWidth', 2, 'color', 'black')
hold on
plot(P(indexOfInterest), pcurves4_m(indexOfInterest,2), ':' , 'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest), pcurves4_m(indexOfInterest,3), '-.' , 'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest), Bound(indexOfInterest), '-' , 'LineWidth', 2, 'color', 'black')
ylim([0 10])
xlim([0 0.25])
saveas(gcf,'Figures/Figures_analytical_examples/two_sided/Fig16a', 'epsc')
close all

% -------------------------- Figure 16 (b) --------------------------------
close all
figure('Visible', 'off')
plot(P, pcurves4_m(:,4), '--' , 'LineWidth', 2, 'color', 'black')
hold on
plot(P, pcurves4_m(:,5), ':' , 'LineWidth', 2, 'color', 'black')
plot(P, pcurves4_m(:,6), '-.' , 'LineWidth', 2, 'color', 'black')
plot(P, Bound, '-' , 'LineWidth', 2, 'color', 'black')
ylim([0 15])
set(gca,'FontSize',18)
xlabel('$$p$$',   'FontSize',25, 'interpreter', 'latex')
ylabel('$$g_4^m(p;K)$$', 'FontSize',25, 'interpreter', 'latex')
title('(b) Dataset selection: minimum, $h=1$', 'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
legend('$K=2$', '$K=5$', '$K=20$', 'Bound', 'Orientation','horizontal', 'interpreter', 'latex')

axes('position',[.35 .34 .45 .45])
box on 
indexOfInterest = (P < 0.4) & (P > 0.0001);
plot(P(indexOfInterest), pcurves4_m(indexOfInterest,4), '--' , 'LineWidth', 2, 'color', 'black')
hold on
plot(P(indexOfInterest), pcurves4_m(indexOfInterest,5), ':' , 'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest), pcurves4_m(indexOfInterest,6), '-.' , 'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest), Bound(indexOfInterest), '-' , 'LineWidth', 2, 'color', 'black')
ylim([0 10])
xlim([0 0.25])
saveas(gcf,'Figures/Figures_analytical_examples/two_sided/Fig16b', 'epsc')
close all