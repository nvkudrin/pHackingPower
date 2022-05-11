% This code replicates Figures 1-7 from "(When) Can we detect p-hacking"
%Authors: G. Elliott, N. Kudrin, K. Wuthrich
%%
close all
clear all
rng(12345)
%% Covariate selection Example
alpha = 0.05;
z = @(h,p) norminv(1 - p) - h;
gamma = [0.5, 0.5, 0.5, 0.1, 0.5, 0.9];
H = [0, 1, 2, 1, 1, 1];

P = linspace(0.0001, 0.9999, 1000);
pcurves = zeros(length(P), 6);
Bound = exp(z(0,P).^2/2).*(P<=0.5) + (P>0.5);
for d = 1:6
    rho = 1 - gamma(d)^2;
    h = H(d);
for j = 1:length(P)
    p = P(j);

Upsilon_covsel_t = (1+normcdf((z(h,alpha) - rho*z(h,p))/sqrt(1-rho^2)))*(p<=alpha)+...
    2*normcdf(z(h,p)*sqrt((1-rho)/(1+rho)))*(p>alpha);

pcurves(j,d) = exp(h*z(0,p) - h^2/2)*Upsilon_covsel_t;

end
end

%Figure 1 Left
figure(1)
plot(P, pcurves(:,1), '--' , 'LineWidth', 2, 'color', 'black')
hold on
plot(P, pcurves(:,2), ':' , 'LineWidth', 2, 'color', 'black')
plot(P, pcurves(:,3), '-.' , 'LineWidth', 2, 'color', 'black')
plot(P, Bound, '-' , 'LineWidth', 2, 'color', 'black')
ylim([0 15])
set(gca,'FontSize',18)
xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('$$g_1^t(p)$$', 'FontSize',25, 'interpreter', 'latex')
title('Covariate selection: threshold, $\gamma=0.5$','fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
lgd = legend('$h=0$', '$h=1$','$h=2$','Bound', 'Orientation','horizontal', 'interpreter', 'latex');
axes('position',[.35 .34 .45 .45])
box on 
indexOfInterest = (P < 0.2) & (P > 0.0001); 
plot(P(indexOfInterest),pcurves(indexOfInterest,1), '--' , 'LineWidth', 2, 'color', 'black') % plot on new axes
hold on
plot(P(indexOfInterest),pcurves(indexOfInterest,2), ':' , 'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest),pcurves(indexOfInterest,3), '-.' , 'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest),Bound(indexOfInterest), '-' , 'LineWidth', 2, 'color', 'black')
ylim([0 12])
xlim([0 0.1])
saveas(gcf,'BW/CovarChoice1_bw', 'epsc')
close all

% Figure 1 Right
figure(2)
plot(P, pcurves(:,4), '--' , 'LineWidth', 2, 'color', 'black')
hold on
plot(P, pcurves(:,5), ':' , 'LineWidth', 2, 'color', 'black')
plot(P, pcurves(:,6), '-.' , 'LineWidth', 2, 'color', 'black')
plot(P, Bound, '-' , 'LineWidth', 2, 'color', 'black')
ylim([0 15])
set(gca,'FontSize',18)
xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('$$g_1^t(p)$$', 'FontSize',25, 'interpreter', 'latex')
title('Covariate selection: threshold, $h=1$','fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
lgd = legend('$\gamma=0.1$', '$\gamma=0.5$','$\gamma=0.9$','Bound', 'Orientation','horizontal', 'interpreter', 'latex');

axes('position',[.35 .34 .45 .45])
box on 
indexOfInterest = (P < 0.6) & (P > 0.0001); 
plot(P(indexOfInterest),pcurves(indexOfInterest,4), '--' , 'LineWidth', 2, 'color', 'black') 
hold on
plot(P(indexOfInterest),pcurves(indexOfInterest,5), ':' , 'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest),pcurves(indexOfInterest,6), '-.' , 'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest),Bound(indexOfInterest), '-' , 'LineWidth', 2, 'color', 'black')
ylim([0 15])
xlim([0 0.1])
saveas(gcf,'BW/CovarChoice2_bw', 'epsc')
close all


h = 1;
gamma = 0.5;
rho = 1 - gamma^2;
pcurves = zeros(length(P), 4);
Bound = exp(z(0,P).^2/2).*(P<=0.5) + (P>0.5);
for j = 1:length(P)
    p = P(j);

Upsilon_covsel_t = (1+normcdf((z(h,alpha) - rho*z(h,p))/sqrt(1-rho^2)))*(p<=alpha)+...
    2*normcdf(z(h,p)*sqrt((1-rho)/(1+rho)))*(p>alpha);
Upsilon_covsel_m = 2*normcdf(z(h,p)*sqrt((1-rho)/(1+rho)));

pcurves(j,1) = exp(h*z(0,p) - h^2/2)*Upsilon_covsel_t;
pcurves(j,2) = exp(h*z(0,p) - h^2/2)*Upsilon_covsel_m;
pcurves(j,3) = exp(h*z(0,p) - h^2/2);
end

% Figure 2
figure(3)
plot(P, pcurves(:,1), '--' , 'LineWidth', 2, 'color', 'black')
hold on
plot(P, pcurves(:,2), ':' , 'LineWidth', 2, 'color', 'black')
plot(P, pcurves(:,3), '-.' , 'LineWidth', 2, 'color', 'black')
plot(P, Bound, '-' , 'LineWidth', 2, 'color', 'black')
ylim([0 15])
xlim([0 1])
set(gca,'FontSize',18)
xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('$$g_1(p)$$', 'FontSize',25, 'interpreter', 'latex')
title('Covariate selection: threshold vs. minimum, $h = 1$','fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
lgd = legend('Threshold', 'Minimum','No $p$-hacking','Bound', 'Orientation','horizontal', 'interpreter', 'latex', 'NumColumns',2);

axes('position',[.35 .34 .45 .45])
box on 
indexOfInterest = (P < 0.2) & (P > 0.0001); 
plot(P(indexOfInterest),pcurves(indexOfInterest,1), '--' , 'LineWidth', 2, 'color', 'black') % plot on new axes
hold on
plot(P(indexOfInterest),pcurves(indexOfInterest,2), ':' , 'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest),pcurves(indexOfInterest,3), '-.' , 'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest),Bound(indexOfInterest), '-' , 'LineWidth', 2, 'color', 'black')
ylim([0 15])
xlim([0 0.1])

saveas(gcf,'BW/CovarChoice3_bw', 'epsc')
close all
%% Covariate selection Size and Bias
clear all
gamma = linspace(0.0001, 0.9999, 1000);
rho = 1-gamma.^2;
alpha=0.05;
z = @(h,p) norminv(1-p)-h;
RR = zeros(length(rho), 1);
for j = 1:length(rho)
    RR(j) = 1 - mvncdf([z(0,alpha), z(0, alpha)], [0, 0], [1 rho(j); rho(j) 1]);
end

%Figure 3 Left
figure(4)
plot(gamma, RR, '-' , 'LineWidth', 2, 'color', 'black')
ylim([0.05 0.1])
xlim([0 1])
set(gca,'FontSize',18)
xlabel('$$\gamma$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('Rejection rate', 'FontSize',25, 'interpreter', 'latex')
title('Covariate selection: size distortion','fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
saveas(gcf,'BW/CovarSize1_bw', 'epsc')
close all
%%
alpha = linspace(0.0001, 0.15, 1000);
gamma = [0.1, 0.5, 0.9];
rho = 1 - gamma.^2;
RR = zeros(length(alpha), 3);

for j = 1:length(alpha)
    RR(j, 1) = 1 - mvncdf([z(0,alpha(j)), z(0, alpha(j))], [0, 0], [1 rho(1); rho(1) 1]);
    RR(j, 2) = 1 - mvncdf([z(0,alpha(j)), z(0, alpha(j))], [0, 0], [1 rho(2); rho(2) 1]);
    RR(j, 3) = 1 - mvncdf([z(0,alpha(j)), z(0, alpha(j))], [0, 0], [1 rho(3); rho(3) 1]);
end
% Figure 3 Right
figure(5)
plot(alpha, RR(:,1), '--' , 'LineWidth', 2, 'color', 'black')
hold on
plot(alpha, RR(:,2), ':' , 'LineWidth', 2, 'color', 'black')
plot(alpha, RR(:,3), '-.' , 'LineWidth', 2, 'color', 'black')
plot(alpha, alpha, '-' , 'LineWidth', 2, 'color', 'black')
ylim([0 0.3])
xlim([0 0.15])
set(gca,'FontSize',18)
xlabel('$$\alpha$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('Rejection rate', 'FontSize',25, 'interpreter', 'latex')
title('Covariate selection: size distortion','fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
lgd = legend('$\gamma=0.1$', '$\gamma=0.5$','$\gamma=0.9$','Nominal size', 'Orientation','vertical', 'interpreter', 'latex', 'location', 'northwest');
saveas(gcf,'BW/CovarSize2_bw', 'epsc')
close all

%% Figure 4
clear all
z = @(h,p) norminv(1-p)-h;
alpha = 0.05;
h = linspace(0, 4, 1000);
rho1 = 1 - 0.1^2;
rho2 = 1 - 0.5^2;
rho3 = 1 - 0.9^2;
Bias_min = zeros(length(h), 3);
Bias_min(:,1) = normpdf(0)*sqrt(2)*ones(length(h), 1)*0.1/sqrt(rho1);
Bias_min(:,2) = normpdf(0)*sqrt(2)*ones(length(h), 1)*0.5/sqrt(rho2);
Bias_min(:,3) = normpdf(0)*sqrt(2)*ones(length(h), 1)*0.9/sqrt(rho3);

Bias_select = zeros(length(h), 3);
Bias_select(:,1) = (1/sqrt(rho1))*(sqrt((1-rho1)/pi)*normcdf(sqrt(2/(1+rho1))*z(h,alpha))+...
    (1-rho1)*normpdf(z(h,alpha)).*normcdf(-z(h,alpha)*sqrt((1-rho1)/(1+rho1))));
Bias_select(:,2) = (1/sqrt(rho2))*(sqrt((1-rho2)/pi)*normcdf(sqrt(2/(1+rho2))*z(h,alpha))+...
    (1-rho2)*normpdf(z(h,alpha)).*normcdf(-z(h,alpha)*sqrt((1-rho2)/(1+rho2))));
Bias_select(:,3) = (1/sqrt(rho3))*(sqrt((1-rho3)/pi)*normcdf(sqrt(2/(1+rho3))*z(h,alpha))+...
    (1-rho3)*normpdf(z(h,alpha)).*normcdf(-z(h,alpha)*sqrt((1-rho3)/(1+rho3))));
Bias_select = Bias_select;

figure(6)
plot(h, Bias_min(:,1), '-' , 'LineWidth', 2, 'color', 'black')
hold on
plot(h, Bias_min(:,2), '--' , 'LineWidth', 2, 'color', 'black')
plot(h, Bias_min(:,3), ':' , 'LineWidth', 2, 'color', 'black')

plot(h, Bias_select(:,1), '-o' ,'MarkerIndices', [1:50:1000], 'LineWidth', 2, 'color', 'black')
hold on
plot(h, Bias_select(:,2), '--o' ,'MarkerIndices', [1:50:1000], 'LineWidth', 2, 'color', 'black')
plot(h, Bias_select(:,3), ':o' ,'MarkerIndices', [1:50:1000], 'LineWidth', 2, 'color', 'black')

xlim([0 4])
set(gca,'FontSize',18)
xlabel('$$h$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('$\sqrt{N}\times$Bias', 'FontSize',25, 'interpreter', 'latex')
title('Covariate selection: bias','fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
lgd = legend('Minimum, $\gamma=0.1$', 'Minimum, $\gamma=0.5$','Minimum, $\gamma=0.9$','Threshold, $\gamma=0.1$', 'Threshold, $\gamma=0.5$','Threshold, $\gamma=0.9$', 'Orientation','vertical', 'interpreter', 'latex', 'location', 'west');
saveas(gcf,'BW/CovarBias_bw', 'epsc') 
close all

%%
%%%%%%%%%%%%%%%% IV Example %%%%%%%%%%%%%%%%%%
clear all
alpha = 0.05;
z = @(h,p) norminv(1 - p) - h;
zeta = @(p) 1 - 2*normcdf((1 - sqrt(2))*z(0,p));
D = @(h,p) sqrt(2)*z(0,p) - 2*h;

H = [0, 1, 2];

P = linspace(0.0001, 0.9999, 1000);
pcurves2_t = zeros(length(P), 3);
pcurves2_m = zeros(length(P), 3);
Bound = exp(z(0,P).^2/2).*(P<=0.5) + (P>0.5);
No_phacking = exp(1*z(0,P) - 1^2/2); %for h = 1
for d = 1:3
    h = H(d);
for j = 1:length(P)
    p = P(j);

    if (p<=alpha)
Upsilon_2_t = normpdf(z(sqrt(2)*h, p))/normpdf(z(h,p)) + 2*normcdf(D(h,alpha) - z(h,p));
Upsilon_2_m = zeta(p)*normpdf(z(sqrt(2)*h, p))/normpdf(z(h,p)) + 2*normcdf(D(h,p) - z(h,p));
    end
    if (p>alpha)&&(p<=0.5)
Upsilon_2_t = zeta(p)*normpdf(z(sqrt(2)*h, p))/normpdf(z(h,p)) + 2*normcdf(D(h,p) - z(h,p));
Upsilon_2_m = zeta(p)*normpdf(z(sqrt(2)*h, p))/normpdf(z(h,p)) + 2*normcdf(D(h,p) - z(h,p));
    end
    if (p>0.5)
        Upsilon_2_t = 2*normcdf(z(h,p));
        Upsilon_2_m = 2*normcdf(z(h,p));
    end

pcurves2_t(j,d) = exp(h*z(0,p) - h^2/2)*Upsilon_2_t;
pcurves2_m(j,d) = exp(h*z(0,p) - h^2/2)*Upsilon_2_m;
end
end

%Figure 5 Left
figure(7)
plot(P, pcurves2_t(:,1), '--' , 'LineWidth', 2, 'color', 'black')
hold on
plot(P, pcurves2_t(:,2), ':' , 'LineWidth', 2, 'color', 'black')
plot(P, pcurves2_t(:,3), '-.' , 'LineWidth', 2, 'color', 'black')
plot(P, Bound, '-' , 'LineWidth', 2, 'color', 'black')
ylim([0 15])
set(gca,'FontSize',18)
xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('$$g_2^t(p)$$', 'FontSize',25, 'interpreter', 'latex')
title('IV selection: threshold','fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
lgd = legend('$h=0$', '$h=1$','$h=2$','Bound', 'Orientation','horizontal', 'interpreter', 'latex');

axes('position',[.35 .34 .45 .45])
box on 
indexOfInterest = (P < 0.2) & (P > 0.0001);
plot(P(indexOfInterest),pcurves2_t(indexOfInterest,1), '--' , 'LineWidth', 2, 'color', 'black') % plot on new axes
hold on
plot(P(indexOfInterest),pcurves2_t(indexOfInterest,2), ':' , 'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest),pcurves2_t(indexOfInterest,3), '-.' , 'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest),Bound(indexOfInterest), '-' , 'LineWidth', 2, 'color', 'black')
ylim([0 10])
xlim([0 0.1])
saveas(gcf,'BW/IVChoice1_bw', 'epsc')
close all

% Figure 5 Right
figure(8)
plot(P, pcurves2_m(:,1), '--' , 'LineWidth', 2, 'color', 'black')
hold on
plot(P, pcurves2_m(:,2), ':' , 'LineWidth', 2, 'color', 'black')
plot(P, pcurves2_m(:,3), '-.' , 'LineWidth', 2, 'color', 'black')
plot(P, Bound, '-' , 'LineWidth', 2, 'color', 'black')
ylim([0 15])
set(gca,'FontSize',18)
xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('$$g_2^m(p)$$', 'FontSize',25, 'interpreter', 'latex')
title('IV selection: minimum','fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
lgd = legend('$h=0$', '$h=1$','$h=2$','Bound', 'Orientation','horizontal', 'interpreter', 'latex');

axes('position',[.35 .34 .45 .45])
box on 
indexOfInterest = (P < 0.6) & (P > 0.0001);
plot(P(indexOfInterest),pcurves2_m(indexOfInterest,1), '--' , 'LineWidth', 2, 'color', 'black') % plot on new axes
hold on
plot(P(indexOfInterest),pcurves2_m(indexOfInterest,2), ':' , 'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest),pcurves2_m(indexOfInterest,3), '-.' , 'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest),Bound(indexOfInterest), '-' , 'LineWidth', 2, 'color', 'black')
ylim([0 10])
xlim([0 0.1])
saveas(gcf,'BW/IVChoice2_bw', 'epsc')
close all

%Figure 6
figure(9)
plot(P, pcurves2_t(:,2), '--' , 'LineWidth', 2, 'color', 'black')
hold on
plot(P, pcurves2_m(:,2), ':' , 'LineWidth', 2, 'color', 'black')
plot(P, No_phacking, '-.' , 'LineWidth', 2, 'color', 'black')
plot(P, Bound, '-' , 'LineWidth', 2, 'color', 'black')
ylim([0 15])
xlim([0 1])
set(gca,'FontSize',18)
xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('$$g_2(p)$$', 'FontSize',25, 'interpreter', 'latex')
title('IV selection: threshold vs. minimum, $h = 1$','fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
lgd = legend('Threshold', 'Minimum','No $p$-hacking','Bound', 'Orientation','horizontal', 'interpreter', 'latex', 'NumColumns',2);

axes('position',[.35 .34 .45 .45])
box on 
indexOfInterest = (P < 0.6) & (P > 0.0001); 
plot(P(indexOfInterest),pcurves2_t(indexOfInterest,2), '--' , 'LineWidth', 2, 'color', 'black') % plot on new axes
hold on
plot(P(indexOfInterest),pcurves2_m(indexOfInterest,2), ':' , 'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest),No_phacking(indexOfInterest), '-.' , 'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest),Bound(indexOfInterest), '-' , 'LineWidth', 2, 'color', 'black')
ylim([0 15])
xlim([0 0.1])

saveas(gcf,'BW/IVChoice3_bw', 'epsc')
close all

%% IV size and Bias
clear all
z = @(h,p) norminv(1-p) - h;

alpha = linspace(0.0001, 0.15, 1000);

RR = zeros(length(alpha), 1);
for j = 1:length(alpha)
    I = integral(@(x) normpdf(x).*normcdf(sqrt(2)*z(0,alpha(j))-x), (sqrt(2)-1)*z(0,alpha(j)), z(0,alpha(j)));
    RR(j) = 1 - normcdf(z(0,alpha(j)))*normcdf((sqrt(2)-1)*z(0,alpha(j)))-I;
end
%Figure 7 Left
figure(10)
plot(alpha, RR, '--' , 'LineWidth', 2, 'color', 'black')
hold on
plot(alpha, alpha, '-' , 'LineWidth', 2, 'color', 'black')
ylim([0 0.3])
xlim([0 0.15])
set(gca,'FontSize',18)
xlabel('$$\alpha$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('Rejection rate', 'FontSize',25, 'interpreter', 'latex')
title('IV selection: size distortion','fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
lgd = legend('Under $p$-hacking', 'Nominal size','Orientation','vertical', 'interpreter', 'latex', 'location', 'northwest');
saveas(gcf,'BW/IVSize1_bw', 'epsc')
close all


% Figure 7 Right
h = linspace(0, 2.5, 1000);
alpha = 0.05;
B2m = (1/sqrt(2-sqrt(2)))*normpdf(h*sqrt((sqrt(2)-1)/sqrt(2))).*(normcdf(h/sqrt(2-sqrt(2))))+...
    sqrt(2)*normpdf(0)*(1-normcdf(sqrt(2)*h));
B2t = B2m - (1/sqrt(2-sqrt(2)))*normpdf(h*sqrt((sqrt(2)-1)/sqrt(2))).*(normcdf(-sqrt(4-2*sqrt(2))*z(0,alpha)+h/sqrt(2-sqrt(2))));
figure(11)
plot(h, B2t, '--' , 'LineWidth', 2, 'color', 'black')
hold on
plot(h, B2m, ':' , 'LineWidth', 2, 'color', 'black')
ylim([0 0.6])
xlim([0 2.5])
set(gca,'FontSize',18)
xlabel('$$h$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('Bias', 'FontSize',25, 'interpreter', 'latex')
title('IV selection: bias','fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
lgd = legend('Threshold', 'Minimum','Orientation','vertical', 'interpreter', 'latex');
saveas(gcf,'BW/IVBias1_bw', 'epsc')
close all

%%
%%%%%%%%% Data set selection %%%%%%%%%%%%
clear all
rng(12345)

alpha = 0.05;
z = @(h,p) norminv(1 - p) - h;
H = [0, 1, 2];

P = linspace(0.0001, 0.9999, 1000);
pcurves3_t = zeros(length(P), 3);
Bound = exp(z(0,P).^2/2).*(P<=0.5) + (P>0.5);
for d = 1:3
    h = H(d);
for j = 1:length(P)
    p = P(j);

Upsilon_3_t = (1+normcdf((z(h,alpha))))*(p<=alpha)+2*normcdf(z(h,p))*(p>alpha);
pcurves3_t(j,d) = exp(h*z(0,p) - h^2/2)*Upsilon_3_t;

end
end
%Figure 8
figure(12)
plot(P, pcurves3_t(:,1), '--' , 'LineWidth', 2, 'color', 'black')
hold on
plot(P, pcurves3_t(:,2), ':' , 'LineWidth', 2, 'color', 'black')
plot(P, pcurves3_t(:,3), '-.' , 'LineWidth', 2, 'color', 'black')
plot(P, Bound, '-' , 'LineWidth', 2, 'color', 'black')
ylim([0 15])
set(gca,'FontSize',18)
xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('$$g_3^t(p)$$', 'FontSize',25, 'interpreter', 'latex')
title('Dataset selection: threshold','fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
lgd = legend('$h=0$', '$h=1$','$h=2$','Bound', 'Orientation','horizontal', 'interpreter', 'latex');

axes('position',[.35 .34 .45 .45])
box on 
indexOfInterest = (P < 0.2) & (P > 0.0001); 
plot(P(indexOfInterest),pcurves3_t(indexOfInterest,1), '--' , 'LineWidth', 2, 'color', 'black') % plot on new axes
hold on
plot(P(indexOfInterest),pcurves3_t(indexOfInterest,2), ':' , 'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest),pcurves3_t(indexOfInterest,3), '-.' , 'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest),Bound(indexOfInterest), '-' , 'LineWidth', 2, 'color', 'black')
ylim([0 10])
xlim([0 0.1])
saveas(gcf,'BW/DataChoice1_bw', 'epsc')
close all

K = [2, 5, 20, 2, 5, 20];
H = [0,0,0,1,1,1];
P = linspace(0.0001, 0.9999, 1000);
pcurves3_m = zeros(length(P), 6);
Bound = exp(z(0,P).^2/2).*(P<=0.5) + (P>0.5);
for d = 1:6
    h = H(d);
    k = K(d);
for j = 1:length(P)
    p = P(j);

Upsilon_3_m = k*normcdf(z(h,p))^(k-1);
pcurves3_m(j,d) = exp(h*z(0,p) - h^2/2)*Upsilon_3_m;

end
end

%Figure 9 Left
figure(13)
plot(P, pcurves3_m(:,1), '--' , 'LineWidth', 2, 'color', 'black')
hold on
plot(P, pcurves3_m(:,2), ':' , 'LineWidth', 2, 'color', 'black')
plot(P, pcurves3_m(:,3), '-.' , 'LineWidth', 2, 'color', 'black')
plot(P, Bound, '-' , 'LineWidth', 2, 'color', 'black')
ylim([0 15])
set(gca,'FontSize',18)
xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('$$g_3^m(p;K)$$', 'FontSize',25, 'interpreter', 'latex')
title('Dataset selection: minimum, $h=0$','fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
lgd = legend('$K=2$', '$K=5$','$K=20$','Bound', 'Orientation','horizontal', 'interpreter', 'latex');

axes('position',[.35 .34 .45 .45])
box on 
indexOfInterest = (P < 0.4) & (P > 0.0001);
plot(P(indexOfInterest),pcurves3_m(indexOfInterest,1), '--' , 'LineWidth', 2, 'color', 'black') % plot on new axes
hold on
plot(P(indexOfInterest),pcurves3_m(indexOfInterest,2), ':' , 'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest),pcurves3_m(indexOfInterest,3), '-.' , 'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest),Bound(indexOfInterest), '-' , 'LineWidth', 2, 'color', 'black')
ylim([0 10])
xlim([0 0.25])
saveas(gcf,'BW/DataChoice2_bw', 'epsc')
close all

%Figure 9 Right
figure(14)
plot(P, pcurves3_m(:,4), '--' , 'LineWidth', 2, 'color', 'black')
hold on
plot(P, pcurves3_m(:,5), ':' , 'LineWidth', 2, 'color', 'black')
plot(P, pcurves3_m(:,6), '-.' , 'LineWidth', 2, 'color', 'black')
plot(P, Bound, '-' , 'LineWidth', 2, 'color', 'black')
ylim([0 15])
set(gca,'FontSize',18)
xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('$$g_3^m(p;K)$$', 'FontSize',25, 'interpreter', 'latex')
title('Dataset selection: minimum, $h=1$','fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
lgd = legend('$K=2$', '$K=5$','$K=20$','Bound', 'Orientation','horizontal', 'interpreter', 'latex');

axes('position',[.35 .34 .45 .45])
box on 
indexOfInterest = (P < 0.4) & (P > 0.0001); 
plot(P(indexOfInterest),pcurves3_m(indexOfInterest,1), '--' , 'LineWidth', 2, 'color', 'black') % plot on new axes
hold on
plot(P(indexOfInterest),pcurves3_m(indexOfInterest,2), ':' , 'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest),pcurves3_m(indexOfInterest,3), '-.' , 'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest),Bound(indexOfInterest), '-' , 'LineWidth', 2, 'color', 'black')
ylim([0 10])
xlim([0 0.25])
saveas(gcf,'BW/DataChoice3_bw', 'epsc')
close all
%%
%%%%%%%%%%%%%%%%%%%%%% Lag Length Selection %%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
rng(12345)
%For Lag Length selection
n = 200;
M = 1000000;
rho = zeros(M,1);
for m = 1:M
X = randn(n,1);
rho(m) = sum((X(2:n) - mean(X)).*(X(1:(n-1)) - mean(X)))/(n-1);
end
alpha = 0.05;
kappa = 1/2;
L = -1/(2*kappa);
z = @(h,p) norminv(1 - p) - h;
omega = @(r) sqrt(max(1+2*kappa*r, 0));
H0 = mean(rho<=0);
HL = mean(rho<=L);


P = linspace(0.0001, 0.9999, 1000);

pcurve_t = zeros(length(P),3);
pcurve_m = zeros(length(P),3);

pcurve_sel_t = zeros(length(P),3);
pcurve_sel_m = zeros(length(P),3);
Bound = exp(z(0,P).^2/2).*(P<=0.5) + (P>0.5);
for h = 0:2
for j = 1:length(P)
    p = P(j);
lp = L*(1 - (z(0, alpha)/z(0,p))^2);
g = omega(rho).*normpdf(z(0,p)*omega(rho)-h);
if (p<=alpha)
I = mean(g.*(rho<=lp).*(rho>L));
Upsilon4_t = 1+I/normpdf(z(h,p));
end
if (p>alpha)&&(p<=0.5)
I = mean(g.*(rho<=0).*(rho>L));
Upsilon4_t = 1-H0+HL+I/normpdf(z(h,p));
end
if (p>0.5)
I = mean(g.*(rho>0));
Upsilon4_t = H0+I/normpdf(z(h,p));
Upsilon4_m = Upsilon4_t;
end

if (p<=0.5)
I = mean(g.*(rho<=0).*(rho>L));
Upsilon4_m = 1-H0+HL+I/normpdf(z(h,p));
end


pcurve_t(j,h+1) = exp(h*z(0,p) - h^2/2)*Upsilon4_t;
pcurve_m(j,h+1) = exp(h*z(0,p) - h^2/2)*Upsilon4_m;
end
end

%Figure 10 Left
figure(15)
plot(P, pcurve_t(:,1), '--' , 'LineWidth', 2, 'color', 'black')
hold on
plot(P, pcurve_t(:,2), ':' , 'LineWidth', 2, 'color', 'black')
plot(P, pcurve_t(:,3), '-.' , 'LineWidth', 2, 'color', 'black')
plot(P, Bound, '-' , 'LineWidth', 2, 'color', 'black')
ylim([0 15])
set(gca,'FontSize',18)
xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('$$g_4^t(p)$$', 'FontSize',25, 'interpreter', 'latex')
title('Lag length selection: threshold','fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
lgd = legend('$h=0$', '$h=1$','$h=2$','Bound', 'Orientation','horizontal', 'interpreter', 'latex');

axes('position',[.35 .34 .45 .45])
box on 
indexOfInterest = (P < 0.2) & (P > 0.0001); 
plot(P(indexOfInterest),pcurve_t(indexOfInterest,1), '--' , 'LineWidth', 2, 'color', 'black') % plot on new axes
hold on
plot(P(indexOfInterest),pcurve_t(indexOfInterest,2), ':' , 'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest),pcurve_t(indexOfInterest,3), '-.' , 'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest),Bound(indexOfInterest), '-' , 'LineWidth', 2, 'color', 'black')
ylim([0 10])
xlim([0 0.1])
saveas(gcf,'BW/LagLengthT_bw', 'epsc')
close all

%Figure 10 Right
figure(16)
plot(P, pcurve_m(:,1), '--' , 'LineWidth', 2, 'color', 'black')
hold on
plot(P, pcurve_m(:,2), ':' , 'LineWidth', 2, 'color', 'black')
plot(P, pcurve_m(:,3), '-.' , 'LineWidth', 2, 'color', 'black')
plot(P, Bound, '-' , 'LineWidth', 2, 'color', 'black')
ylim([0 15])
set(gca,'FontSize',18)
xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('$$g_4^m(p)$$', 'FontSize',25, 'interpreter', 'latex')
title('Lag length selection: minimum','fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
lgd = legend('$h=0$', '$h=1$','$h=2$','Bound', 'Orientation','horizontal', 'interpreter', 'latex');

axes('position',[.35 .34 .45 .45])
box on
indexOfInterest = (P < 0.6) & (P > 0.0001); 
plot(P(indexOfInterest),pcurve_m(indexOfInterest,1), '--' , 'LineWidth', 2, 'color', 'black') % plot on new axes
hold on
plot(P(indexOfInterest),pcurve_m(indexOfInterest,2), ':' , 'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest),pcurve_m(indexOfInterest,3), '-.' , 'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest),Bound(indexOfInterest), '-' , 'LineWidth', 2, 'color', 'black')
ylim([0 2])
xlim([0 0.6])
saveas(gcf,'BW/LagLengthM_bw', 'epsc')
close all

%% rejection rate
alpha = linspace(0.0001, 0.15, 1000);
RR = zeros(length(alpha), 1);
for j = 1:length(alpha)
g = normcdf(z(0,alpha(j))*omega(rho));
RR(j) = alpha(j) + (1-alpha(j))*(H0-HL) - mean(g.*(rho<=0).*(rho>L));
end

%Figure 11
figure(17)
plot(alpha, RR, '--' , 'LineWidth', 2, 'color', 'black')
hold on
plot(alpha, alpha, '-' , 'LineWidth', 2, 'color', 'black')
ylim([0 0.15])
xlim([0 0.15])
set(gca,'FontSize',18)
xlabel('$$\alpha$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('Rejection rate', 'FontSize',25, 'interpreter', 'latex')
title('Lag length selection: size distortion','fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
lgd = legend('Under $p$-hacking', 'Nominal size', 'Orientation','vertical', 'interpreter', 'latex','location','northwest');
saveas(gcf,'BW/LagLengthSize1_bw', 'epsc')
close all