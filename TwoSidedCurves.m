clear all
rng(12345)


%bound
h_grid = linspace(0, 6, 10000);
cv2 = @(p) norminv(1-p/2);
p = linspace(0,1, 203);
for j = 1:length(p)
bound(j) = max((exp(h_grid*cv2(p(j)) - h_grid.^2/2) + exp(-h_grid*cv2(p(j)) - h_grid.^2/2))/2);
end

fig=0;
alpha = 0.05;
N = 10^6;
gamma = [0.1, 0.5, 0.9,1];
Rho  =1-gamma.^2;
H = [0,1,2];
P0_s = struct();
Pr_s = struct();
Pr_min_s = struct();
M = 100;

for k = 1:4
    rho = Rho(k);
    %T2 = rho*T1 + sqrt(1-rho^2)*randn(N,1);
for j = 1:3

    for m = 0:M
T1 = randn(N,1);
T2 = rho*T1 + sqrt(1-rho^2)*randn(N,1);
h = H(j);
T1h=abs(T1+h);
T2h = abs(T2+h);
p1 = 2*normcdf(-T1h);
p2 = 2*normcdf(-T2h);
% T1h=(T1+h);
% T2h = (T2+h);
% p1 = normcdf(-T1h);
% p2 = normcdf(-T2h);
P0 = p1;
Pr_min = (p2.*(p2<=p1) + p1.*(p1<p2));
Pr = p1.*(p1<=alpha) + (p1>alpha).*Pr_min;
[edges0, N0] = EdgesN(P0);
[edges1, N1] = EdgesN(Pr);
[edges1min, N1min] = EdgesN(Pr_min);
if m == 0
    K0 = 400;
P0_s.(['h' num2str(H(j)) 'g' num2str(gamma(k)*10)]) = zeros(K0, 2);
Pr_s.(['h' num2str(H(j)) 'g' num2str(gamma(k)*10)]) = zeros(K0, 2);
Pr_min_s.(['h' num2str(H(j)) 'g' num2str(gamma(k)*10)]) = zeros(K0, 2);

end
if m>0
P0_s.(['h' num2str(H(j)) 'g' num2str(gamma(k)*10)]) = P0_s.(['h' num2str(H(j)) 'g' num2str(gamma(k)*10)])+[edges0', N0']/M;
Pr_s.(['h' num2str(H(j)) 'g' num2str(gamma(k)*10)]) = Pr_s.(['h' num2str(H(j)) 'g' num2str(gamma(k)*10)])+ [edges1', N1']/M;
Pr_min_s.(['h' num2str(H(j)) 'g' num2str(gamma(k)*10)]) =Pr_min_s.(['h' num2str(H(j)) 'g' num2str(gamma(k)*10)])+ [edges1min', N1min']/M;
end
end
end
end


close all
fig = fig+1;
figure(fig)
hold on
plot(Pr_s.h0g5(:,1), Pr_s.h0g5(:,2), 'LineStyle','--', 'LineWidth',2, 'Color', 'black')
plot(Pr_s.h1g5(:,1), Pr_s.h1g5(:,2), 'LineStyle',':', 'LineWidth',2, 'Color', 'black')
plot(Pr_s.h2g5(:,1), Pr_s.h2g5(:,2), 'LineStyle','-.', 'LineWidth',2, 'Color', 'black')
plot(p, bound, 'LineStyle','-', 'LineWidth',2, 'Color', 'black')
ylim([0 15])
set(gca,'FontSize',18)
xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('$$g_1^t(p)$$', 'FontSize',25, 'interpreter', 'latex')
title('Covariate selection: threshold, $\gamma=0.5$','fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
lgd = legend('$h=0$', '$h=1$','$h=2$','Bound', 'Orientation','horizontal', 'interpreter', 'latex');

axes('position',[.35 .34 .45 .45])
box on
hold on
indexOfInterest = (edges0 < 0.2) & (edges0 > 0.0001);
plot(Pr_s.h0g5(indexOfInterest,1), Pr_s.h0g5(indexOfInterest,2), 'LineStyle','--', 'LineWidth',2, 'Color', 'black')
plot(Pr_s.h1g5(indexOfInterest,1), Pr_s.h1g5(indexOfInterest,2), 'LineStyle',':', 'LineWidth',2, 'Color', 'black')
plot(Pr_s.h2g5(indexOfInterest,1), Pr_s.h2g5(indexOfInterest,2), 'LineStyle','-.', 'LineWidth',2, 'Color', 'black')
plot(p(indexOfInterest), bound(indexOfInterest), 'LineStyle','-', 'LineWidth',2, 'Color', 'black')
xlim([0 0.1])
ylim([0, 12])
saveas(gcf,'BW/CovarChoice1_bw_2sided', 'epsc')
close all

%%
close all
fig = fig+1;
figure(fig)
hold on
plot(Pr_s.h1g1(:,1), Pr_s.h1g1(:,2), 'LineStyle','--', 'LineWidth',2, 'Color', 'black')
plot(Pr_s.h1g5(:,1), Pr_s.h1g5(:,2), 'LineStyle',':', 'LineWidth',2, 'Color', 'black')
plot(Pr_s.h1g9(:,1), Pr_s.h1g9(:,2), 'LineStyle','-.', 'LineWidth',2, 'Color', 'black')
plot(p, bound, 'LineStyle','-', 'LineWidth',2, 'Color', 'black')
ylim([0 15])
set(gca,'FontSize',18)
xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('$$g_1^t(p)$$', 'FontSize',25, 'interpreter', 'latex')
title('Covariate selection: threshold, $h=1$','fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
lgd = legend('$\gamma=0.1$', '$\gamma=0.5$','$\gamma=0.9$','Bound', 'Orientation','horizontal', 'interpreter', 'latex');

axes('position',[.35 .34 .45 .45])
box on
hold on
indexOfInterest = (edges0 < 0.2) & (edges0 > 0.0001);
plot(Pr_s.h1g1(indexOfInterest,1), Pr_s.h1g1(indexOfInterest,2), 'LineStyle','--', 'LineWidth',2, 'Color', 'black')
plot(Pr_s.h1g5(indexOfInterest,1), Pr_s.h1g5(indexOfInterest,2), 'LineStyle',':', 'LineWidth',2, 'Color', 'black')
plot(Pr_s.h1g9(indexOfInterest,1), Pr_s.h1g9(indexOfInterest,2), 'LineStyle','-.', 'LineWidth',2, 'Color', 'black')
plot(p(indexOfInterest), bound(indexOfInterest), 'LineStyle','-', 'LineWidth',2, 'Color', 'black')
xlim([0 0.1])
ylim([0, 15])
saveas(gcf,'BW/CovarChoice2_bw_2sided', 'epsc')
close all
%%
close all
fig = fig+1;
figure(fig)
hold on
plot(Pr_s.h1g5(:,1), Pr_s.h1g5(:,2), 'LineStyle','--', 'LineWidth',2, 'Color', 'black')
plot(Pr_min_s.h1g5(:,1), Pr_min_s.h1g5(:,2), 'LineStyle',':', 'LineWidth',2, 'Color', 'black')
plot(P0_s.h1g5(:,1), P0_s.h1g5(:,2), 'LineStyle','-.', 'LineWidth',2, 'Color', 'black')
plot(p, bound, 'LineStyle','-', 'LineWidth',2, 'Color', 'black')
ylim([0 15])
set(gca,'FontSize',18)
xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('$$g_1(p)$$', 'FontSize',25, 'interpreter', 'latex')
title('Covariate selection: threshold vs. minimum, $h = 1$','fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
lgd = legend('Threshold', 'Minimum','No $p$-hacking','Bound', 'Orientation','horizontal', 'interpreter', 'latex', 'NumColumns',2);


axes('position',[.35 .34 .45 .45])
box on
hold on
indexOfInterest = (edges0 < 0.2) & (edges0 > 0.0001);
plot(Pr_s.h1g5(indexOfInterest,1), Pr_s.h1g5(indexOfInterest,2), 'LineStyle','--', 'LineWidth',2, 'Color', 'black')
plot(Pr_min_s.h1g5(indexOfInterest,1), Pr_min_s.h1g5(indexOfInterest,2), 'LineStyle',':', 'LineWidth',2, 'Color', 'black')
plot(P0_s.h1g5(indexOfInterest,1), P0_s.h1g5(indexOfInterest,2), 'LineStyle','-.', 'LineWidth',2, 'Color', 'black')
plot(p(indexOfInterest), bound(indexOfInterest), 'LineStyle','-', 'LineWidth',2, 'Color', 'black')
xlim([0 0.1])
ylim([0, 15])
saveas(gcf,'BW/CovarChoice3_bw_2sided', 'epsc')
close all


%% IV

M = 100;
for j = 1:3
    for m = 0:M
        W1 = randn(N, 1);
W2 = randn(N, 1);
W12 = (W1+W2)/sqrt(2);

h = H(j);
T1h=abs(W1+h);
T2h = abs(W2+h);
T12h = abs(W12+sqrt(2)*h);

p1 = 2*normcdf(-T1h);
p2 = 2*normcdf(-T2h);
p12 = 2*normcdf(-T12h);
% T1h=(W1+h);
% T2h = (W2+h);
% T12h = (W12+sqrt(2)*h);
% 
% p1 = normcdf(-T1h);
% p2 = normcdf(-T2h);
% p12 = normcdf(-T12h);

P0 = p12;
Pr_min = min([p1, p2, p12]')';
Pr = P0.*(P0<=alpha) + (P0>alpha).*Pr_min;
[edges0, N0] = EdgesN(P0);
[edges1, N1] = EdgesN(Pr);
[edges1min, N1min] = EdgesN(Pr_min);
if m == 0
P0_s.(['h' num2str(H(j)) 'iv']) = zeros(K0, 2);
Pr_s.(['h' num2str(H(j)) 'iv']) = zeros(K0, 2);
Pr_min_s.(['h' num2str(H(j)) 'iv']) = zeros(K0, 2);

end
if m>0
P0_s.(['h' num2str(H(j)) 'iv']) = P0_s.(['h' num2str(H(j)) 'iv'])+[edges0', N0']/M;
Pr_s.(['h' num2str(H(j)) 'iv']) = Pr_s.(['h' num2str(H(j)) 'iv'])+ [edges1', N1']/M;
Pr_min_s.(['h' num2str(H(j)) 'iv']) =Pr_min_s.(['h' num2str(H(j)) 'iv'])+ [edges1min', N1min']/M;
end
    end 
end


% for j = 1:length(p)
% g2(j) = mean(normpdf((Pr - p(j))/0.0001))/0.0001;
% end

%%
close all
fig = fig+1;
figure(fig)
hold on
plot(Pr_s.h0iv(:,1), Pr_s.h0iv(:,2), 'LineStyle','--', 'LineWidth',2, 'Color', 'black')
plot(Pr_s.h1iv(:,1), Pr_s.h1iv(:,2), 'LineStyle',':', 'LineWidth',2, 'Color', 'black')
plot(Pr_s.h2iv(:,1), Pr_s.h2iv(:,2), 'LineStyle','-.', 'LineWidth',2, 'Color', 'black')
plot(p, bound, 'LineStyle','-', 'LineWidth',2, 'Color', 'black')
ylim([0 15])
set(gca,'FontSize',18)
xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('$$g_2^t(p)$$', 'FontSize',25, 'interpreter', 'latex')
title('IV selection: threshold','fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
lgd = legend('$h=0$', '$h=1$','$h=2$','Bound', 'Orientation','horizontal', 'interpreter', 'latex');
axes('position',[.35 .34 .45 .45])
box on
hold on
indexOfInterest = (edges0 < 0.2) & (edges0 > 0.0001);
plot(Pr_s.h0iv(indexOfInterest,1), Pr_s.h0iv(indexOfInterest,2), 'LineStyle','--', 'LineWidth',2, 'Color', 'black')
plot(Pr_s.h1iv(indexOfInterest,1), Pr_s.h1iv(indexOfInterest,2), 'LineStyle',':', 'LineWidth',2, 'Color', 'black')
plot(Pr_s.h2iv(indexOfInterest,1), Pr_s.h2iv(indexOfInterest,2), 'LineStyle','-.', 'LineWidth',2, 'Color', 'black')
plot(p(indexOfInterest), bound(indexOfInterest), 'LineStyle','-', 'LineWidth',2, 'Color', 'black')
xlim([0 0.1])
ylim([0, 10])
saveas(gcf,'BW/IVChoice1_bw_2sided', 'epsc')
close all

%%
close all
fig = fig+1;
figure(fig)
hold on
plot(Pr_min_s.h0iv(:,1), Pr_min_s.h0iv(:,2), 'LineStyle','--', 'LineWidth',2, 'Color', 'black')
plot(Pr_min_s.h1iv(:,1), Pr_min_s.h1iv(:,2), 'LineStyle',':', 'LineWidth',2, 'Color', 'black')
plot(Pr_min_s.h2iv(:,1), Pr_min_s.h2iv(:,2), 'LineStyle','-.', 'LineWidth',2, 'Color', 'black')
plot(p, bound, 'LineStyle','-', 'LineWidth',2, 'Color', 'black')
ylim([0 15])
set(gca,'FontSize',18)
xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('$$g_2^m(p)$$', 'FontSize',25, 'interpreter', 'latex')
title('IV selection: minimum','fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
lgd = legend('$h=0$', '$h=1$','$h=2$','Bound', 'Orientation','horizontal', 'interpreter', 'latex');

axes('position',[.35 .34 .45 .45])
box on
hold on
indexOfInterest = (edges0 < 0.2) & (edges0 > 0.0001);
plot(Pr_min_s.h0iv(indexOfInterest,1), Pr_min_s.h0iv(indexOfInterest,2), 'LineStyle','--', 'LineWidth',2, 'Color', 'black')
plot(Pr_min_s.h1iv(indexOfInterest,1), Pr_min_s.h1iv(indexOfInterest,2), 'LineStyle',':', 'LineWidth',2, 'Color', 'black')
plot(Pr_min_s.h2iv(indexOfInterest,1), Pr_min_s.h2iv(indexOfInterest,2), 'LineStyle','-.', 'LineWidth',2, 'Color', 'black')
plot(p(indexOfInterest), bound(indexOfInterest), 'LineStyle','-', 'LineWidth',2, 'Color', 'black')
xlim([0 0.1])
ylim([0, 10])
saveas(gcf,'BW/IVChoice2_bw_2sided', 'epsc')
close all

%%
close all
fig = fig+1;
figure(fig)
hold on
plot(Pr_s.h1iv(:,1), Pr_s.h1iv(:,2), 'LineStyle','--', 'LineWidth',2, 'Color', 'black')
plot(Pr_min_s.h1iv(:,1), Pr_min_s.h1iv(:,2), 'LineStyle',':', 'LineWidth',2, 'Color', 'black')
plot(P0_s.h1iv(:,1), P0_s.h1iv(:,2), 'LineStyle','-.', 'LineWidth',2, 'Color', 'black')
plot(p, bound, 'LineStyle','-', 'LineWidth',2, 'Color', 'black')
ylim([0 15])
set(gca,'FontSize',18)
xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('$$g_2(p)$$', 'FontSize',25, 'interpreter', 'latex')
title('IV selection: threshold vs. minimum, $h = 1$','fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
lgd = legend('Threshold', 'Minimum','No $p$-hacking','Bound', 'Orientation','horizontal', 'interpreter', 'latex', 'NumColumns',2);


axes('position',[.35 .34 .45 .45])
box on
hold on
indexOfInterest = (edges0 < 0.2) & (edges0 > 0.0001);
plot(Pr_s.h1iv(indexOfInterest,1), Pr_s.h1iv(indexOfInterest,2), 'LineStyle','--', 'LineWidth',2, 'Color', 'black')
plot(Pr_min_s.h1iv(indexOfInterest,1), Pr_min_s.h1iv(indexOfInterest,2), 'LineStyle',':', 'LineWidth',2, 'Color', 'black')
plot(P0_s.h1iv(indexOfInterest,1), P0_s.h1iv(indexOfInterest,2), 'LineStyle','-.', 'LineWidth',2, 'Color', 'black')
plot(p(indexOfInterest), bound(indexOfInterest), 'LineStyle','-', 'LineWidth',2, 'Color', 'black')
xlim([0 0.1])
ylim([0, 15])
saveas(gcf,'BW/IVChoice3_bw_2sided', 'epsc')
close all

%%
close all
fig = fig+1;
figure(fig)
hold on
plot(Pr_s.h0g10(:,1), Pr_s.h0g10(:,2), 'LineStyle','--', 'LineWidth',2, 'Color', 'black')
plot(Pr_s.h1g10(:,1), Pr_s.h1g10(:,2), 'LineStyle',':', 'LineWidth',2, 'Color', 'black')
plot(Pr_s.h2g10(:,1), Pr_s.h2g10(:,2), 'LineStyle','-.', 'LineWidth',2, 'Color', 'black')
plot(p, bound, 'LineStyle','-', 'LineWidth',2, 'Color', 'black')
ylim([0 15])
set(gca,'FontSize',18)
xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('$$g_3^t(p)$$', 'FontSize',25, 'interpreter', 'latex')
title('Dataset selection: threshold','fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
lgd = legend('$h=0$', '$h=1$','$h=2$','Bound', 'Orientation','horizontal', 'interpreter', 'latex');

axes('position',[.35 .34 .45 .45])
box on
hold on
indexOfInterest = (edges0 < 0.2) & (edges0 > 0.0001);
plot(Pr_s.h0g10(indexOfInterest,1), Pr_s.h0g10(indexOfInterest,2), 'LineStyle','--', 'LineWidth',2, 'Color', 'black')
plot(Pr_s.h1g10(indexOfInterest,1), Pr_s.h1g10(indexOfInterest,2), 'LineStyle',':', 'LineWidth',2, 'Color', 'black')
plot(Pr_s.h2g10(indexOfInterest,1), Pr_s.h2g10(indexOfInterest,2), 'LineStyle','-.', 'LineWidth',2, 'Color', 'black')
plot(p(indexOfInterest), bound(indexOfInterest), 'LineStyle','-', 'LineWidth',2, 'Color', 'black')
ylim([0 10])
xlim([0 0.1])
saveas(gcf,'BW/DataChoice1_bw_2sided', 'epsc')
close all


%minimum
z = @(h,p) norminv(1 - p) - h;
K = [2, 5, 20, 2, 5, 20];
H = [0,0,0,1,1,1];
P = linspace(0.0001, 0.9999, 1000);
pcurves3_m = zeros(length(P), 6);
Bound = zeros(size(P));
for j = 1:length(P)
Bound(j) = max((exp(h_grid*cv2(P(j)) - h_grid.^2/2) + exp(-h_grid*cv2(P(j)) - h_grid.^2/2))/2);
end
for d = 1:6
    h = H(d);
    k = K(d);
for j = 1:length(P)
    p = P(j);

Upsilon_3_m = k*(normcdf(z(-h,p/2)) - normcdf(-z(h,p/2)))^(k-1);
pcurves3_m(j,d) = 0.5*(exp(h*z(0,p/2) - h^2/2)+exp(-h*z(0,p/2) - h^2/2))*Upsilon_3_m;

end
end

%Figure 9 Left
close all
fig = fig+1;
figure(fig)
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
saveas(gcf,'BW/DataChoice2_bw_2sided', 'epsc')
close all

%Figure 9 Right
close all
fig = fig+1;
figure(fig)
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
plot(P(indexOfInterest),pcurves3_m(indexOfInterest,4), '--' , 'LineWidth', 2, 'color', 'black') % plot on new axes
hold on
plot(P(indexOfInterest),pcurves3_m(indexOfInterest,5), ':' , 'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest),pcurves3_m(indexOfInterest,6), '-.' , 'LineWidth', 2, 'color', 'black')
plot(P(indexOfInterest),Bound(indexOfInterest), '-' , 'LineWidth', 2, 'color', 'black')
ylim([0 10])
xlim([0 0.25])
saveas(gcf,'BW/DataChoice3_bw_2sided', 'epsc')
close all


%%%Lag length

Mc = 10^7;
Pr_min = zeros(Mc, 3);
Pr = zeros(Mc, 3);
N = 200;
H = [0,1,2];
for m = 1:(Mc)
U = randn(200,1);
for j = 1:3
    h = H(j);
T0 = abs(h + sqrt(N)*mean(U));
%T0 = (h + sqrt(N)*mean(U));
rho_hat = mean((U(2:end)-mean(U)).*(U(1:(end-1))-mean(U)));
omega2_hat = (1 + rho_hat);
if (omega2_hat<=0)
T1 = T0;
end
if omega2_hat>0
T1 = T0/sqrt(omega2_hat);
end
p0 = 2*normcdf(-T0);
p1 = 2*normcdf(-T1);
%p0 = normcdf(-T0);
%p1 = normcdf(-T1);


Pr_min(m,j) = min(p0, p1);
Pr(m,j) = p0*(p0<alpha) + (p0>alpha)*min(p0,p1);
end
end


[edges10, N10] = EdgesN(Pr(:,1));
[edges11, N11] = EdgesN(Pr(:,2));
[edges12, N12] = EdgesN(Pr(:,3));
[edges1min0, N1min0] = EdgesN(Pr_min(:,1));
[edges1min1, N1min1] = EdgesN(Pr_min(:,2));
[edges1min2, N1min2] = EdgesN(Pr_min(:,3));


edges0 = edges10;
h_grid = linspace(0, 6, 10000);
cv2 = @(p) norminv(1-p/2);
p = edges10;
bound = zeros(size(p));
for j = 1:length(p)
bound(j) = max((exp(h_grid*cv2(p(j)) - h_grid.^2/2) + exp(-h_grid*cv2(p(j)) - h_grid.^2/2))/2);
end


close all
fig = fig+1;
figure(fig)
hold on
plot(edges10, N10, 'LineStyle','--', 'LineWidth',2, 'Color', 'black')
plot(edges11, N11, 'LineStyle',':', 'LineWidth',2, 'Color', 'black')
plot(edges12, N12, 'LineStyle','-.', 'LineWidth',2, 'Color', 'black')
plot(p, bound, 'LineStyle','-', 'LineWidth',2, 'Color', 'black')
ylim([0 15])
set(gca,'FontSize',18)
xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('$$g_4^t(p)$$', 'FontSize',25, 'interpreter', 'latex')
title('Lag length selection: threshold','fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
lgd = legend('$h=0$', '$h=1$','$h=2$','Bound', 'Orientation','horizontal', 'interpreter', 'latex');

axes('position',[.35 .34 .45 .45])
box on
hold on
indexOfInterest = (edges0 < 0.2) & (edges0 > 0.0001);
plot(edges10(indexOfInterest), N10(indexOfInterest), 'LineStyle','--', 'LineWidth',2, 'Color', 'black')
plot(edges11(indexOfInterest), N11(indexOfInterest), 'LineStyle',':', 'LineWidth',2, 'Color', 'black')
plot(edges12(indexOfInterest), N12(indexOfInterest), 'LineStyle','-.', 'LineWidth',2, 'Color', 'black')

plot(p(indexOfInterest), bound(indexOfInterest), 'LineStyle','-', 'LineWidth',2, 'Color', 'black')
ylim([0 10])
xlim([0 0.1])
saveas(gcf,'BW/LagLengthT_bw_2sided', 'epsc')
close all

%%
close all
fig = fig+1;
figure(fig)
hold on
plot(edges1min0, N1min0, 'LineStyle','--', 'LineWidth',2, 'Color', 'black')
plot(edges1min1, N1min1, 'LineStyle',':', 'LineWidth',2, 'Color', 'black')
plot(edges1min2, N1min2, 'LineStyle','-.', 'LineWidth',2, 'Color', 'black')
plot(p, bound, 'LineStyle','-', 'LineWidth',2, 'Color', 'black')
ylim([0 15])
set(gca,'FontSize',18)
xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('$$g_4^m(p)$$', 'FontSize',25, 'interpreter', 'latex')
title('Lag length selection: minimum','fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
lgd = legend('$h=0$', '$h=1$','$h=2$','Bound', 'Orientation','horizontal', 'interpreter', 'latex');

axes('position',[.35 .34 .45 .45])
box on
hold on
indexOfInterest = (edges0 < 0.6) & (edges0 > 0.0001);
plot(edges1min0(indexOfInterest), N1min0(indexOfInterest), 'LineStyle','--', 'LineWidth',2, 'Color', 'black')
plot(edges1min1(indexOfInterest), N1min1(indexOfInterest), 'LineStyle',':', 'LineWidth',2, 'Color', 'black')
plot(edges1min2(indexOfInterest), N1min2(indexOfInterest), 'LineStyle','-.', 'LineWidth',2, 'Color', 'black')

plot(p(indexOfInterest), bound(indexOfInterest), 'LineStyle','-', 'LineWidth',2, 'Color', 'black')
ylim([0 2])
xlim([0 0.6])
saveas(gcf,'BW/LagLengthM_bw_2sided', 'epsc')
close all