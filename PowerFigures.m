clear all
S = importdata('PowerCurves/RejectionRates_April6.csv');
V = S.data;

tau = [0:0.05:1];

H = [0, 1, 2, 3];
K = [3,5,7];
Up = [0,1];

fig = 0;

for minimum = 0:1
for j = 1:4
for k = 1:3
    fig=fig+1;

figure (fig)
plot(tau, V(:, 2+(fig-1)*9),'LineWidth',3,'color' , 'black')
hold on
%plot(tau, V(:,3+ (fig-1)*9),'-x','LineWidth',3, 'color', 'black','MarkerSize',12)
plot(tau, V(:, 4+ (fig-1)*9),'-.','LineWidth',3, 'color', 'blue','MarkerSize',12)
plot(tau, V(:, 5+ (fig-1)*9),'--','LineWidth',3, 'color', [0 0.5 0],'MarkerSize',12)
plot(tau, V(:,6+ (fig-1)*9),'-.o','LineWidth',3, 'color', [0 0.5 0],'MarkerSize',12)
plot(tau, V(:,7+ (fig-1)*9),'-.*','LineWidth',3, 'color', [0 0.5 0],'MarkerSize',12)
plot(tau, V(:,8+ (fig-1)*9),'-.x','LineWidth',3, 'color', 'red','MarkerSize',12)

lgd = legend('Binomial', 'Discontinuity',...
'CS1','CSUB', 'CS2B','LCM', 'Location','northeast');
if (k==1)&&(j==1)
lgd = legend('Binomial', 'Discontinuity',...
'CS1','CSUB', 'CS2B','LCM', 'Location','northwest');
end
if (k==1)&&(j==2)
lgd = legend('Binomial', 'Discontinuity',...
'CS1','CSUB', 'CS2B','LCM', 'Location','northwest');
end
if (k==1)&&(j==4)
lgd = legend('Binomial', 'Discontinuity',...
'CS1','CSUB', 'CS2B','LCM', 'Location','northwest');
end

lgd.FontSize = 16;

set(gca,'FontSize',18)
xlabel('$$\tau$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('Rejection Rate', 'FontSize',25, 'interpreter', 'latex')
ylim([0 1])
xlim([0 1])
title(append('$K = $ ', int2str(K(k)), ', $h = $ ', int2str(H(j))),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')

if (j==4)
    title(append('$K = $ ', int2str(K(k)), ', $h \sim \chi^2(1)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end
if (minimum == 0)
saveas(gcf,append('PowerCurves/CovariateSelection/','power_sel', int2str(H(j)), int2str(K(k))), 'epsc')
close all
end
if (minimum == 1)
saveas(gcf,append('PowerCurves/CovariateSelection/','power_sel', int2str(H(j)), int2str(K(k)), 'min'), 'epsc')
close all
end
end
end

for j = 1:4
for k = 1:2
    fig=fig+1;

figure (fig)
plot(tau, V(:, 2+(fig-1)*9),'LineWidth',3,'color' , 'black')
hold on
%plot(tau, V(:,3+ (fig-1)*9),'-x','LineWidth',3, 'color', 'black','MarkerSize',12)
plot(tau, V(:, 4+ (fig-1)*9),'-.','LineWidth',3, 'color', 'blue','MarkerSize',12)
plot(tau, V(:, 5+ (fig-1)*9),'--','LineWidth',3, 'color', [0 0.5 0],'MarkerSize',12)
plot(tau, V(:,6+ (fig-1)*9),'-.o','LineWidth',3, 'color', [0 0.5 0],'MarkerSize',12)
plot(tau, V(:,7+ (fig-1)*9),'-.*','LineWidth',3, 'color', [0 0.5 0],'MarkerSize',12)
plot(tau, V(:,8+ (fig-1)*9),'-.x','LineWidth',3, 'color', 'red','MarkerSize',12)

lgd = legend('Binomial', 'Discontinuity',...
'CS1','CSUB', 'CS2B','LCM', 'Location','northeast');
if (j==4)
lgd = legend('Binomial', 'Discontinuity',...
'CS1','CSUB', 'CS2B','LCM', 'Location','northwest');
end
if (k==1)
lgd = legend('Binomial', 'Discontinuity',...
'CS1','CSUB', 'CS2B','LCM', 'Location','northwest');
end



lgd.FontSize = 16;

set(gca,'FontSize',18)
xlabel('$$\tau$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('Rejection Rate', 'FontSize',25, 'interpreter', 'latex')
ylim([0 1])
xlim([0 1])
title(append('$K = $ ', int2str(K(k)), ', $h = $ ', int2str(H(j))),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')

if (j==4)
    title(append('$K = $ ', int2str(K(k)), ', $h \sim \chi^2(1)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end

if (minimum==0)
saveas(gcf,append('PowerCurves/IVSelection/','power_iv', int2str(H(j)), int2str(K(k))), 'epsc')
close all
end
if (minimum==1)
saveas(gcf,append('PowerCurves/IVSelection/','power_iv', int2str(H(j)), int2str(K(k)), 'min'), 'epsc')
close all
end
end
end

for j = 1:4
for k = 1:1
    fig=fig+1;

figure (fig)
plot(tau, V(:, 2+(fig-1)*9),'LineWidth',3,'color' , 'black')
hold on
%plot(tau, V(:,3+ (fig-1)*9),'-x','LineWidth',3, 'color', 'black','MarkerSize',12)
plot(tau, V(:, 4+ (fig-1)*9),'-.','LineWidth',3, 'color', 'blue','MarkerSize',12)
plot(tau, V(:, 5+ (fig-1)*9),'--','LineWidth',3, 'color', [0 0.5 0],'MarkerSize',12)
plot(tau, V(:,6+ (fig-1)*9),'-.o','LineWidth',3, 'color', [0 0.5 0],'MarkerSize',12)
plot(tau, V(:,7+ (fig-1)*9),'-.*','LineWidth',3, 'color', [0 0.5 0],'MarkerSize',12)
plot(tau, V(:,8+ (fig-1)*9),'-.x','LineWidth',3, 'color', 'red','MarkerSize',12)

lgd = legend('Binomial', 'Discontinuity',...
'CS1','CSUB', 'CS2B','LCM', 'Location','northwest');



lgd.FontSize = 16;

set(gca,'FontSize',18)
xlabel('$$\tau$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('Rejection Rate', 'FontSize',25, 'interpreter', 'latex')
ylim([0 1])
xlim([0 1])
title(append('$h = $ ', int2str(H(j))),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
if (j==4)
    title(append('$h \sim \chi^2(1)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end

if (minimum == 0)
saveas(gcf,append('PowerCurves/LagLengthSelection/','power_var', int2str(H(j))), 'epsc')
close all
end
if (minimum == 1)
saveas(gcf,append('PowerCurves/LagLengthSelection/','power_var', int2str(H(j)), 'min'), 'epsc')
close all
end
end
end

end