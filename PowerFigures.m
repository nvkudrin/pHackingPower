clear all
%S = importdata('PowerCurves/RejectionRates_Apr5.csv');
%V = S.data;
S = readtable("PowerCurves/RejectionRates_Apr5.csv");
Vfull=table2array(S);
Vfull = Vfull(:,2:end);

tau = [0:0.05:1];
H = [0, 1, 2, 3];
K = [3,5,7];
num_per_dgp = 11;


%% Tables
A = Vfull([11;32;53], :);
column_names = ["Binomial", "Discontinuity", "CS1", "CSUB", "CS2B", "LCM"];
row_names = ["No Pub Bias", "Sharp Pub Bias", "Smooth Pub Bias"];
row_names = [" ", "Cov Selection (K = 3,thresholding)", row_names, "Cov Selection (K = 3, minimum)", row_names,...
              "IV Selection (K = 3,thresholding)", row_names, "IV Selection (K = 3, minimum)", row_names,...
              "Lag Selection (thresholding)", row_names, "Lag Selection (minimum)", row_names];
for sided = [1,2]
for h = H
  Tablel = NaN(24,6);
  row_num = 1;
for j = ([1+3*h,13+2*h,21+h]+48*(sided-1))
  sub_table_sel = A(:,((j-1)*num_per_dgp+1):(j*num_per_dgp));
  sub_table_sel(:,[1,3,9,10,11])=[];
  sub_table_min = A(:,((j-1)*num_per_dgp+1+24*num_per_dgp):(j*num_per_dgp+24*num_per_dgp));
  sub_table_min(:,[1,3,9,10,11])=[];
  Tablel((row_num+1):(row_num+3),:) = sub_table_sel;
  Tablel((row_num+5):(row_num+7),:) = sub_table_min;
  row_num = row_num+8;
end
Table_final = [row_names', [column_names; round(Tablel,3)]];
writematrix(Table_final,append('PowerCurves/Table_',int2str(sided), 'sided_',int2str(h),'.csv'))
end
end
%% Figures
%fig = 0;
for sided = [1,2]
    V = Vfull(:, (1+(sided-1)*1056/2): (1056/(3-sided)));
    if sided == 1
        side = "1sided";
    end
    if sided == 2
        side = "2sided";
    end
for pub_b = [0,1,2]
    fig = 0;
    start = pub_b*21+1;
    finish = 21*(pub_b+1);
for minimum = 0:1
for j = 1:4
for k = 1:3
    fig=fig+1;

figure (fig)
plot(tau, V(start:finish, 2+(fig-1)*num_per_dgp),'LineWidth',3,'color' , 'black')
hold on
%plot(tau, V(:,3+ (fig-1)*9),'-x','LineWidth',3, 'color', 'black','MarkerSize',12)
plot(tau, V(start:finish, 4+ (fig-1)*num_per_dgp),'-.','LineWidth',3, 'color', 'blue','MarkerSize',12)
plot(tau, V(start:finish, 5+ (fig-1)*num_per_dgp),'--','LineWidth',3, 'color', [0 0.5 0],'MarkerSize',12)
plot(tau, V(start:finish,6+ (fig-1)*num_per_dgp),'-.o','LineWidth',3, 'color', [0 0.5 0],'MarkerSize',12)
plot(tau, V(start:finish,7+ (fig-1)*num_per_dgp),'-.*','LineWidth',3, 'color', [0 0.5 0],'MarkerSize',12)
plot(tau, V(start:finish,8+ (fig-1)*num_per_dgp),'-.x','LineWidth',3, 'color', 'red','MarkerSize',12)

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
saveas(gcf,append('PowerCurves/CovariateSelection/','power_sel', int2str(H(j)), int2str(K(k)), side ,int2str(pub_b)), 'epsc')
close all
end
if (minimum == 1)
saveas(gcf,append('PowerCurves/CovariateSelection/','power_sel', int2str(H(j)), int2str(K(k)), 'min', side,int2str(pub_b)), 'epsc')
close all
end
end
end

for j = 1:4
for k = 1:2
    fig=fig+1;

figure (fig)
plot(tau, V(start:finish, 2+(fig-1)*num_per_dgp),'LineWidth',3,'color' , 'black')
hold on
%plot(tau, V(:,3+ (fig-1)*9),'-x','LineWidth',3, 'color', 'black','MarkerSize',12)
plot(tau, V(start:finish, 4+ (fig-1)*num_per_dgp),'-.','LineWidth',3, 'color', 'blue','MarkerSize',12)
plot(tau, V(start:finish, 5+ (fig-1)*num_per_dgp),'--','LineWidth',3, 'color', [0 0.5 0],'MarkerSize',12)
plot(tau, V(start:finish,6+ (fig-1)*num_per_dgp),'-.o','LineWidth',3, 'color', [0 0.5 0],'MarkerSize',12)
plot(tau, V(start:finish,7+ (fig-1)*num_per_dgp),'-.*','LineWidth',3, 'color', [0 0.5 0],'MarkerSize',12)
plot(tau, V(start:finish,8+ (fig-1)*num_per_dgp),'-.x','LineWidth',3, 'color', 'red','MarkerSize',12)

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
saveas(gcf,append('PowerCurves/IVSelection/','power_iv', int2str(H(j)), int2str(K(k)), side,int2str(pub_b)), 'epsc')
close all
end
if (minimum==1)
saveas(gcf,append('PowerCurves/IVSelection/','power_iv', int2str(H(j)), int2str(K(k)), 'min', side,int2str(pub_b)), 'epsc')
close all
end
end
end

for j = 1:4
for k = 1:1
    fig=fig+1;

figure (fig)
plot(tau, V(start:finish, 2+(fig-1)*num_per_dgp),'LineWidth',3,'color' , 'black')
hold on
%plot(tau, V(:,3+ (fig-1)*9),'-x','LineWidth',3, 'color', 'black','MarkerSize',12)
plot(tau, V(start:finish, 4+ (fig-1)*num_per_dgp),'-.','LineWidth',3, 'color', 'blue','MarkerSize',12)
plot(tau, V(start:finish, 5+ (fig-1)*num_per_dgp),'--','LineWidth',3, 'color', [0 0.5 0],'MarkerSize',12)
plot(tau, V(start:finish,6+ (fig-1)*num_per_dgp),'-.o','LineWidth',3, 'color', [0 0.5 0],'MarkerSize',12)
plot(tau, V(start:finish,7+ (fig-1)*num_per_dgp),'-.*','LineWidth',3, 'color', [0 0.5 0],'MarkerSize',12)
plot(tau, V(start:finish,8+ (fig-1)*num_per_dgp),'-.x','LineWidth',3, 'color', 'red','MarkerSize',12)

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
saveas(gcf,append('PowerCurves/LagLengthSelection/','power_var', int2str(H(j)), side,int2str(pub_b)), 'epsc')
close all
end
if (minimum == 1)
saveas(gcf,append('PowerCurves/LagLengthSelection/','power_var', int2str(H(j)), 'min', side,int2str(pub_b)), 'epsc')
close all
end
end
end
end
end
end