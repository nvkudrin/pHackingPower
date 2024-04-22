clear all

S = readtable("PowerCurves/RejectionRates_March18.csv");
Vfull=table2array(S);
Vfull = Vfull(:,2:end);


tau = [0:0.05:1];
H = [0, 1, 2, 3];
K = [3,5,7];
num_per_dgp = 11;


% Tables
A = Vfull([11;32;53], :);
column_names = ["Binomial", "Discontinuity", "CS1", "CSUB", "CS2B", "LCM"];
row_names = ["No Pub Bias", "Sharp Pub Bias", "Smooth Pub Bias"];
row_names = [" ", "Cov Selection (K = 3,thresholding)", row_names, "Cov Selection (K = 3, minimum)", row_names,...
              "IV Selection (K = 3,thresholding)", row_names, "IV Selection (K = 3, minimum)", row_names,...
              "Lag Selection (thresholding)", row_names, "Lag Selection (minimum)", row_names,...
              "Cluster Selection (thresholding)", row_names, "Cluster Selection (minimum)", row_names];
    M0 = [10 50 25 65 41 81 45 85];
    M1 = [14 54 27 67 42 82 46 86];
    M2 = [18 58 29 69 43 83 47 87];
    M3 = [22 62 31 71 44 84 48 88];
    Mall = [M0;M1;M2;M3];

for h = H+1
  Tablel = NaN(32,6);
  row_num = 1;
for j = Mall(h,:)
  sub_table = A(:,((j-1)*num_per_dgp+1):(j*num_per_dgp));
  sub_table(:,[1,3,9,10,11])=[];
  Tablel((row_num+1):(row_num+3),:) = sub_table;
  %Tablel((row_num+5):(row_num+7),:) = sub_table_min;
  row_num = row_num+4;
end
Table_final = [row_names', [column_names; round(Tablel,3)]];
writematrix(Table_final,append('PowerCurves/Table_', '2', 'sided_',int2str(h),'.csv'))
end
%% Figures
fig = 0;
for sided = [1,2]
    %V = Vfull(:, (1+(sided-1)*num_col/2): (num_col/(3-sided)));
    V = Vfull;


    if sided == 1
        side = "1sided";
    end
    if sided == 2
        side = "2sided";
    end
%for pub_b = [0,1,2]
for pub_b = [0]
    %fig = 0;
    start = pub_b*21+1;
    finish = 21*(pub_b+1);
for minimum = 0:1
for j = 1:4
    for g2s = 0:1
for k = 1:3
    if (((g2s==1)&&(sided == 2))||((g2s==0)&&(k==1)&&(sided==2))||((g2s==1)&&(k==1)&&(sided==1)))
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

if (k == 1)
if (j==1)&&(minimum==0)
lgd = legend('Binomial', 'Discontinuity',...
'CS1','CSUB', 'CS2B','LCM', 'Location','northwest');
lgd.FontSize = 16;
end
else
if (j==1)&&(minimum==1)
lgd = legend('Binomial', 'Discontinuity',...
'CS1','CSUB', 'CS2B','LCM', 'Location','northwest');
lgd.FontSize = 16;
end
end


set(gca,'FontSize',18)
xlabel('$$\tau$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('Rejection Rate', 'FontSize',25, 'interpreter', 'latex')
ylim([0 1])
xlim([0 1])
title(append('$K = $ ', int2str(K(k)), ', $h = $ ', int2str(H(j))),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')


if (j==4)
    title(append('$K = $ ', int2str(K(k)), ', $h \sim \widehat \Pi$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end
if (minimum == 0)
saveas(gcf,append('PowerCurves/CovariateSelection/','power_sel', int2str(H(j)), int2str(K(k)), side ,int2str(pub_b),int2str(g2s)), 'epsc')
close all
end
if (minimum == 1)
saveas(gcf,append('PowerCurves/CovariateSelection/','power_sel', int2str(H(j)), int2str(K(k)), 'min', side,int2str(pub_b),int2str(g2s)), 'epsc')
close all
end
end
    end
end
end
%%IV
for j = 1:4
    for g2s = 0:1
for k = 1:2
    if ((sided == 2)&&(g2s==1))

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

if (k == 1)
if (j==1)&&(minimum==0)
lgd = legend('Binomial', 'Discontinuity',...
'CS1','CSUB', 'CS2B','LCM', 'Location','northwest');
lgd.FontSize = 16;
end
else
if (j==1)&&(minimum==1)
lgd = legend('Binomial', 'Discontinuity',...
'CS1','CSUB', 'CS2B','LCM', 'Location','northwest');
lgd.FontSize = 16;
end
end


set(gca,'FontSize',18)
xlabel('$$\tau$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('Rejection Rate', 'FontSize',25, 'interpreter', 'latex')
ylim([0 1])
xlim([0 1])
title(append('$K = $ ', int2str(K(k)), ', $h = $ ', int2str(H(j))),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')

if (j==4)
    title(append('$K = $ ', int2str(K(k)), ', $h \sim \widehat \Pi$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end

if (minimum==0)
saveas(gcf,append('PowerCurves/IVSelection/','power_iv', int2str(H(j)), int2str(K(k)), side,int2str(pub_b),int2str(g2s)), 'epsc')
close all
end
if (minimum==1)
saveas(gcf,append('PowerCurves/IVSelection/','power_iv', int2str(H(j)), int2str(K(k)), 'min', side,int2str(pub_b),int2str(g2s)), 'epsc')
close all
end
end
    end
end
end
%%IV_F
for j = 1:4
    for g2s=0:1
for k = 1:2
    if ((sided == 2)&&(g2s==1))
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

if (k == 1)
if (j==1)&&(minimum==0)
lgd = legend('Binomial', 'Discontinuity',...
'CS1','CSUB', 'CS2B','LCM', 'Location','northwest');
lgd.FontSize = 16;
end
else
if (j==1)&&(minimum==1)
lgd = legend('Binomial', 'Discontinuity',...
'CS1','CSUB', 'CS2B','LCM', 'Location','northwest');
lgd.FontSize = 16;
end
end


set(gca,'FontSize',18)
xlabel('$$\tau$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('Rejection Rate', 'FontSize',25, 'interpreter', 'latex')
ylim([0 1])
xlim([0 1])
title(append('$K = $ ', int2str(K(k)), ', $h = $ ', int2str(H(j))),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')

if (j==4)
    title(append('$K = $ ', int2str(K(k)), ', $h \sim \widehat \Pi$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end

if (minimum==0)
saveas(gcf,append('PowerCurves/IV_FSelection/','power_iv_F', int2str(H(j)), int2str(K(k)), side,int2str(pub_b),int2str(g2s)), 'epsc')
close all
end
if (minimum==1)
saveas(gcf,append('PowerCurves/IV_FSelection/','power_iv_F', int2str(H(j)), int2str(K(k)), 'min', side,int2str(pub_b),int2str(g2s)), 'epsc')
close all
end
end
    end
end
end

%%LagLength

for j = 1:4
for g2s=0:1
for k = 1:1
    if ((sided == 2)&&(g2s==1))
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

if (j==1)&&(minimum==0)
lgd = legend('Binomial', 'Discontinuity',...
'CS1','CSUB', 'CS2B','LCM', 'Location','northwest');
lgd.FontSize = 16;
end

set(gca,'FontSize',18)
xlabel('$$\tau$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('Rejection Rate', 'FontSize',25, 'interpreter', 'latex')
ylim([0 1])
xlim([0 1])
title(append('$h = $ ', int2str(H(j))),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')

if (j==4)
    title(append('$h \sim \widehat \Pi$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end

if (minimum == 0)
saveas(gcf,append('PowerCurves/LagLengthSelection/','power_var', int2str(H(j)), side,int2str(pub_b),int2str(g2s)), 'epsc')
close all
end
if (minimum == 1)
saveas(gcf,append('PowerCurves/LagLengthSelection/','power_var', int2str(H(j)), 'min', side,int2str(pub_b),int2str(g2s)), 'epsc')
close all
end
end
end
end
end

%Clustering
for j = 1:4
    for g2s=0:1
for k = 1:1
    if ((sided == 2)&&(g2s==1))
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

if (j==1)&&(minimum==0)
lgd = legend('Binomial', 'Discontinuity',...
'CS1','CSUB', 'CS2B','LCM', 'Location','northwest');
lgd.FontSize = 16;
end


set(gca,'FontSize',18)
xlabel('$$\tau$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('Rejection Rate', 'FontSize',25, 'interpreter', 'latex')
ylim([0 1])
xlim([0 1])
title(append('$h = $ ', int2str(H(j))),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')

if (j==4)
    title(append('$h \sim \widehat \Pi$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end

if (minimum == 0)
saveas(gcf,append('PowerCurves/ClusterSelection/','power_clust', int2str(H(j)), side,int2str(pub_b),int2str(g2s)), 'epsc')
close all
end
if (minimum == 1)
saveas(gcf,append('PowerCurves/ClusterSelection/','power_clust', int2str(H(j)), 'min', side,int2str(pub_b),int2str(g2s)), 'epsc')
close all
end
end
    end
end
end
end
end
end