clear all
load('DGPs/MC_distributions_March18.mat') %use the file generated by DataGeneration.m
rng(12345)
tstat = readmatrix("Brodeur_data.csv",'NumHeaderLines',1);
threshold = 10;
data = tstat(abs(tstat)<threshold);
fun = @(theta) -PiGamma(data, theta(1), theta(2));
res = fmincon(fun, [1, 1]);
h_brodeur_RCT = gamrnd(res(1), res(2), 1, 10^6);

H = [0, 1, 2, 3, 4];
K = [3,5,7];
Up = [0,1];
nobs = 200;

BIAS = struct;
BIAS.var = zeros(2,5,2);
BIAS.sel = zeros(2, 12,2);
BIAS.iv = zeros(2, 8, 2);
%%
fig = 0;
for s = 1:2
b_ind_sel = 0;
b_ind_iv = 0;
for j = 1:4
    close all
    if (j<4)
b = (H(j)/sqrt(nobs))*ones(size(mcout_bias));
b_iv3 = (H(j)/sqrt(nobs)/3)*ones(size(mcout_bias_iv3));
b_iv5 = (H(j)/sqrt(nobs)/3)*ones(size(mcout_bias_iv5));
b_var = (H(j)/sqrt(nobs))*ones(size(mcout_bias_var));
b_clust = (H(j)/sqrt(nobs))*ones(size(mcout_bias_cluster));
    end


    if (j==4)
b = repmat(h_brodeur_RCT, 128, 1)/sqrt(nobs);
b_iv3 = repmat(h_brodeur_RCT, 7, 1)/sqrt(nobs)/3;
b_iv5 = repmat(h_brodeur_RCT, 31, 1)/sqrt(nobs)/3;
b_var = repmat(h_brodeur_RCT, 5, 1)/sqrt(nobs);
b_clust = repmat(h_brodeur_RCT, 5, 1)/sqrt(nobs);
    end

for GeneralToSpecific = 0:1   
 %LagLengthSelection   
    [P0, P1, P1min] = NullAndAlt_var_bic(b_var, mcout_bias_var, mcout_se_var, s, 0.05, mcout_bic_var);
csvwrite(append('DGPs/LagLengthSelection/','DGPs','P0var', int2str(H(j)),int2str(s),'sided', int2str(GeneralToSpecific), '.csv'), P0)
% fig = fig+1;
% figure(fig)
% histogram(P0, 100, 'Normalization', 'probability')
% set(gca,'FontSize',18)
% xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
% ylabel('Probability', 'FontSize',25, 'interpreter', 'latex')
% if (j<4)
% title(append('No $p$-hacking: $h = $ ', int2str(H(j))),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% if (j>4)
%     title(append('No $p$-hacking: $h \sim \chi^2(1)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% if (j==4)
%     title(append('No $p$-hacking: $h \sim \widehat\Pi(h)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% saveas(gcf,append('DGPs/LagLengthSelection/','P0var', int2str(H(j)),int2str(s),'sided', int2str(GeneralToSpecific)), 'epsc')

csvwrite(append('DGPs/LagLengthSelection/','DGPs','P1var', int2str(H(j)),int2str(s),'sided', int2str(GeneralToSpecific), '.csv'), P1)

% fig = fig+1;
% figure(fig)
% histogram(P1, 100, 'Normalization', 'probability')
% set(gca,'FontSize',18)
% xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
% ylabel('Probability', 'FontSize',25, 'interpreter', 'latex')
% if (j<4)
% title(append('p-hacked (threshold): $h = $ ', int2str(H(j))),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% if (j>4)
%     title(append('p-hacked (threshold): $h \sim \chi^2(1)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% if (j==4)
%     title(append('p-hacked (threshold): $h \sim \widehat\Pi(h)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% saveas(gcf,append('DGPs/LagLengthSelection/','P1var', int2str(H(j)),int2str(s),'sided', int2str(GeneralToSpecific)), 'epsc')

% fig = fig+1;
% figure(fig)
% histogram(P1min, 100, 'Normalization', 'probability')
% set(gca,'FontSize',18)
% xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
% ylabel('Probability', 'FontSize',25, 'interpreter', 'latex')
% if (j<4)
% title(append('p-hacked (minimum): $h = $ ', int2str(H(j))),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% if (j>4)
%     title(append('p-hacked (minimum): $h \sim \chi^2(1)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% if (j==4)
%     title(append('p-hacked (minimum): $h \sim \widehat\Pi(h)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% saveas(gcf,append('DGPs/LagLengthSelection/','P1var', int2str(H(j)), 'min',int2str(s),'sided', int2str(GeneralToSpecific)), 'epsc')
csvwrite(append('DGPs/LagLengthSelection/','DGPs','P1var', int2str(H(j)), 'min',int2str(s),'sided', int2str(GeneralToSpecific), '.csv'), P1min) 
[edges0, N0] = EdgesN(P0);
[edges1, N1] = EdgesN(P1);
[edges1min, N1min] = EdgesN(P1min);
fig = fig+1;
figure(fig)
hold on
plot(edges0, N0, 'LineStyle','-', 'LineWidth',2, 'Color', 'black')
plot(edges1, N1, 'LineStyle','--', 'LineWidth',2, 'Color', 'black')
plot(edges1min, N1min, 'LineStyle',':', 'LineWidth',2, 'Color', 'black')
if j==1
    ylim([0 10])
end
set(gca,'FontSize',18)
xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('PDF', 'FontSize',25, 'interpreter', 'latex')
if (j<4)
title(append('Null and $p$-hacked distributions: $h = $ ', int2str(H(j))),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end
if (j>4)
    title(append('Null and $p$-hacked distributions: $h \sim \chi^2(1)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end
if (j==4)
    title(append('Null and $p$-hacked distributions: $h \sim \widehat\Pi(h)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end
lgd = legend('No $p$-hacking','Threshold', 'Minimum', 'Orientation','horizontal', 'interpreter', 'latex', 'NumColumns',3);
axes('position',[.35 .34 .45 .45])
box on 
indexOfInterest = (edges0 < 0.2) & (edges0 > 0.0001); 
plot(edges0(indexOfInterest),N0(indexOfInterest), '-' , 'LineWidth', 2, 'color', 'black') % plot on new axes
hold on
plot(edges0(indexOfInterest),N1(indexOfInterest), '--' , 'LineWidth', 2, 'color', 'black')
plot(edges0(indexOfInterest),N1min(indexOfInterest), ':' , 'LineWidth', 2, 'color', 'black')
%ylim([0 40])
xlim([0 0.1])
saveas(gcf,append('DGPs/LagLengthSelection/','Pvar', int2str(H(j)), 'all',int2str(s),'sided', int2str(GeneralToSpecific)), 'epsc')
%%
 %Clustering   
    [P0, P1, P1min] = NullAndAlt_var_clust(b_clust, mcout_bias_cluster, mcout_se_cluster, s, 0.05, GeneralToSpecific);
csvwrite(append('DGPs/ClusterSelection/','DGPs','P0clust', int2str(H(j)),int2str(s),'sided', int2str(GeneralToSpecific), '.csv'), P0)
% fig = fig+1;
% figure(fig)
% histogram(P0, 100, 'Normalization', 'probability')
% set(gca,'FontSize',18)
% xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
% ylabel('Probability', 'FontSize',25, 'interpreter', 'latex')
% if (j<4)
% title(append('No $p$-hacking: $h = $ ', int2str(H(j))),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% if (j==4)
%     title(append('No $p$-hacking: $h \sim \chi^2(1)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% if (j==5)
%     title(append('No $p$-hacking: $h \sim \widehat\Pi(h)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% saveas(gcf,append('DGPs/ClusterSelection/','P0clust', int2str(H(j)),int2str(s),'sided', int2str(GeneralToSpecific)), 'epsc')

csvwrite(append('DGPs/ClusterSelection/','DGPs','P1clust', int2str(H(j)),int2str(s),'sided', int2str(GeneralToSpecific), '.csv'), P1)

% fig = fig+1;
% figure(fig)
% histogram(P1, 100, 'Normalization', 'probability')
% set(gca,'FontSize',18)
% xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
% ylabel('Probability', 'FontSize',25, 'interpreter', 'latex')
% if (j<4)
% title(append('p-hacked (threshold): $h = $ ', int2str(H(j))),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% if (j>4)
%     title(append('p-hacked (threshold): $h \sim \chi^2(1)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% if (j==4)
%     title(append('p-hacked (threshold): $h \sim \widehat\Pi(h)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% saveas(gcf,append('DGPs/ClusterSelection/','P1clust', int2str(H(j)),int2str(s),'sided', int2str(GeneralToSpecific)), 'epsc')
% 
% fig = fig+1;
% figure(fig)
% histogram(P1min, 100, 'Normalization', 'probability')
% set(gca,'FontSize',18)
% xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
% ylabel('Probability', 'FontSize',25, 'interpreter', 'latex')
% if (j<4)
% title(append('p-hacked (minimum): $h = $ ', int2str(H(j))),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% if (j==5)
%     title(append('p-hacked (minimum): $h \sim \chi^2(1)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% if (j==4)
%     title(append('p-hacked (minimum): $h \sim \widehat\Pi(h)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% saveas(gcf,append('DGPs/ClusterSelection/','P1clust', int2str(H(j)), 'min',int2str(s),'sided', int2str(GeneralToSpecific)), 'epsc')
csvwrite(append('DGPs/ClusterSelection/','DGPs','P1clust', int2str(H(j)), 'min',int2str(s),'sided', int2str(GeneralToSpecific), '.csv'), P1min) 
[edges0, N0] = EdgesN(P0);
[edges1, N1] = EdgesN(P1);
[edges1min, N1min] = EdgesN(P1min);
fig = fig+1;
figure(fig)
hold on
plot(edges0, N0, 'LineStyle','-', 'LineWidth',2, 'Color', 'black')
plot(edges1, N1, 'LineStyle','--', 'LineWidth',2, 'Color', 'black')
plot(edges1min, N1min, 'LineStyle',':', 'LineWidth',2, 'Color', 'black')
if j==1
    ylim([0 10])
end
set(gca,'FontSize',18)
xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('PDF', 'FontSize',25, 'interpreter', 'latex')
if (j<4)
title(append('Null and $p$-hacked distributions: $h = $ ', int2str(H(j))),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end
if (j==5)
    title(append('Null and $p$-hacked distributions: $h \sim \chi^2(1)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end
if (j==4)
    title(append('Null and $p$-hacked distributions: $h \sim \widehat\Pi(h)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end
lgd = legend('No $p$-hacking','Threshold', 'Minimum', 'Orientation','horizontal', 'interpreter', 'latex', 'NumColumns',3);
axes('position',[.35 .34 .45 .45])
box on 
indexOfInterest = (edges0 < 0.2) & (edges0 > 0.0001); 
plot(edges0(indexOfInterest),N0(indexOfInterest), '-' , 'LineWidth', 2, 'color', 'black') % plot on new axes
hold on
plot(edges0(indexOfInterest),N1(indexOfInterest), '--' , 'LineWidth', 2, 'color', 'black')
plot(edges0(indexOfInterest),N1min(indexOfInterest), ':' , 'LineWidth', 2, 'color', 'black')
%ylim([0 40])
xlim([0 0.1])
saveas(gcf,append('DGPs/ClusterSelection/','Pclust', int2str(H(j)), 'all',int2str(s),'sided', int2str(GeneralToSpecific)), 'epsc')
%%
for k = 1:3
    %fig
%CovariateSelection
[P0, P1,P1min, B0, B1, B1min] = NullAndAlt(b, mcout_bias, mcout_se, s, 0.05, 7, K(k), GeneralToSpecific);
if ((s==2)&&(GeneralToSpecific==1))
b_ind_sel = b_ind_sel+1;
BIAS.sel(1, b_ind_sel,s) = mean(B1);
BIAS.sel(2, b_ind_sel,s) = mean(B1min);
end

csvwrite(append('DGPs/CovariateSelection/','DGPs','P0sel', int2str(H(j)), int2str(K(k)),int2str(s),'sided', int2str(GeneralToSpecific), '.csv'), P0)
csvwrite(append('DGPs/CovariateSelection/','DGPs','P1sel', int2str(H(j)), int2str(K(k)),int2str(s),'sided', int2str(GeneralToSpecific), '.csv'), P1)
csvwrite(append('DGPs/CovariateSelection/','DGPs','P1sel', int2str(H(j)), int2str(K(k)),'min',int2str(s),'sided', int2str(GeneralToSpecific), '.csv'), P1min)
% fig = fig+1;
% figure(fig)
% histogram(P0, 100, 'Normalization', 'probability')
% set(gca,'FontSize',18)
% xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
% ylabel('Probability', 'FontSize',25, 'interpreter', 'latex')
% if (j<4)
% title(append('No $p$-hacking: $K = $ ', int2str(K(k)), ', $h = $ ', int2str(H(j))),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% if (j==5)
%     title(append('No $p$-hacking: $K = $', int2str(K(k)), ', $h \sim \chi^2(1)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% if (j==4)
%     title(append('No $p$-hacking: $K = $', int2str(K(k)), ', $h \sim \widehat\Pi(h)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% saveas(gcf,append('DGPs/CovariateSelection/','P0sel', int2str(H(j)), int2str(K(k)),int2str(s),'sided', int2str(GeneralToSpecific)), 'epsc')
% 
% fig = fig+1;
% figure(fig)
% histogram(P1, 100, 'Normalization', 'probability')
% set(gca,'FontSize',18)
% xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
% ylabel('Probability', 'FontSize',25, 'interpreter', 'latex')
% if (j<4)
% title(append('p-hacked (threshold): $K = $ ', int2str(K(k)), ', $h = $ ', int2str(H(j))),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% if (j==5)
%     title(append('p-hacked (threshold): $K = $ ', int2str(K(k)), ', $h \sim \chi^2(1)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% if (j==4)
%     title(append('p-hacked (threshold): $K = $ ', int2str(K(k)), ', $h \sim \widehat\Pi(h)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% saveas(gcf,append('DGPs/CovariateSelection/','P1sel', int2str(H(j)), int2str(K(k)),int2str(s),'sided', int2str(GeneralToSpecific)), 'epsc')
% 
% fig = fig+1;
% figure(fig)
% histogram(P1min, 100, 'Normalization', 'probability')
% set(gca,'FontSize',18)
% xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
% ylabel('Probability', 'FontSize',25)
% if (j<4)
% title(append('p-hacked (minimum): $K = $ ', int2str(K(k)), ', $h = $ ', int2str(H(j))),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% if (j==5)
%     title(append('p-hacked (minimum): $K = $ ', int2str(K(k)), ', $h \sim \chi^2(1)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% if (j==4)
%     title(append('p-hacked (minimum): $K = $ ', int2str(K(k)), ', $h \sim \widehat\Pi(h)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% saveas(gcf,append('DGPs/CovariateSelection/','P1sel', int2str(H(j)), int2str(K(k)), 'min',int2str(s),'sided', int2str(GeneralToSpecific)), 'epsc')

[edges0, N0] = EdgesN(P0);
[edges1, N1] = EdgesN(P1);
[edges1min, N1min] = EdgesN(P1min);
fig = fig+1;
figure(fig)
hold on
plot(edges0, N0, 'LineStyle','-', 'LineWidth',2, 'Color', 'black')
plot(edges1, N1, 'LineStyle','--', 'LineWidth',2, 'Color', 'black')
plot(edges1min, N1min, 'LineStyle',':', 'LineWidth',2, 'Color', 'black')
if j==1
    ylim([0 10])
end
set(gca,'FontSize',18)
xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('PDF', 'FontSize',25, 'interpreter', 'latex')
if (j<4)
title(append('Null and $p$-hacked distributions: $h = $ ', int2str(H(j))),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end
if (j==5)
    title(append('Null and $p$-hacked distributions: $h \sim \chi^2(1)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end
if (j==4)
    title(append('Null and $p$-hacked distributions: $h \sim \widehat\Pi(h)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end
lgd = legend('No $p$-hacking','Threshold', 'Minimum', 'Orientation','horizontal', 'interpreter', 'latex', 'NumColumns',3);
axes('position',[.35 .34 .45 .45])
box on 
indexOfInterest = (edges0 < 0.2) & (edges0 > 0.0001); 
plot(edges0(indexOfInterest),N0(indexOfInterest), '-' , 'LineWidth', 2, 'color', 'black') % plot on new axes
hold on
plot(edges0(indexOfInterest),N1(indexOfInterest), '--' , 'LineWidth', 2, 'color', 'black')
plot(edges0(indexOfInterest),N1min(indexOfInterest), ':' , 'LineWidth', 2, 'color', 'black')
%ylim([0 40])
xlim([0 0.1])
saveas(gcf,append('DGPs/CovariateSelection/','Psel', int2str(H(j)), int2str(K(k)), 'all',int2str(s),'sided', int2str(GeneralToSpecific)), 'epsc')

%IVSelection
if (k<3)
    %%
    if (k==1)
[P0, P1,P1min, B0, B1, B1min] = NullAndAlt(b_iv3, mcout_bias_iv3, mcout_se_iv3, s, 0.05, 3, 3,GeneralToSpecific);
    end
    if (k==2)
[P0, P1,P1min, B0, B1, B1min] = NullAndAlt(b_iv5, mcout_bias_iv5, mcout_se_iv5, s, 0.05, 5, 5, GeneralToSpecific);  
    end
    if ((s==2)&&(GeneralToSpecific==1))&&(k<3)
        b_ind_iv = b_ind_iv+1; 
BIAS.iv(1, b_ind_iv,s) = mean(B1);
BIAS.iv(2, b_ind_iv,s) = mean(B1min);
    end
csvwrite(append('DGPs/IVSelection/','DGPs','P0iv', int2str(H(j)), int2str(K(k)),int2str(s),'sided', int2str(GeneralToSpecific), '.csv'), P0)
csvwrite(append('DGPs/IVSelection/','DGPs','P1iv', int2str(H(j)), int2str(K(k)),int2str(s),'sided', int2str(GeneralToSpecific), '.csv'), P1)
csvwrite(append('DGPs/IVSelection/','DGPs','P1iv', int2str(H(j)), int2str(K(k)),'min',int2str(s),'sided', int2str(GeneralToSpecific), '.csv'), P1min)
% fig = fig+1;
% figure(fig)
% histogram(P0, 100, 'Normalization', 'probability')
% set(gca,'FontSize',18)
% xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
% ylabel('Probability', 'FontSize',25, 'interpreter', 'latex')
% if (j<4)
% title(append('No $p$-hacking: $K = $ ', int2str(K(k)), ', $h = $ ', int2str(H(j))),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% if (j==5)
%     title(append('No $p$-hacking: $K = $ ', int2str(K(k)), ', $h \sim \chi^2(1)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% if (j==4)
%     title(append('No $p$-hacking: $K = $ ', int2str(K(k)), ', $h \sim \widehat\Pi(h)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% saveas(gcf,append('DGPs/IVSelection/','P0iv', int2str(H(j)), int2str(K(k)),int2str(s),'sided', int2str(GeneralToSpecific)), 'epsc')
% 
% 
% fig = fig+1;
% figure(fig)
% histogram(P1, 100, 'Normalization', 'probability')
% set(gca,'FontSize',18)
% xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
% ylabel('Probability', 'FontSize',25, 'interpreter', 'latex')
% if (j<4)
% title(append('p-hacked (threshold): $K = $ ', int2str(K(k)), ', $h = $ ', int2str(H(j))),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% if (j==5)
%     title(append('p-hacked (threshold): $K = $ ', int2str(K(k)), ', $h \sim \chi^2(1)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% if (j==4)
%     title(append('p-hacked (threshold): $K = $ ', int2str(K(k)), ', $h \sim \widehat\Pi(h)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% saveas(gcf,append('DGPs/IVSelection/','P1iv', int2str(H(j)), int2str(K(k)),int2str(s),'sided', int2str(GeneralToSpecific)), 'epsc')
% 
% 
% 
% fig = fig+1;
% figure(fig)
% histogram(P1min, 100, 'Normalization', 'probability')
% set(gca,'FontSize',18)
% xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
% ylabel('Probability', 'FontSize',25, 'interpreter', 'latex')
% if (j<4)
% title(append('p-hacked (minimum): $K = $ ', int2str(K(k)), ', $h = $ ', int2str(H(j))),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% if (j==5)
%     title(append('p-hacked (minimum): $K = $ ', int2str(K(k)), ', $h \sim \chi^2(1)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% if (j==4)
%     title(append('p-hacked (minimum): $K = $ ', int2str(K(k)), ', $h \sim \widehat\Pi(h)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% saveas(gcf,append('DGPs/IVSelection/','P1iv', int2str(H(j)), int2str(K(k)), 'min',int2str(s),'sided', int2str(GeneralToSpecific)), 'epsc')

[edges0, N0] = EdgesN(P0);
[edges1, N1] = EdgesN(P1);
[edges1min, N1min] = EdgesN(P1min);
fig = fig+1;
figure(fig)
hold on
plot(edges0, N0, 'LineStyle','-', 'LineWidth',2, 'Color', 'black')
plot(edges1, N1, 'LineStyle','--', 'LineWidth',2, 'Color', 'black')
plot(edges1min, N1min, 'LineStyle',':', 'LineWidth',2, 'Color', 'black')
if j==1
    ylim([0 10])
end
set(gca,'FontSize',18)
xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('PDF', 'FontSize',25, 'interpreter', 'latex')
if (j<4)
title(append('Null and $p$-hacked distributions: $h = $ ', int2str(H(j))),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end
if (j==5)
    title(append('Null and $p$-hacked distributions: $h \sim \chi^2(1)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end
if (j==4)
    title(append('Null and $p$-hacked distributions: $h \sim \widehat\Pi(h)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end
lgd = legend('No $p$-hacking','Threshold', 'Minimum', 'Orientation','horizontal', 'interpreter', 'latex', 'NumColumns',3);
axes('position',[.35 .34 .45 .45])
box on 
indexOfInterest = (edges0 < 0.2) & (edges0 > 0.0001); 
plot(edges0(indexOfInterest),N0(indexOfInterest), '-' , 'LineWidth', 2, 'color', 'black') % plot on new axes
hold on
plot(edges0(indexOfInterest),N1(indexOfInterest), '--' , 'LineWidth', 2, 'color', 'black')
plot(edges0(indexOfInterest),N1min(indexOfInterest), ':' , 'LineWidth', 2, 'color', 'black')
%ylim([0 40])
xlim([0 0.1])
saveas(gcf,append('DGPs/IVSelection/','Piv', int2str(H(j)), int2str(K(k)), 'all',int2str(s),'sided', int2str(GeneralToSpecific)), 'epsc')

%%
    %%
    if (k==1)
[P0, P1,P1min, B0, B1, B1min] = NullAndAlt_iv_Fstat(b_iv3, mcout_bias_iv3, mcout_se_iv3, s, 0.05, 3, 3,GeneralToSpecific, mcout_F_iv3);
%b_ind_iv = b_ind_iv+1;
    end
    if (k==2)
[P0, P1,P1min, B0, B1, B1min] = NullAndAlt_iv_Fstat(b_iv5, mcout_bias_iv5, mcout_se_iv5, s, 0.05, 5, 5, GeneralToSpecific,mcout_F_iv5);
%b_ind_iv = b_ind_iv+1;  
    end
%BIAS.iv(1, b_ind_iv,s) = mean(B1);
%BIAS.iv(2, b_ind_iv,s) = mean(B1min);
csvwrite(append('DGPs/IVSelection/','DGPs','P0iv_F', int2str(H(j)), int2str(K(k)),int2str(s),'sided', int2str(GeneralToSpecific), '.csv'), P0)
csvwrite(append('DGPs/IVSelection/','DGPs','P1iv_F', int2str(H(j)), int2str(K(k)),int2str(s),'sided', int2str(GeneralToSpecific), '.csv'), P1)
csvwrite(append('DGPs/IVSelection/','DGPs','P1iv_F', int2str(H(j)), int2str(K(k)),'min',int2str(s),'sided', int2str(GeneralToSpecific), '.csv'), P1min)

% fig = fig+1;
% figure(fig)
% histogram(P0, 100, 'Normalization', 'probability')
% set(gca,'FontSize',18)
% xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
% ylabel('Probability', 'FontSize',25, 'interpreter', 'latex')
% if (j<4)
% title(append('No $p$-hacking: $K = $ ', int2str(K(k)), ', $h = $ ', int2str(H(j))),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% if (j==5)
%     title(append('No $p$-hacking: $K = $ ', int2str(K(k)), ', $h \sim \chi^2(1)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% if (j==4)
%     title(append('No $p$-hacking: $K = $ ', int2str(K(k)), ', $h \sim \widehat\Pi(h)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% saveas(gcf,append('DGPs/IVSelection/','P0ivF', int2str(H(j)), int2str(K(k)),int2str(s),'sided', int2str(GeneralToSpecific)), 'epsc')
% 
% 
% fig = fig+1;
% figure(fig)
% histogram(P1, 100, 'Normalization', 'probability')
% set(gca,'FontSize',18)
% xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
% ylabel('Probability', 'FontSize',25, 'interpreter', 'latex')
% if (j<4)
% title(append('p-hacked (threshold): $K = $ ', int2str(K(k)), ', $h = $ ', int2str(H(j))),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% if (j==5)
%     title(append('p-hacked (threshold): $K = $ ', int2str(K(k)), ', $h \sim \chi^2(1)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% if (j==4)
%     title(append('p-hacked (threshold): $K = $ ', int2str(K(k)), ', $h \sim \widehat\Pi(h)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% saveas(gcf,append('DGPs/IVSelection/','P1ivF', int2str(H(j)), int2str(K(k)),int2str(s),'sided', int2str(GeneralToSpecific)), 'epsc')
% 
% 
% fig = fig+1;
% figure(fig)
% histogram(P1min, 100, 'Normalization', 'probability')
% set(gca,'FontSize',18)
% xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
% ylabel('Probability', 'FontSize',25, 'interpreter', 'latex')
% if (j<4)
% title(append('p-hacked (minimum): $K = $ ', int2str(K(k)), ', $h = $ ', int2str(H(j))),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% if (j==5)
%     title(append('p-hacked (minimum): $K = $ ', int2str(K(k)), ', $h \sim \chi^2(1)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% if (j==4)
%     title(append('p-hacked (minimum): $K = $ ', int2str(K(k)), ', $h \sim \widehat\Pi(h)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
% end
% saveas(gcf,append('DGPs/IVSelection/','P1ivF', int2str(H(j)), int2str(K(k)), 'min',int2str(s),'sided', int2str(GeneralToSpecific)), 'epsc')

[edges0, N0] = EdgesN(P0);
[edges1, N1] = EdgesN(P1);
[edges1min, N1min] = EdgesN(P1min);
fig = fig+1;
figure(fig)
hold on
plot(edges0, N0, 'LineStyle','-', 'LineWidth',2, 'Color', 'black')
plot(edges1, N1, 'LineStyle','--', 'LineWidth',2, 'Color', 'black')
plot(edges1min, N1min, 'LineStyle',':', 'LineWidth',2, 'Color', 'black')
if j==1
    ylim([0 10])
end
set(gca,'FontSize',18)
xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('PDF', 'FontSize',25, 'interpreter', 'latex')
if (j<4)
title(append('Null and $p$-hacked distributions: $h = $ ', int2str(H(j))),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end
if (j==5)
    title(append('Null and $p$-hacked distributions: $h \sim \chi^2(1)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end
if (j==4)
    title(append('Null and $p$-hacked distributions: $h \sim \widehat\Pi(h)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end
lgd = legend('No $p$-hacking','Threshold', 'Minimum', 'Orientation','horizontal', 'interpreter', 'latex', 'NumColumns',3);
axes('position',[.35 .34 .45 .45])
box on 
indexOfInterest = (edges0 < 0.2) & (edges0 > 0.0001); 
plot(edges0(indexOfInterest),N0(indexOfInterest), '-' , 'LineWidth', 2, 'color', 'black') % plot on new axes
hold on
plot(edges0(indexOfInterest),N1(indexOfInterest), '--' , 'LineWidth', 2, 'color', 'black')
plot(edges0(indexOfInterest),N1min(indexOfInterest), ':' , 'LineWidth', 2, 'color', 'black')
%ylim([0 40])
xlim([0 0.1])
saveas(gcf,append('DGPs/IVSelection/','PivF', int2str(H(j)), int2str(K(k)), 'all',int2str(s),'sided', int2str(GeneralToSpecific)), 'epsc')

%%
end
end
end
close all
end
end
save('DGPs/Bias_struct.mat', 'BIAS')

