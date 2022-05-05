clear all
load('DGPs/MC_distributions_Apr6.mat')
rng(12345)
 Xi = chi2rnd(1, 1, 1000000);
H = [0, 1, 2, 3];
K = [3,5,7];
Up = [0,1];
nobs = 200;

BIAS = struct;
BIAS.var = zeros(2,5);
BIAS.sel = zeros(2, 12);
BIAS.iv = zeros(2, 8);
%%
fig = 0;
b_ind_sel = 0;
b_ind_iv = 0;
for j = 1:4
    close all
    if (j<4)
b = (H(j)/sqrt(nobs))*ones(size(mcout_bias));
b_iv3 = (H(j)/sqrt(nobs)/3)*ones(size(mcout_bias_iv3));
b_iv5 = (H(j)/sqrt(nobs)/3)*ones(size(mcout_bias_iv5));
b_var = (H(j)/sqrt(nobs))*ones(size(mcout_bias_var));
    end
    if (j==4)
b = repmat(Xi, 127, 1)/sqrt(nobs);
b_iv3 = repmat(Xi, 7, 1)/sqrt(nobs)/3;
b_iv5 = repmat(Xi, 31, 1)/sqrt(nobs)/3;
b_var = repmat(Xi, 5, 1)/sqrt(nobs);
    end
    
 %LagLengthSelection   
    [P0, P1, P1min] = NullAndAlt_var_bic(b_var, mcout_bias_var, mcout_se_var, 1, 0.05, mcout_bic_var);
csvwrite(append('DGPs','P0var', int2str(H(j)), '.csv'), P0)
fig = fig+1;
figure(fig)
histogram(P0, 100, 'Normalization', 'probability')
set(gca,'FontSize',18)
xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('Probability', 'FontSize',25, 'interpreter', 'latex')
if (j<4)
title(append('No $p$-hacking: $h = $ ', int2str(H(j))),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end
if (j==4)
    title(append('No $p$-hacking: $h \sim \chi^2(1)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end
saveas(gcf,append('DGPs/LagLengthSelection/','P0var', int2str(H(j))), 'epsc')

csvwrite(append('DGPs','P1var', int2str(H(j)), '.csv'), P1)

fig = fig+1;
figure(fig)
histogram(P1, 100, 'Normalization', 'probability')
set(gca,'FontSize',18)
xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('Probability', 'FontSize',25, 'interpreter', 'latex')
if (j<4)
title(append('p-hacked (threshold): $h = $ ', int2str(H(j))),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end
if (j==4)
    title(append('p-hacked (threshold): $h \sim \chi^2(1)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end
saveas(gcf,append('DGPs/LagLengthSelection/','P1var', int2str(H(j))), 'epsc')

fig = fig+1;
figure(fig)
histogram(P1min, 100, 'Normalization', 'probability')
set(gca,'FontSize',18)
xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('Probability', 'FontSize',25, 'interpreter', 'latex')
if (j<4)
title(append('p-hacked (minimum): $h = $ ', int2str(H(j))),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end
if (j==4)
    title(append('p-hacked (minimum): $h \sim \chi^2(1)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end
saveas(gcf,append('DGPs/LagLengthSelection/','P1var', int2str(H(j)), 'min'), 'epsc')
csvwrite(append('DGPs','P1var', int2str(H(j)), 'min', '.csv'), P1min) 


for k = 1:3
    fig
%CovariateSelection
[P0, P1,P1min, B0, B1, B1min] = NullAndAlt(b, mcout_bias, mcout_se, 1, 0.05, 7, K(k));
b_ind_sel = b_ind_sel+1;
BIAS.sel(1, b_ind_sel) = mean(B1);
BIAS.sel(2, b_ind_sel) = mean(B1min);

csvwrite(append('DGPs','P0sel', int2str(H(j)), int2str(K(k)), '.csv'), P0)
csvwrite(append('DGPs','P1sel', int2str(H(j)), int2str(K(k)), '.csv'), P1)
csvwrite(append('DGPs','P1sel', int2str(H(j)), int2str(K(k)),'min', '.csv'), P1min)
fig = fig+1;
figure(fig)
histogram(P0, 100, 'Normalization', 'probability')
set(gca,'FontSize',18)
xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('Probability', 'FontSize',25, 'interpreter', 'latex')
if (j<4)
title(append('No $p$-hacking: $K = $ ', int2str(K(k)), ', $h = $ ', int2str(H(j))),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end
if (j==4)
    title(append('No $p$-hacking: $K = $', int2str(K(k)), ', $h \sim \chi^2(1)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end
saveas(gcf,append('DGPs/CovariateSelection/','P0sel', int2str(H(j)), int2str(K(k))), 'epsc')

fig = fig+1;
figure(fig)
histogram(P1, 100, 'Normalization', 'probability')
set(gca,'FontSize',18)
xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('Probability', 'FontSize',25, 'interpreter', 'latex')
if (j<4)
title(append('p-hacked (threshold): $K = $ ', int2str(K(k)), ', $h = $ ', int2str(H(j))),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end
if (j==4)
    title(append('p-hacked (threshold): $K = $ ', int2str(K(k)), ', $h \sim \chi^2(1)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end
saveas(gcf,append('DGPs/CovariateSelection/','P1sel', int2str(H(j)), int2str(K(k))), 'epsc')

fig = fig+1;
figure(fig)
histogram(P1min, 100, 'Normalization', 'probability')
set(gca,'FontSize',18)
xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('Probability', 'FontSize',25)
if (j<4)
title(append('p-hacked (minimum): $K = $ ', int2str(K(k)), ', $h = $ ', int2str(H(j))),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end
if (j==4)
    title(append('p-hacked (minimum): $K = $ ', int2str(K(k)), ', $h \sim \chi^2(1)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end
saveas(gcf,append('DGPs/CovariateSelection/','P1sel', int2str(H(j)), int2str(K(k)), 'min'), 'epsc')

%IVSelection
if (k<3)
    if (k==1)
[P0, P1,P1min, B0, B1, B1min] = NullAndAlt(b_iv3, mcout_bias_iv3, mcout_se_iv3, 1, 0.05, 3, 3);
b_ind_iv = b_ind_iv+1;
    end
    if (k==2)
[P0, P1,P1min, B0, B1, B1min] = NullAndAlt(b_iv5, mcout_bias_iv5, mcout_se_iv5, 1, 0.05, 5, 5);
b_ind_iv = b_ind_iv+1;  
    end
BIAS.iv(1, b_ind_iv) = mean(B1);
BIAS.iv(2, b_ind_iv) = mean(B1min);
csvwrite(append('DGPs','P0iv', int2str(H(j)), int2str(K(k)), '.csv'), P0)
fig = fig+1;
figure(fig)
histogram(P0, 100, 'Normalization', 'probability')
set(gca,'FontSize',18)
xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('Probability', 'FontSize',25, 'interpreter', 'latex')
if (j<4)
title(append('No $p$-hacking: $K = $ ', int2str(K(k)), ', $h = $ ', int2str(H(j))),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end
if (j==4)
    title(append('No $p$-hacking: $K = $ ', int2str(K(k)), ', $h \sim \chi^2(1)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end
saveas(gcf,append('DGPs/IVSelection/','P0iv', int2str(H(j)), int2str(K(k))), 'epsc')

csvwrite(append('DGPs','P1iv', int2str(H(j)), int2str(K(k)), '.csv'), P1)

fig = fig+1;
figure(fig)
histogram(P1, 100, 'Normalization', 'probability')
set(gca,'FontSize',18)
xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('Probability', 'FontSize',25, 'interpreter', 'latex')
if (j<4)
title(append('p-hacked (threshold): $K = $ ', int2str(K(k)), ', $h = $ ', int2str(H(j))),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end
if (j==4)
    title(append('p-hacked (threshold): $K = $ ', int2str(K(k)), ', $h \sim \chi^2(1)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end
saveas(gcf,append('DGPs/IVSelection/','P1iv', int2str(H(j)), int2str(K(k))), 'epsc')


csvwrite(append('DGPs','P1iv', int2str(H(j)), int2str(K(k)),'min', '.csv'), P1min)

fig = fig+1;
figure(fig)
histogram(P1min, 100, 'Normalization', 'probability')
set(gca,'FontSize',18)
xlabel('$$p$$', 'FontSize',25, 'interpreter', 'latex')
ylabel('Probability', 'FontSize',25, 'interpreter', 'latex')
if (j<4)
title(append('p-hacked (minimum): $K = $ ', int2str(K(k)), ', $h = $ ', int2str(H(j))),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end
if (j==4)
    title(append('p-hacked (minimum): $K = $ ', int2str(K(k)), ', $h \sim \chi^2(1)$'),'fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
end
saveas(gcf,append('DGPs/IVSelection/','P1iv', int2str(H(j)), int2str(K(k)), 'min'), 'epsc')

end
end
close all
end
save('DGPs/Bias_struct.mat', 'BIAS')