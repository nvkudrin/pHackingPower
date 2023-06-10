clear all

%S = importdata('DGPs/RejectionRates_April6.csv');
%V = S.data;


S = readtable("PowerCurves/RejectionRates_Apr5.csv");
Vfull=table2array(S);
Vfull = Vfull(:,2:end);

tau = [0:0.05:1];
tau_s = [0.25, 0.75];
num_per_dgp = 11;

load('DGPs/Bias_struct.mat')

H = [0, 1, 2, 3, 4];
K = [3,5,7];

for sided = [1,2]
    V = Vfull(1:21, (1+(sided-1)*1056/2): (1056/(3-sided)));
    if sided == 1
        side = "1sided";
    end
    if sided == 2
        side = "2sided";
    end
for sel = 0:1
for mnm=0:1
tau_j = find(tau==tau_s(mnm+1));
Bin = [];
Fisher = [];
rdd = [];
CS1 = [];
CS1B = [];
CS2B = [];
LCM = [];
for fig = (24*mnm+1):(24*mnm+20)
    Bin = [Bin, V(tau_j, 2+(fig-1)*num_per_dgp)];
    Fisher = [Fisher, V(tau_j, 3+(fig-1)*num_per_dgp)];
    rdd = [rdd, V(tau_j, 4+(fig-1)*num_per_dgp)];
    CS1 = [CS1, V(tau_j, 5+(fig-1)*num_per_dgp)];
    CS1B = [CS1B, V(tau_j, 6+(fig-1)*num_per_dgp)];
    CS2B = [CS2B, V(tau_j, 7+(fig-1)*num_per_dgp)];
    LCM = [LCM, V(tau_j, 8+(fig-1)*num_per_dgp)];
end

    if (sel == 0)
        Ind = 13:20;
    end
     if (sel == 1)
       Ind = 1:12;
     end
if (mnm==0)
     B = [BIAS.sel(1, 1:12,sided), BIAS.iv(1,1:8,sided)];
end
if (mnm==1)
B = [BIAS.sel(2, 1:12,sided), BIAS.iv(2,1:8,sided)];
end

figure(1)
hold on
h1 = scatter(B(Ind), Bin(Ind), 60, 'MarkerEdgeColor',[0 0 0],...
              'MarkerFaceColor',[0 0 0],...
              'LineWidth',2)

lgd.FontSize = 16;

h2=scatter(B(Ind), rdd(Ind), 60, '+', 'MarkerEdgeColor',[0 0 1],...
              'MarkerFaceColor',[0 0 1],...
              'LineWidth',2)

lgd.FontSize = 16;

h3 = scatter(B(Ind ), CS2B(Ind), 60,'s', 'MarkerEdgeColor',[0 0.5 0],...
              'MarkerFaceColor',[0 0.5 0],...
              'LineWidth',2)
line = lsline
line(1).Color = [0 0.5 0];
line(2).Color = [0 0 1];
line(3).Color = [0 0 0];

line(1).LineWidth = 2;
line(2).LineWidth = 2;
line(3).LineWidth = 2;

legend([h1 h2 h3],{'Binomial','Discontinuity', 'CS2B'})
lgd.FontSize = 16;
legend('Location','northwest')

set(gca,'FontSize',18)
xlabel('Average Bias', 'FontSize',25, 'interpreter', 'latex')
ylabel('Power', 'FontSize',25, 'interpreter', 'latex')

if (mnm==0)
    ylim([0,1])
    if (sel==0)
        if sided == 1
        xlim([min(B),0.05])
        end
        if sided == 2
                xlim([min(B),0.015])
            end
    title('IV Selection, Thresholding ($\tau = 0.25$)','fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
    saveas(gcf,append('Scatters/','IV_scatter_t', side), 'epsc')
    end
    if (sel==1)
        if sided == 1
        xlim([0.03,0.1])
        end
        if sided == 2
                xlim([0.03,0.07])
            end
    title('Covariate Selection, Thresholding ($\tau = 0.25$)','fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
    saveas(gcf,append('Scatters/','CovSel_scatter_t',side), 'epsc')
    end

end

if (mnm==1)
    ylim([0,0.8])
    if (sel==0)
        if sided == 1
        xlim([min(B),0.05])
        end
          if sided == 2
                xlim([min(B),0.027])
            end
    title('IV Selection, Minimum ($\tau = 0.75$)','fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
    saveas(gcf,append('Scatters/','IV_scatter_m',side), 'epsc')
    end
        if (sel==1)
            if sided == 1
            xlim([0.03,0.1])
            end
            if sided == 2
                xlim([0.03,0.07])
            end
    title('Covariate Selection, Minimum ($\tau = 0.75$)','fontweight','bold', 'FontSize',20, 'interpreter', 'latex')
        saveas(gcf,append('Scatters/','CovSel_scatter_m', side), 'epsc')
        end

end
close all
end
end
end