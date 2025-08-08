%% This file runs all m-files that generate remaining figures and tables from the paper
clear all
addpath('Functions')
addpath('Scripts')
%% Figures
mkdir('Figures')
% Analytical figures
for k = [1, 2, 3, 7:17]
    filename = sprintf('Figure_%d.m', k);
    run(filename);
end
% Main power curves
run Figures_Power_main
% Power curves for Appendices F, G, H
run Figures_Appendix_FGH
% Bias-Power plots in Appendix I
run Figures_Appendix_I

%% Tables
mkdir('LaTeX_Tables')
% Empirical Application (Table 3)
run Table_Empirical_Application
% Effects of publication bias (Tables 2 and 4)
run Tables_Publication_Bias