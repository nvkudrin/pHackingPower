clear all
mkdir('Figures/Extensions/Appendix_F')
mkdir('Figures/Extensions/Appendix_G')
mkdir('Figures/Extensions/Appendix_H')
for k = 1:3
    if k ==1
filename = 'csvFiles/Power_Calculations/RejectionRates_different_rho.csv';
    end
    if k ==2
filename = 'csvFiles/Power_Calculations/RejectionRates_different_pmis.csv';
    end
    if k ==3
filename = 'csvFiles/Power_Calculations/RejectionRates_different_N.csv';
    end
% Read all data 
raw_data = readcell(filename);

headers = raw_data(1, :);
tests = raw_data(2, :);

keep_columns = contains(headers, {'Threshold', 'Minimum'});
filtered_headers = unique(headers(keep_columns));
% Extract the data (starting from row 3)
data_rows = raw_data(3:end, :);
x_values = cell2mat(data_rows(:, 1));

for j = 1:numel(filtered_headers)
    dgp = cell2mat(filtered_headers(j));
    dgp_parts = split(dgp);

    h = str2double(dgp_parts{1});
    
    dgp_columns = contains(headers, dgp);
    dgp_tests = contains(tests, {'LocBin', 'Disc', 'CS_1', 'CS_UB', 'CS_2B', 'LCM'});
    idx = (dgp_columns & dgp_tests);
    data_fig = cell2mat(data_rows(:, idx));

    % Create the plot with exact specifications
    fig = figure('Visible','off'); hold on;
    
    % Plot 1: Binomial (column 1) - solid line, black
    plot(x_values, data_fig(:, 1), 'LineWidth', 3, 'color', 'black');
    hold on;
    
    % Plot 2: Discontinuity (column 2) - dash-dot line, blue
    plot(x_values, data_fig(:, 2), '-.', 'LineWidth', 3, 'color', 'blue', 'MarkerSize', 12);
    
    % Plot 3: CS1 (column 3) - dashed line, dark green
    plot(x_values, data_fig(:, 3), '--', 'LineWidth', 3, 'color', [0 0.5 0], 'MarkerSize', 12);
    
    % Plot 4: CSUB (column 4) - dash-dot with circles, dark green
    plot(x_values, data_fig(:, 4), '-.o', 'LineWidth', 3, 'color', [0 0.5 0], 'MarkerSize', 12);
    
    % Plot 5: CS2B (column 5) - dash-dot with stars, dark green
    plot(x_values, data_fig(:, 5), '-.*', 'LineWidth', 3, 'color', [0 0.5 0], 'MarkerSize', 12);
    
    % Plot 6: LCM (column 6) - dash-dot with x, red
    plot(x_values, data_fig(:, 6), '-.x', 'LineWidth', 3, 'color', 'red', 'MarkerSize', 12);
    
    % Add legend only if fig == 1
    if strcmp(dgp, '0 Threshold')
        lgd = legend('Binomial', 'Discontinuity', 'CS1', 'CSUB', 'CS2B', 'LCM', 'Location', 'northwest');
        lgd.FontSize = 16;
    end

    if contains(dgp, 'Threshold')
        title_item = 'Thresholding';
    elseif contains(dgp, 'Minimum')
        title_item = 'Minimum';
    end
    
    if k==1
    % Set axes and labels with exact specifications
    set(gca, 'FontSize', 18);
    xlabel('$$\rho$$', 'FontSize', 25, 'interpreter', 'latex');
    ylabel('Rejection Rate', 'FontSize', 25, 'interpreter', 'latex');
    ylim([0 1]);
    xlim([0 1]);
    xticks([0, 0.25, 0.5, 0.75, 1]);
    xticklabels({'0', '0.25', '0.50', '0.75', '1'});

    if h==3
        title(append('$h\sim \widehat{\Pi}$, ', title_item, ' $(\tau = 0.5)$'), 'fontweight', 'bold', 'FontSize', 20, 'interpreter', 'latex');
    else
        title(append('$h = ',int2str(h),'$, ', title_item, ' $(\tau = 0.5)$'), 'fontweight', 'bold', 'FontSize', 20, 'interpreter', 'latex');
    end
    saveas(fig, append('Figures/Extensions/Appendix_F/Appendix_F_',dgp) , 'epsc')
       hold off;
    end

        if k==2
    % Set axes and labels with exact specifications
    set(gca, 'FontSize', 18);
    xlabel('$$p_{mis}$$', 'FontSize', 25, 'interpreter', 'latex');
    ylabel('Rejection Rate', 'FontSize', 25, 'interpreter', 'latex');
    ylim([0 1]);
    xlim([0 1]);
    xticks([0, 0.25, 0.5, 0.75, 1]);
    xticklabels({'0', '0.25', '0.50', '0.75', '1'});

    if h==3
        title(append('$h\sim \widehat{\Pi}$, ', title_item, ' $(\tau = 0.5)$'), 'fontweight', 'bold', 'FontSize', 20, 'interpreter', 'latex');
    else
        title(append('$h = ',int2str(h),'$, ', title_item, ' $(\tau = 0.5)$'), 'fontweight', 'bold', 'FontSize', 20, 'interpreter', 'latex');
    end
    saveas(fig, append('Figures/Extensions/Appendix_G/Appendix_G_',dgp) , 'epsc')
       hold off;
    end

        if k==3
    % Set axes and labels with exact specifications
    set(gca, 'FontSize', 18);
    xlabel('Sample Size', 'FontSize', 25, 'interpreter', 'latex');
    ylabel('Rejection Rate', 'FontSize', 25, 'interpreter', 'latex');
    xlim([400 5000]);
    ylim([0 1]);
    xticks([400, 1000, 2000, 3000, 4000, 5000]);
    xticklabels({'400', '1000', '2000', '3000', '4000', '5000'});

    if h==3
        title(append('$h\sim \widehat{\Pi}$, ', title_item, ' $(\tau = 0.5)$'), 'fontweight', 'bold', 'FontSize', 20, 'interpreter', 'latex');
    else
        title(append('$h = ',int2str(h),'$, ', title_item, ' $(\tau = 0.5)$'), 'fontweight', 'bold', 'FontSize', 20, 'interpreter', 'latex');
    end
    saveas(fig, append('Figures/Extensions/Appendix_H/Appendix_H_',dgp) , 'epsc')
       hold off;
    end

end
end
