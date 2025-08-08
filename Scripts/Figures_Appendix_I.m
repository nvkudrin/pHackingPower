clear all
%% Load Relative_Bias.csv
bias_data = readtable('csvFiles/ForAppendixI/Relative_Bias.csv', 'Delimiter', ',', 'ReadVariableNames', false);

% Extract test names and values
test_names_bias = bias_data.Var1;
bias_values = bias_data.Var2;

% Convert to proper data types if needed
if iscell(test_names_bias)
    test_names_bias = string(test_names_bias);
end
if iscell(bias_values)
    bias_values = str2double(bias_values);
end

% Create a table with proper column names
bias_table = table(test_names_bias, bias_values, ...
    'VariableNames', {'DGP', 'AverageAbsoluteBias'});

%% Load RejectionRates_many_h.csv

% Read the rejection rates data - read all as text first to handle the header properly
rejection_data = readtable('csvFiles/Power_Calculations/RejectionRates_many_h.csv', 'Delimiter', ',', 'ReadVariableNames', true);
rejection_data.Properties.VariableNames{1} = 'DGP';


% Remove the columns we don't want: fail_cs1, fail_csub, fail_cs2b, and FM
columns_to_remove = {'fail_cs1', 'fail_csub', 'fail_cs2b', 'FM'};
for i = 1:length(columns_to_remove)
    if any(strcmp(rejection_data.Properties.VariableNames, columns_to_remove{i}))
        rejection_data.(columns_to_remove{i}) = [];
        fprintf('Removed column: %s\n', columns_to_remove{i});
    end
end

% Merge
merged_data = innerjoin(bias_table, rejection_data, 'Keys', 'DGP');



%%
%% Plot Bias vs Rejection Rates for Threshold and Minimum Data
mkdir('Figures/Extensions/Appendix_I')
% Define line specifications
lineSpecs = {
    {'-', 3, 'black', 'none'}      % LocBin
    {'-.', 3, 'blue', 'none'}      % Disc
    {'--', 3, [0 0.5 0], 'none'}   % CS_1
    {'-.', 3, [0 0.5 0], 'o'}      % CS_UB
    {'-', 3, [0 0.5 0], '*'}       % CS_2B
    {'-.', 3, 'red', 'x'}          % LCM
};

% Test names for legend
test_names = {'Binomial', 'Discontinuity', 'CS1', 'CSUB', 'CS2B', 'LCM'};

if istable(merged_data)
    test_names_col = merged_data{:, 1};
    bias_values = merged_data{:, 2};
    test_results = merged_data{:, 3:8};
else
    test_names_col = merged_data(:, 1);
    bias_values = merged_data(:, 2);
    test_results = merged_data(:, 3:8);
end

% Separate Threshold and Minimum data
threshold_idx = contains(string(test_names_col), 'Threshold');
minimum_idx = contains(string(test_names_col), 'Minimum');

threshold_bias = bias_values(threshold_idx);
threshold_results = test_results(threshold_idx, :);

minimum_bias = bias_values(minimum_idx);
minimum_results = test_results(minimum_idx, :);

% Sort data by bias values for smooth curves
[threshold_bias_sorted, thresh_sort_idx] = sort(threshold_bias);
threshold_results_sorted = threshold_results(thresh_sort_idx, :);

[minimum_bias_sorted, min_sort_idx] = sort(minimum_bias);
minimum_results_sorted = minimum_results(min_sort_idx, :);

%% Figure 1: Threshold Data
fig = figure('Visible','off'); hold on;

for i = 1:6  % 6 tests
    lineStyle = lineSpecs{i}{1};
    lineWidth = lineSpecs{i}{2};
    lineColor = lineSpecs{i}{3};
    markerStyle = lineSpecs{i}{4};
   
    plot(threshold_bias_sorted, threshold_results_sorted(:, i), ...
            'LineStyle', lineStyle, 'LineWidth', lineWidth, 'Color', lineColor, ...
            'Marker', markerStyle, 'MarkerSize', 12, 'MarkerFaceColor', 'none');

end

% Legend and formatting
lgd = legend(test_names, 'Location', 'northeast');
lgd.FontSize = 16;

% Axes & Labels
set(gca, 'FontSize', 18);
xlabel('Relative Bias', 'FontSize', 25, 'interpreter', 'latex');
ylabel('Rejection Rate', 'FontSize', 25, 'interpreter', 'latex');
ylim([0 1]);
xticks([0, 0.1, 0.2, 0.3, 0.4])
xticklabels({'0%','10%','20%', '30%', '40%'})
grid on;
title('Covariate Selection: $K = 3, \tau = 0.25$, Thresholding', 'FontSize', 20, 'interpreter', 'latex');
saveas(fig, 'Figures/Extensions/Appendix_I/PowerBias_Thresholding', 'epsc')
hold off;

%% Figure 2: Minimum Data
fig = figure('Visible','off'); hold on;

for i = 1:6  % 6 tests
    lineStyle = lineSpecs{i}{1};
    lineWidth = lineSpecs{i}{2};
    lineColor = lineSpecs{i}{3};
    markerStyle = lineSpecs{i}{4};
    
    plot(minimum_bias_sorted, minimum_results_sorted(:, i), ...
            'LineStyle', lineStyle, 'LineWidth', lineWidth, 'Color', lineColor, ...
            'Marker', markerStyle, 'MarkerSize', 12, 'MarkerFaceColor', 'none');

end

% Legend and formatting
lgd = legend(test_names, 'Location', 'northeast');
lgd.FontSize = 16;

% Axes & Labels
set(gca, 'FontSize', 18);
xlabel('Relative Bias', 'FontSize', 25, 'interpreter', 'latex');
ylabel('Rejection Rate', 'FontSize', 25, 'interpreter', 'latex');
ylim([0 1]);
xticks([0, 0.1, 0.2, 0.3, 0.4])
xticklabels({'0%','10%','20%', '30%', '40%'})
grid on;
title('Covariate Selection: $K = 3, \tau = 0.9$, Minimum', 'FontSize', 20,'interpreter', 'latex');
saveas(fig, 'Figures/Extensions/Appendix_I/PowerBias_Minimum', 'epsc')
hold off;