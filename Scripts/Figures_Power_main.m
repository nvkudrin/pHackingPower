% This program creates Figures showing the power of tests for different DGPS
clear all

% Read all data 
RejectionRates = readcell('csvFiles/Power_Calculations/RejectionRates_main.csv');
RejectionRates(:, 1:21:size(RejectionRates, 2))=[];

folders = {'csvFiles/Distributions/CovariateSelection', 'csvFiles/Distributions/IVSelection', 'csvFiles/Distributions/LagLengthSelection', 'csvFiles/Distributions/ClusterSelection'};
processedNames = {};


for i = 1:length(folders)
    files = dir(fullfile(folders{i}, 'P0_*.csv'));  % Only files starting with 'P0_' and ending in '.csv'
    for j = 1:length(files)
        name = files(j).name;

        % Remove 'P0_' prefix
        if startsWith(name, 'P0_')
            name = extractAfter(name, 'P0_');
        end

        % Remove '.csv' extension
        name = erase(name, '.csv');

        processedNames{end+1} = name; %#ok<SAGROW>
    end
end

FilteredByName = cell(1, numel(processedNames));

for k = 1:numel(processedNames)
    name = processedNames{k};
    colMatches = cellfun(@(x) contains(string(x), name), RejectionRates(1, :));
    FilteredByName{k} = RejectionRates(:, colMatches);
end

targetTests = ["LocBin", "Disc", "CS_1", "CS_UB", "CS_2B", "LCM"];
variants = ["Threshold", "Minimum"];
tau = 0:0.05:1;
H = [0, 1, 2, 3];
K_vals = [3, 5, 7];

lineSpecs = {
    {'-',  3, 'black',       'none'}     % LocBin
    {'-.', 3, 'blue',        'none'}     % Disc
    {'--', 3, [0 0.5 0],     'none'}     % CS_1
    {'-.', 3, [0 0.5 0],     'o'}        % CS_UB
    {'-',  3, [0 0.5 0],     '*'}        % CS_2B
    {'-.', 3, 'red',         'x'}        % LCM
};

% Define output folders
folderMap = containers.Map(...
    ["Covariate", "IV", "IVF", "LagLength", "Cluster"], ...
    ["Figures/PowerCurves/CovariateSelection", ...
     "Figures/PowerCurves/IVSelection", ...
     "Figures/PowerCurves/IVSelection", ...
     "Figures/PowerCurves/LagLengthSelection", ...
     "Figures/PowerCurves/ClusterSelection"]);

for k = 1:numel(FilteredByName)
    data = FilteredByName{k};
    headers = string(data(1, :));
    tests = string(data(2, :));
    numericData = cell2mat(data(3:23, :));
    baseName = processedNames{k};

    % Parse baseName
    underscoreSplit = split(baseName, "_");
    designType = underscoreSplit(1);
    code = underscoreSplit(2);

    if designType == "Covariate" || designType == "IV" || designType == "IVF"
        h = str2double(extractBetween(code, 1, 1));
        Kval = str2double(extractBetween(code, 2, 2));
        j = find(H == h);
        k_index = find(K_vals == Kval);
    elseif designType == "LagLength" || designType == "Cluster"
        h = str2double(code);
        j = find(H == h);
        k_index = NaN;
    else
        error("Unknown design type: " + designType);
    end

    % Create figure folder
    figFolder = folderMap(char(designType));
    if ~exist(figFolder, 'dir')
        mkdir(figFolder);
    end

    for v = 1:2
        variant = variants(v);
        minimum = double(variant == "Minimum");
        fullPrefix = baseName + "_" + variant;
        isMatch = startsWith(headers, fullPrefix);
        matchingTests = tests(isMatch);
        idxMatch = find(isMatch);

        % Select desired methods
        idxSelected = [];
        for s = targetTests
            match = find(matchingTests == s, 1);
            if ~isempty(match)
                idxSelected(end+1) = idxMatch(match); %#ok<SAGROW>
            end
        end

        if length(idxSelected) == 6
            fig = figure('Visible','off'); hold on;
            for i = 1:6
                col = idxSelected(i);
                spec = lineSpecs{i};
                plot(tau, numericData(:, col), ...
                    spec{1}, 'LineWidth', spec{2}, ...
                    'Color', spec{3}, ...
                    'Marker', spec{4}, 'MarkerSize', 12);
            end

            % --- Legend Logic ---
            if j == 1
                if (designType == "Covariate" || designType == "IV" || designType == "IVF")
                    if (k_index == 1 && minimum == 0) || (k_index > 1 && minimum == 1)
                        lgd = legend('Binomial', 'Discontinuity', ...
                            'CS1','CSUB', 'CS2B','LCM', 'Location','northwest');
                        lgd.FontSize = 16;
                    end
                elseif designType == "LagLength" || designType == "Cluster"
                    if minimum == 0
                        lgd = legend('Binomial', 'Discontinuity', ...
                            'CS1','CSUB', 'CS2B','LCM', 'Location','northwest');
                        lgd.FontSize = 16;
                    end
                end
            end

            % Axes & Labels
            set(gca, 'FontSize', 18)
            xlabel('$$\tau$$ (fraction of $p$-hackers)', ...
                   'FontSize', 25, 'interpreter', 'latex')
            ylabel('Rejection Rate', ...
                   'FontSize', 25, 'interpreter', 'latex')
            ylim([0 1])
            xlim([0 1])

            % Title
            if designType == "Covariate" || designType == "IV" || designType == "IVF"
                if j == 4
                    titleText = "$K = $ " + int2str(Kval) + ", $h \sim \widehat \Pi$";
                else
                    titleText = "$K = $ " + int2str(Kval) + ", $h = $ " + int2str(h);
                end
            else
                if j == 4
                    titleText = "$h \sim \widehat \Pi$";
                else
                    titleText = "$h = $ " + int2str(h);
                end
            end
            title(titleText, 'FontSize', 20, 'FontWeight', 'bold', 'Interpreter', 'latex');

            % Save figure
            saveas(fig, fullfile(figFolder, fullPrefix + ".eps"), 'epsc')
            close(fig);
        end
    end
end