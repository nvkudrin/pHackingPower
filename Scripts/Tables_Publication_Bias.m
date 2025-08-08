% Clear workspace and command window
clear all; clc
mkdir('LaTeX_Tables');
% Setup
% Define filtering criteria for selecting simulation results
% - targetDGP: the design type (e.g., 'Covariate')
% - targetSided: '1sided' or '2sided' tests
% - targetGTS: '1' = General-to-Specific; '0' = Specific-to-General
% - variants: reporting method ('Threshold' or 'Minimum')
% - methods: names of tests to report in the table header
% - bias_labels: rows are grouped into No Pub Bias, Sharp Pub Bias, Smooth Pub Bias
% - tau_indices: indices (1-based) for tau = 0, 0.5, and 1 within each block
% - tau_label: LaTeX label for the tau header

% Note: Each block in the data corresponds to a different pub bias scenario (21 rows each)
targetDGP = "Covariate";
targetSided = "2sided";
targetGTS = "1";
variants = ["Threshold", "Minimum"];
methods = ["LocBin", "Disc", "CS_1", "CS_UB", "CS_2B", "LCM"];
bias_labels = ["No Pub Bias", "Sharp Pub Bias", "Smooth Pub Bias"];
tau_indices = [1, 11, 21:27];
tau_label = 'Frac. of $p$-hackers';

% Read data from CSV file
% The file contains two header rows:
% - First row: full parameter ID (e.g., Covariate_03_2sided_1_Threshold)
% - Second row: short test name (e.g., CS_1)
% Then 27 rows of numerical values (21 tau levels for no pub bias and 3 tau levels × 2,, bias types)

%opts = detectImportOptions('csvFiles/Power_Calculations/RejectionRates_main_July21.csv');
%opts = setvartype(opts, 'char');
%T = readtable('csvFiles/Power_Calculations/RejectionRates_main_July21.csv', opts);
%T = table2cell(T);

RejectionRates = readcell('csvFiles/Power_Calculations/RejectionRates_main.csv');
RejectionRates(:, 1:21:size(RejectionRates, 2))=[];

headers     = string(RejectionRates(1, :));   % Full parameter names
testNames  = string(RejectionRates(2, :));   % Method names (e.g., CS_1)
data        = cell2mat(RejectionRates(3:end, :));  % Data matrix: 27 × num_methods
data        = data(tau_indices, :);

% for csv tables
csv_tab2 = [];
csv_tab4 = [];
% ========= SMALL TABLE FOR h=0 ,,
latex = {};
latex{end+1} = '\begin{table}[H]';
latex{end+1} = '\begin{center}';
latex{end+1} = '\caption{The effect of publication bias}';
latex{end+1} = '\label{tab:publication_bias_1}';
latex{end+1} = '\onehalfspacing';
latex{end+1} = '\footnotesize';
latex{end+1} = '\begin{tabular}{lccccccccccccccccccccccc}';
latex{end+1} = '\toprule';
latex{end+1} = ' & \multicolumn{18}{c}{Test} \\ \cline{2-19}';
latex{end+1} = [' & \multicolumn{3}{c}{Binomial}' ...
                ' & \multicolumn{3}{c}{Discontinuity}' ...
                ' & \multicolumn{3}{c}{CS1}' ...
                ' & \multicolumn{3}{c}{CSUB}' ...
                ' & \multicolumn{3}{c}{CS2B}' ...
                ' & \multicolumn{3}{c}{LCM} \\'];
latex{end+1} = [tau_label ...
                ' & 0 & 0.5 & 1 & 0 & 0.5 & 1' ...
                ' & 0 & 0.5 & 1 & 0 & 0.5 & 1' ...
                ' & 0 & 0.5 & 1 & 0 & 0.5 & 1 \\ \hline'];

% Loop over reporting variants (Thresholding and Minimum)
for v = 1:2
    variant = variants(v);
    variantLabel = variant;
    if variant == "Threshold"
        variantLabel = "Thresholding";
    end
    latex{end+1} = sprintf('\\textit{} & \\multicolumn{18}{c}{%s} \\\\ \\cline{2-19}', variantLabel);

    for b = 1:3
        row_offset = (b - 1) * 3;
        line = bias_labels(b);
        for m = 1:length(methods)
            hcode = "03"; % Encodes h = 0, K = 3 (2 digits, not 'h0K3')
            method = methods(m);
            pattern = targetDGP + "_" + hcode + "_" + targetSided + "_" + targetGTS + "_" + variant;
            colIdx = find(contains(headers, pattern) & testNames == method);            
            for t = 1:3
                val = data(row_offset + t, colIdx);
                line = line + " & " + sprintf('%.2f', val);
                csv_tab2 = [csv_tab2, val];
            end
        end
        latex{end+1} = line + " \\";
    end
end

latex{end+1} = '\bottomrule';
latex{end+1} = '\end{tabular}';
latex{end+1} = '\end{center}';
latex{end+1} = ['\doublespacing', ...
    '\textit{Notes:} Table shows the impact of publication bias on the power of the tests when $p$-hacking is based on covariate selection with $K=3$ (general-to-specific, two-sided tests) and $h=0$. The results are based on 5,000 simulation repetitions.'];
latex{end+1} = '\end{table}';

% Write LaTeX code to file
fid = fopen('LaTeX_Tables/Table_2.tex', 'w');
for i = 1:length(latex)
    fprintf(fid, '%s\n', latex{i});
end
fclose(fid);
disp('Saved: Table_2.tex');

% Also create a csv table
csv_tab2 = reshape(csv_tab2, 18, 6)';
csv_tab2 = num2cell(csv_tab2);
col1 = {'Test', 'fraction of p-hackers', 'Thresholding', 'No pub bias', 'Sharp pub bias', 'Smooth pub bias', 'Minimum', 'No pub bias', 'Sharp pub bias', 'Smooth pub bias'}';
gap = repmat({' '}, 1, 18);
tau_row = repmat({'0', '0.5', '1'}, 1, 6);
test_row = reshape(repmat({'Binomial', 'Discontinuity', 'CS1', 'CSUB', 'CS2B', 'LCM'}, 3, 1), 1, 18);
csv_tab2 = [col1, [test_row; tau_row; gap; csv_tab2(1:3,:); gap; csv_tab2(4:6,:)]];
writecell(csv_tab2, 'csvFiles/Power_Calculations/Table2.csv');

% ========= PANEL TABLE FOR h = 1,2,'3' ==========
latex = {};
HKcodes = ["13", "23", "33"];  % Encoded as 2-digit string: h = 1,2,'3' with K = 3
Hvals = [1, 2, 3];

latex{end+1} = '\begin{table}[H]';
latex{end+1} = '\begin{center}';
latex{end+1} = '\caption{The effect of publication bias for $h=1$, $h=2$, and $h\sim \widehat{\Pi}$}';
latex{end+1} = '\label{tab:publication_bias_appendix}';
latex{end+1} = '\footnotesize';
latex{end+1} = '\begin{tabular}{lccccccccccccccccccccccc}';
latex{end+1} = '\toprule';
latex{end+1} = ' & \multicolumn{18}{c}{Test} \\ \cline{2-19}';
latex{end+1} = [' & \multicolumn{3}{c}{Binomial}' ...
                ' & \multicolumn{3}{c}{Discontinuity}' ...
                ' & \multicolumn{3}{c}{CS1}' ...
                ' & \multicolumn{3}{c}{CSUB}' ...
                ' & \multicolumn{3}{c}{CS2B}' ...
                ' & \multicolumn{3}{c}{LCM} \\'];
latex{end+1} = [tau_label ...
                ' & 0 & 0.5 & 1 & 0 & 0.5 & 1' ...
                ' & 0 & 0.5 & 1 & 0 & 0.5 & 1' ...
                ' & 0 & 0.5 & 1 & 0 & 0.5 & 1 \\ \hline'];

% Loop over h = 1, 2, 3 (or h ~ distribution)
for hi = 1:3
    h = Hvals(hi);
    hcode = HKcodes(hi);  % string encoding of h and K = 3
    if h == 3
        hlabel = '$h\sim\widehat{\Pi}$';
    else
        hlabel = "$h=" + int2str(h) + "$";
    end
    latex{end+1} = sprintf('\\multicolumn{1}{c}{} & \\multicolumn{18}{c}{%s} \\\\', hlabel);

    for v = 1:2
        variant = variants(v);
        variantLabel = variant;
        if variant == "Threshold"
            variantLabel = "Thresholding";
        end
        latex{end+1} = sprintf('\\textit{} & \\multicolumn{18}{c}{%s} \\\\ \\cline{2-19}', variantLabel);

        for b = 1:3
            row_offset = (b - 1) * 3;
            line = bias_labels(b);
            for m = 1:length(methods)
                method = methods(m);
                pattern = targetDGP + "_" + hcode + "_" + targetSided + "_" + targetGTS + "_" + variant;
                colIdx = find(contains(headers, pattern) & testNames == method);
                for t = 1:3
                    val = data(row_offset + t, colIdx);
                    line = line + " & " + sprintf('%.2f', val);
                    csv_tab4 = [csv_tab4, val];
                end
            end
            latex{end+1} = line + " \\";
        end
    end
    if hi<3
        latex{end+1} = '\hline';
    end
end

latex{end+1} = '\bottomrule';
latex{end+1} = '\end{tabular}';
latex{end+1} = '\end{center}';
latex{end+1} = ['\footnotesize{\textit{Notes:} Table shows the impact of publication bias on the power of the tests when $p$-hacking is based on covariate selection (general-to-specific, two-sided tests) with $K=3$ and $h=1$, $h=2$, and $h\sim \widehat{\Pi}$. The results are based on 5,000 simulation repetitions.}'];
latex{end+1} = '\end{table}';

fid = fopen('LaTeX_Tables/Table_4.tex', 'w');
for i = 1:length(latex)
    fprintf(fid, '%s\n', latex{i});
end
fclose(fid);
disp('Saved: Table_4.tex');

% Also create a csv table
csv_tab4 = reshape(csv_tab4, 18, 18)';
csv_tab4 = num2cell(csv_tab4);
col1 = {'Test', 'fraction of p-hackers', 'Thresholding, h = 1', 'No pub bias', 'Sharp pub bias', 'Smooth pub bias', 'Minimum, h = 1', 'No pub bias', 'Sharp pub bias', 'Smooth pub bias',...
    'Thresholding, h = 2', 'No pub bias', 'Sharp pub bias', 'Smooth pub bias', 'Minimum, h = 2', 'No pub bias', 'Sharp pub bias', 'Smooth pub bias',...
    'Thresholding, h ~ Pi_hat', 'No pub bias', 'Sharp pub bias', 'Smooth pub bias', 'Minimum, h ~ Pi_hat', 'No pub bias', 'Sharp pub bias', 'Smooth pub bias'}';
gap = repmat({' '}, 1, 18);
tau_row = repmat({'0', '0.5', '1'}, 1, 6);
test_row = reshape(repmat({'Binomial', 'Discontinuity', 'CS1', 'CSUB', 'CS2B', 'LCM'}, 3, 1), 1, 18);
csv_tab4 = [col1, [test_row; tau_row; gap; csv_tab4(1:3,:); gap; csv_tab4(4:6,:); gap; csv_tab4(7:9,:); gap; csv_tab4(10:12,:);gap; csv_tab4(13:15,:); gap; csv_tab4(16:18,:)]];
writecell(csv_tab4, 'csvFiles/Power_Calculations/Table4.csv');

