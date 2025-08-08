% Load CSV data
panelA = readtable('csvFiles/Empirical_Application/Table3a.csv', 'ReadRowNames', true, 'ReadVariableNames',true, VariableNamingRule='preserve');
panelB = readtable('csvFiles/Empirical_Application/Table3b.csv', 'ReadRowNames', true, 'ReadVariableNames',true, VariableNamingRule='preserve');
panelA = panelA(1:(end-1),1:13);
panelB = panelB(1:(end-1),1:13);
%panelB.Properties.RowNames = panelB{:,1};  % Set row names from first column
%panelB(:,1) = [];                     % Remove first column from the table

% Observation counts (hardcoded based on table structure)
obs = [5853, 7569, 3148, 5170, 679, 760, 3954, 17786, ...
       11211, 10529, 3341, 17772, 21740];

% Row and column labels
tests = panelA.Properties.RowNames;
subsamples = {'DID', 'RCT', 'RDD', 'IV', 'F<30', 'F>=30', ...
              'Top 5', '!Top 5', '2015', '2018', 'AJQ', ...
              'Star Wars', 'All'};
nTests = numel(tests);
nSubs = numel(subsamples);

% Open file for writing
mkdir('LaTeX_Tables')
fid = fopen('LaTeX_Tables/Table_3.tex', 'w');

% Begin LaTeX table output
fprintf(fid, '\\begin{table}[H]\n');
fprintf(fid, '\\footnotesize\n');
fprintf(fid, '\\caption{Results empirical application}\n');
fprintf(fid, '\\label{tab:application}\n');
fprintf(fid, '\\onehalfspacing\n');
fprintf(fid, '\\begin{tabular}{l%s c}\n', repmat('c', 1, nSubs));
fprintf(fid, '\\toprule\n\n');

% Header rows
fprintf(fid, '\\multicolumn{1}{c}{\\multirow{2}{*}{Test}} & ');
fprintf(fid, '\\multicolumn{%d}{c}{Subsample} & \\multicolumn{1}{c}{\\multirow{2}{*}{Overall RR}} \\\\\n', nSubs);
fprintf(fid, '\\cline{2-%d}\n', nSubs+1);
fprintf(fid, '\\multicolumn{1}{c}{} ');
fprintf(fid, '&%s ', subsamples{:});
fprintf(fid, '\\multicolumn{1}{c}{} ');
fprintf(fid, '\\\\\n\\midrule\n');

% Panel A
fprintf(fid, '& \\multicolumn{%d}{c}{(a) Original (rounded) data: $p$-values} \\\\\n', nSubs+1);
fprintf(fid, '\\cline{2-%d}\n', nSubs+1);
for i = 1:nTests
    row = panelA{i, :};
    fprintf(fid, '%-10s ', tests{i});
    fprintf(fid, '& %.3f ', row);
    fprintf(fid, '& %.3f \\\\\n', mean(row<=0.05)); % Overall RR
end

% Panel B
fprintf(fid, '\\midrule\n');
fprintf(fid, '& \\multicolumn{%d}{c}{(b) Average rejection rates across 1,000 deroundings, 5\\%% significance level} \\\\\n', nSubs+1);
fprintf(fid, '\\cline{2-%d}\n', nSubs+1);
for i = 1:nTests
    row = panelB{i, :};
    fprintf(fid, '%-10s ', tests{i});
    fprintf(fid, '& %.3f ', row);
    fprintf(fid, '& %.3f \\\\\n', mean(row)); % Overall RR
end

% Observation counts row
fprintf(fid, '\\midrule\n');
fprintf(fid, 'Obs ');
fprintf(fid, '& %d ', obs);
fprintf(fid, '& \\\\\n');

% Footer
fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n\n');
fprintf(fid, '\\vspace{3mm}\n\\normalsize\n\\doublespacing\n');
fprintf(fid, '\\textit{Notes:} COMPLETE THE NOTES MANUALLY');  % Replace or complete this as needed
fprintf(fid, '\n\\end{table}\n');

% Close file
fclose(fid);

disp('LaTeX table written to Table_3.tex');