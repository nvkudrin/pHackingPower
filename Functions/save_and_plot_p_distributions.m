function save_and_plot_p_distributions(P0, P1, P1min, fig_folder, csv_folder, tag, j, k)
    
    
    % Save CSVs
    csvwrite(fullfile(csv_folder, ['P0_', tag, '.csv']), P0);
    csvwrite(fullfile(csv_folder, ['P1_', tag, '.csv']), P1);
    csvwrite(fullfile(csv_folder, ['P1min_', tag, '.csv']), P1min);

    % Plot distributions (invisible)
    fig_handle = figure('Visible','off');
    hold on
    [centers0, N0] = Edges_N(P0);
    [centers1, N1] = Edges_N(P1);
    [centers1min, N1min] = Edges_N(P1min);
    AvePowerH0 = mean(P0 <= 0.05);

    plot(centers0, N0, '-',  'LineWidth', 2, 'Color', 'black');
    plot(centers1, N1, '--', 'LineWidth', 2, 'Color', 'black');
    plot(centers1min, N1min, ':', 'LineWidth', 2, 'Color', 'black');
    ylim([0 10]);
    if k == 3
        ylim([0 15]);
    end
    set(gca, 'FontSize', 18);
    xlabel('$$p$$', 'FontSize', 25, 'Interpreter', 'latex');
    ylabel('PDF', 'FontSize', 25, 'Interpreter', 'latex');
    if j < 4
        title(['Null and $p$-hacked distributions: $h = $ ', int2str(j - 1)], ...
            'FontWeight', 'bold', 'FontSize', 20, 'Interpreter', 'latex');
    else
        title('Null and $p$-hacked distributions: $h \sim \widehat\Pi(h)$', ...
            'FontWeight', 'bold', 'FontSize', 20, 'Interpreter', 'latex');
    end

    lgd = legend(["No $p$-hacking" + newline + ...
                 "(avg power: " + string(round(AvePowerH0 * 100)) + "$\%$)", ...
                 "Threshold", "Minimum"], ...
                 'Orientation', 'horizontal', ...
                 'Interpreter', 'latex', 'NumColumns', 3);
    lgd.FontSize = 16;

    % Inset (zoom) axis
    axes('position', [0.35 0.34 0.45 0.45]); box on;
    plot(centers0, N0, '-',  'LineWidth', 2, 'Color', 'black'); hold on
    plot(centers1, N1, '--', 'LineWidth', 2, 'Color', 'black');
    plot(centers1min, N1min, ':', 'LineWidth', 2, 'Color', 'black');
    xlim([0 0.1]);
    ylim([0 15]);
    if k == 3
        ylim([0 20]);
    end

    saveas(fig_handle, fullfile(fig_folder, ['P_', tag]), 'epsc');
    close(fig_handle);

    % --- Local function for computing histogram as PDF ---
    function [centers, N] = Edges_N(P)
    %EDGESN  Compute histogram for Monte Carlo p-curves
        %   [centers, N] = Edges_N(P)
        %   Computes a normalized histogram (PDF) of p-values using 400 bins on [0,1].
        %   Returns the bin centers and bin counts (N).
    
        edges = linspace(0, 1, 401);  % 400 bins between 0 and 1
        [N, edges] = histcounts(P, edges, 'Normalization', 'pdf');
    
        % Convert bin edges to bin centers
        centers = edges(2:end) - (edges(2) - edges(1))/2;
    end
end

