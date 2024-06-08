%% Auxiliary plot function

function plot_band(x, y, lower, upper, plot_title, plot_ylim, plot_xtick)

    % Plot IRF and error bands

    plot(x, y, '-k', 'LineWidth', 2);
    hold on;
    plot(x, [lower; upper], 'LineStyle', '--', 'Color', 'k');
    line([min(x) max(x)], [0 0], 'Color', 'k');
    hold off;
    
    xlim([min(x) max(x)]);
    ylim(plot_ylim);
    set(gca, 'XTick', plot_xtick);
    set(gca, 'FontSize', 12);
    title(plot_title, 'FontSize', 14);

end