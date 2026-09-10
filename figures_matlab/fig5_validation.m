function fig5_validation()
% fig5_validation  数值验证 (解析解对比 + 网格收敛)
    NM = nm_setup();
    d = load('results/matlab/validation.mat');

    fig = figure('Units', 'centimeters', 'Position', [2 2 16.5 6.5], 'Color', 'w');
    tiled = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    % ---- a) 数值 vs Bessel 级数解析解 (阶跃边界) ----
    ax = nexttile; hold(ax, 'on');
    for k = 1:numel(d.t_targets)
        plot(ax, d.r_cm, d.T_num(k, :), 'o', 'MarkerSize', 3, ...
            'Color', NM.palette(k, :), 'LineWidth', 0.8, ...
            'MarkerFaceColor', 'none', 'DisplayName', sprintf('数值 t=%.0f s', d.t_targets(k)));
        plot(ax, d.r_cm, d.T_ref(k, :), '-', 'Color', NM.palette(k, :), ...
            'LineWidth', 1.2, 'DisplayName', sprintf('解析 t=%.0f s', d.t_targets(k)));
    end
    style_ax(ax, 'r / cm', 'T / °C');
    legend(ax, 'Box', 'off', 'FontSize', 6.5, 'Location', 'southeast');
    % 误差子图 (对数坐标)
    ax2 = axes('Position', [ax.Position(1)+0.10*ax.Position(3), ...
                            ax.Position(2)+0.10*ax.Position(4), ...
                            0.32*ax.Position(3), 0.32*ax.Position(4)]);
    hold(ax2, 'on');
    for k = 1:numel(d.t_targets)
        semilogy(ax2, d.r_cm, abs(d.T_num(k, :) - d.T_ref(k, :)), ...
            'Color', NM.palette(k, :), 'LineWidth', 0.9);
    end
    ax2.FontSize = 6; ax2.Box = 'on'; ax2.TickDir = 'out';
    ax2.XGrid = 'on'; ax2.YGrid = 'on'; ax2.GridAlpha = 0.12;
    xlabel(ax2, 'r / cm'); ylabel(ax2, '|误差| / °C');
    title(ax2, '误差分布', 'FontSize', 6.5, 'FontWeight', 'normal');
    panel_label(ax, 'a');

    % ---- b) 网格收敛 (真实问题1参数, t=1800s) ----
    ax = nexttile; hold(ax, 'on');
    loglog(ax, d.Ns, d.err_T, 'o-', 'Color', NM.palette(1, :), ...
        'MarkerSize', 4, 'MarkerFaceColor', NM.palette(1, :));
    loglog(ax, d.Ns, d.err_C, 's-', 'Color', NM.palette(2, :), ...
        'MarkerSize', 4, 'MarkerFaceColor', NM.palette(2, :));
    ref = 1e-4 * (400 ./ d.Ns).^2;
    loglog(ax, d.Ns, ref, 'k--', 'LineWidth', 0.9, 'DisplayName', '二阶参考线');
    style_ax(ax, '网格数 N', '相邻网格最大差');
    legend(ax, {'温度', '水分浓度', '二阶参考线'}, ...
        'Box', 'off', 'FontSize', 6.5, 'Location', 'southwest');
    text(ax, 0.97, 0.14, {'守恒残差 (Q1):', '水分 7.6×10^{-7}', ...
        '能量 7.8×10^{-8}', '时间收敛: dt减半差 < 10^{-7}'}, ...
        'Units', 'normalized', 'FontSize', 6.5, ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
        'BackgroundColor', 'w', 'EdgeColor', [0.6 0.6 0.6]);
    panel_label(ax, 'b');

    export_clean(fig, 'fig5_validation', 'figures/matlab');
    fprintf('fig5_validation 已输出\n');
end
