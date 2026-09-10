function fig3_q3()
% fig3_q3  问题3: 烘干时长判定 (三维曲面 + 0.15等值线 + 结束平面 + 判据跨越放大图)
    NM = nm_setup();
    d = load('results/matlab/q2.mat');
    t_end = 207960 / 3600;                      % h

    fig = figure('Units', 'centimeters', 'Position', [2 2 16.5 6.8], 'Color', 'w');
    tiled = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    % ---- a) 全流程水分浓度三维曲面 (surfc: 曲面 + 底部等高线投影) ----
    ax = nexttile;
    step_t = 4; step_r = 2;
    r4 = d.r_cm(1:step_r:end);
    t4 = d.t60(1:step_t:end) / 3600;
    Z4 = d.C60(1:step_t:end, 1:step_r:end);
    surfc(ax, r4, t4, Z4, 'EdgeColor', 'none');
    shading(ax, 'interp');
    colormap(ax, parula);
    hold(ax, 'on');
    % 0.15 等值线 (贴合曲面)
    contour3(ax, r4, t4, Z4, [0.15 0.15], 'r', 'LineWidth', 1.2);
    % 烘干结束时刻半透明平面
    [Rp, Zp] = meshgrid(linspace(0, 2, 2), linspace(0, max(Z4(:)), 8));
    surf(ax, Rp, t_end * ones(size(Zp)), Zp, 'FaceColor', [0.1 0.55 0.3], ...
        'FaceAlpha', 0.15, 'EdgeColor', 'none');
    hold(ax, 'off');
    light('Position', [-0.6 -0.4 1.0], 'Style', 'infinite');
    lighting(ax, 'gouraud'); material(ax, 'dull');
    view(ax, [-42 22]);
    xlabel(ax, 'r / cm'); ylabel(ax, 't / h'); zlabel(ax, 'C / (kg·kg^{-1})');
    ax.FontSize = NM.fs;
    cb = colorbar(ax); cb.Label.String = 'C / (kg·kg^{-1})'; cb.FontSize = NM.fs;
    panel_label(ax, 'a');

    % ---- b) 中心/表面历程 + 判据跨越 ----
    ax = nexttile; hold(ax, 'on');
    plot(ax, d.t60 / 3600, d.C60(:, 1),   'Color', NM.palette(1, :), 'LineWidth', 1.4);
    plot(ax, d.t60 / 3600, d.C60(:, end), 'Color', NM.palette(2, :), 'LineWidth', 1.4);
    plot(ax, [0 t_end], [0.15 0.15], '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 0.9);
    plot(ax, [t_end t_end], [0 2.55], ':', 'Color', [0.1 0.55 0.3], 'LineWidth', 1.2);
    style_ax(ax, 't / h', 'C / (kg·kg^{-1})');
    xlim(ax, [0 72]); ylim(ax, [0 2.6]);
    legend(ax, {'中心 r = 0', '表面 r = 2 cm', '判据 C = 0.15', ...
        sprintf('t_{end} = %.2f h', t_end)}, ...
        'Box', 'off', 'FontSize', 7, 'Location', 'northeast');
    % 判据跨越放大子图 (55-59h)
    ax2 = axes('Position', [ax.Position(1)+0.09*ax.Position(3), ...
                            ax.Position(2)+0.40*ax.Position(4), ...
                            0.30*ax.Position(3), 0.42*ax.Position(4)]);
    hold(ax2, 'on');
    ii = d.t60 >= 55*3600 & d.t60 <= 59*3600;
    plot(ax2, d.t60(ii)/3600, d.C60(ii, 1), 'Color', NM.palette(1, :), 'LineWidth', 1.3);
    plot(ax2, [t_end t_end], [0.148 0.152], ':', 'Color', [0.1 0.55 0.3], 'LineWidth', 1.2);
    plot(ax2, [55 59], [0.15 0.15], '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 0.8);
    ax2.FontSize = 6; ax2.Box = 'on'; ax2.TickDir = 'out';
    xlabel(ax2, 't / h'); ylabel(ax2, 'C(0,t)');
    title(ax2, '判据跨越 (60 s 分辨率)', 'FontSize', 6.5, 'FontWeight', 'normal');
    panel_label(ax, 'b');

    export_clean(fig, 'fig3_q3', 'figures/matlab');
    fprintf('fig3_q3 已输出\n');
end
