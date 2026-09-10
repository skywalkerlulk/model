function fig6_env()
% fig6_env  烘房环境条件与平台期检测 (附件1)
    NM = nm_setup();
    d = load('results/matlab/env.mat');

    fig = figure('Units', 'centimeters', 'Position', [2 2 16.5 5.6], 'Color', 'w');
    tiled = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    % ---- a) 烘房温度 ----
    ax = nexttile; hold(ax, 'on');
    tmin = min(d.T_env) - 1; tmax = max(d.T_env) + 1;
    patch(ax, [7200 14400 14400 7200], [tmin tmin tmax tmax], ...
        [0.6 0.6 0.6], 'FaceAlpha', 0.12, 'EdgeColor', 'none');
    % ±σ 波动带 (平台期统计证据)
    patch(ax, [7200 14400 14400 7200], ...
        [d.T_bar-d.T_std d.T_bar-d.T_std d.T_bar+d.T_std d.T_bar+d.T_std], ...
        NM.palette(2, :), 'FaceAlpha', 0.20, 'EdgeColor', 'none');
    plot(ax, d.t / 3600, d.T_env, '-', 'Color', NM.palette(1, :), 'LineWidth', 1.3);
    plot(ax, [7200 14400], [d.T_bar d.T_bar], '--', ...
        'Color', NM.palette(2, :), 'LineWidth', 1.0);
    plot(ax, [0.5 0.5], [tmin tmax], ':', 'Color', [0.2 0.2 0.2], 'LineWidth', 0.9);
    style_ax(ax, 't / h', '烘房温度 / °C');
    legend(ax, {'平台期 (t ≥ 2 h)', '±1σ 波动带', '附件1 实测', ...
        sprintf('均值 %.2f ± %.2f °C', d.T_bar, d.T_std), '问题1窗口 (0.5 h)'}, ...
        'Box', 'off', 'FontSize', 6.5, 'Location', 'southeast');
    panel_label(ax, 'a');

    % ---- b) 烘房水分浓度 ----
    ax = nexttile; hold(ax, 'on');
    cmin = min(d.C_env) - 0.001; cmax = max(d.C_env) + 0.001;
    patch(ax, [7200 14400 14400 7200], [cmin cmin cmax cmax], ...
        [0.6 0.6 0.6], 'FaceAlpha', 0.12, 'EdgeColor', 'none');
    patch(ax, [7200 14400 14400 7200], ...
        [d.C_bar-d.C_std d.C_bar-d.C_std d.C_bar+d.C_std d.C_bar+d.C_std], ...
        NM.palette(2, :), 'FaceAlpha', 0.20, 'EdgeColor', 'none');
    plot(ax, d.t / 3600, d.C_env, '-', 'Color', NM.palette(1, :), 'LineWidth', 1.3);
    plot(ax, [7200 14400], [d.C_bar d.C_bar], '--', ...
        'Color', NM.palette(2, :), 'LineWidth', 1.0);
    plot(ax, [0.5 0.5], [cmin cmax], ':', 'Color', [0.2 0.2 0.2], 'LineWidth', 0.9);
    style_ax(ax, 't / h', '烘房水分浓度 / (kg·kg^{-1})');
    legend(ax, {'平台期 (t ≥ 2 h)', '±1σ 波动带', '附件1 实测', ...
        sprintf('均值 %.4f ± %.4f', d.C_bar, d.C_std), '问题1窗口 (0.5 h)'}, ...
        'Box', 'off', 'FontSize', 6.5, 'Location', 'southeast');
    panel_label(ax, 'b');

    export_clean(fig, 'fig6_env', 'figures/matlab');
    fprintf('fig6_env 已输出\n');
end
