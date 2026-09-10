function fig1_q1()
% fig1_q1  问题1: 预热阶段温度/水分时空演化 (三维网格 + 环境驱动线 + 剖面曲线)
    NM = nm_setup();
    d = load('results/matlab/q1.mat');
    e = load('results/matlab/env.mat');
    r = d.r_cm; t = d.t1;                       % t: 1..1800 s, r: 401 节点
    idx_t = interp1(t, 1:numel(t), d.TABLE_TIMES, 'nearest');

    fig = figure('Units', 'centimeters', 'Position', [2 2 16.5 11.5], 'Color', 'w');
    tiled = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    step = 5;
    t_m = t(1:step:end) / 60;                   % 时间线降采样 1800→360 条

    % ---- a) 温度三维网格 (按高度着色) + 烘房温度驱动线 ----
    ax = nexttile;
    mesh(ax, r, t_m, d.T1(1:step:end, :), 'EdgeColor', [0.28 0.36 0.52], ...
        'EdgeAlpha', 0.40, 'LineWidth', 0.20);
    colormap(ax, parula); view(ax, [-38 26]);
    hold(ax, 'on');
    T_env_i = interp1(e.t, e.T_env, t(1:step:end), 'pchip');
    plot3(ax, 2 * ones(size(t_m)), t_m, T_env_i, 'Color', NM.palette(2, :), ...
        'LineWidth', 1.4);
    hold(ax, 'off');
    xlabel(ax, 'r / cm'); ylabel(ax, 't / min'); zlabel(ax, 'T / °C');
    ax.FontSize = NM.fs;
    cb = colorbar(ax); cb.Label.String = 'T / °C'; cb.FontSize = NM.fs;
    legend(ax, {'烘房温度 T_{env}(t)'}, 'Box', 'off', 'FontSize', 6.5, ...
        'Location', 'northwest');
    panel_label(ax, 'a');

    % ---- b) 水分浓度三维网格 ----
    ax = nexttile;
    mesh(ax, r, t_m, d.C1(1:step:end, :), 'EdgeColor', [0.30 0.50 0.42], ...
        'EdgeAlpha', 0.40, 'LineWidth', 0.20);
    colormap(ax, parula); view(ax, [-38 26]);
    xlabel(ax, 'r / cm'); ylabel(ax, 't / min'); zlabel(ax, 'C / (kg·kg^{-1})');
    ax.FontSize = NM.fs;
    cb = colorbar(ax); cb.Label.String = 'C / (kg·kg^{-1})'; cb.FontSize = NM.fs;
    panel_label(ax, 'b');

    % ---- c) 温度剖面曲线 ----
    ax = nexttile; hold(ax, 'on');
    for k = 1:numel(d.TABLE_TIMES)
        plot(ax, r, d.T1(idx_t(k), :), 'Color', NM.palette(mod(k-1, 6)+1, :));
    end
    style_ax(ax, 'r / cm', 'T / °C');
    legend(ax, arrayfun(@(x) sprintf('%d s', x), d.TABLE_TIMES, 'UniformOutput', false), ...
        'Box', 'off', 'FontSize', 6.5, 'Location', 'southeast');
    panel_label(ax, 'c');

    % ---- d) 水分浓度剖面曲线 ----
    ax = nexttile; hold(ax, 'on');
    for k = 1:numel(d.TABLE_TIMES)
        plot(ax, r, d.C1(idx_t(k), :), 'Color', NM.palette(mod(k-1, 6)+1, :));
    end
    style_ax(ax, 'r / cm', 'C / (kg·kg^{-1})');
    legend(ax, arrayfun(@(x) sprintf('%d s', x), d.TABLE_TIMES, 'UniformOutput', false), ...
        'Box', 'off', 'FontSize', 6.5, 'Location', 'northeast');
    panel_label(ax, 'd');

    export_clean(fig, 'fig1_q1', 'figures/matlab');
    fprintf('fig1_q1 已输出\n');
end
