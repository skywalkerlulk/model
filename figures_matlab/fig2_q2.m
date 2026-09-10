function fig2_q2()
% fig2_q2  问题2: 变物性耦合模型 (3h 剖面 + 72h 历程 + 多温度 D-C 机理曲线)
    NM = nm_setup();
    d = load('results/matlab/q2.mat');
    t_end3 = 207960 / 3600;

    fig = figure('Units', 'centimeters', 'Position', [2 2 16.5 11.5], 'Color', 'w');
    tiled = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nT = size(d.T3, 1);
    iT = interp1(d.t3, 1:nT, 1800:1800:10800, 'nearest');
    lab = arrayfun(@(x) sprintf('%.1f h', x), d.t3(iT) / 3600, 'UniformOutput', false);

    % ---- a) 温度剖面 (0.5-3h) ----
    ax = nexttile; hold(ax, 'on');
    for k = 1:numel(iT)
        plot(ax, d.r_out_cm, d.T3(iT(k), :), 'Color', NM.palette(mod(k-1, 6)+1, :));
    end
    style_ax(ax, 'r / cm', 'T / °C');
    legend(ax, lab, 'Box', 'off', 'FontSize', 6.5, 'Location', 'southeast');
    panel_label(ax, 'a');

    % ---- b) 水分浓度剖面 (0.5-3h) ----
    ax = nexttile; hold(ax, 'on');
    for k = 1:numel(iT)
        plot(ax, d.r_out_cm, d.C3(iT(k), :), 'Color', NM.palette(mod(k-1, 6)+1, :));
    end
    style_ax(ax, 'r / cm', 'C / (kg·kg^{-1})');
    legend(ax, lab, 'Box', 'off', 'FontSize', 6.5, 'Location', 'northeast');
    panel_label(ax, 'b');

    % ---- c) 72h 中心/表面水分历程 + 烘干结束标记 ----
    ax = nexttile; hold(ax, 'on');
    plot(ax, d.t60 / 3600, d.C60(:, 1),   'Color', NM.palette(1, :), 'LineWidth', 1.4);
    plot(ax, d.t60 / 3600, d.C60(:, end), 'Color', NM.palette(2, :), 'LineWidth', 1.4);
    plot(ax, [0 72], [0.15 0.15], '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 0.8);
    plot(ax, [t_end3 t_end3], [0 2.55], ':', 'Color', [0.1 0.55 0.3], 'LineWidth', 1.1);
    style_ax(ax, 't / h', 'C / (kg·kg^{-1})');
    xlim(ax, [0 72]);
    legend(ax, {'中心 r = 0', '表面 r = 2 cm', '判据 0.15', ...
        sprintf('t_{end} = %.2f h', t_end3)}, ...
        'Box', 'off', 'FontSize', 7, 'Location', 'northeast');
    panel_label(ax, 'c');

    % ---- d) 多温度 D-C 机理曲线 (附录3, Arrhenius 项) ----
    ax = nexttile; hold(ax, 'on');
    Cg = linspace(0.05, 2.55, 300);
    Tc_list = [28 40 50];
    c_idx = [5 4 1];
    for k = 1:numel(Tc_list)
        Dg = 2.4e-3 * exp(-0.45 ./ Cg) * exp(-3850 / (Tc_list(k) + 273.15));
        semilogy(ax, Cg, Dg, 'Color', NM.palette(c_idx(k), :), 'LineWidth', 1.3, ...
            'DisplayName', sprintf('T = %.0f °C', Tc_list(k)));
    end
    plot(ax, [0.15 0.15], [1e-12 1e-7], '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 0.8);
    plot(ax, [2.55 2.55], [1e-12 1e-7], ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.8);
    text(ax, 0.30, 3e-8, '判据 C = 0.15', 'FontSize', 6.5, ...
        'Color', [0.2 0.2 0.2], 'Rotation', 90);
    text(ax, 2.28, 3e-8, '湿态冻结点', 'FontSize', 6.5, 'Color', [0.4 0.4 0.4]);
    style_ax(ax, 'C / (kg·kg^{-1})', 'D / (m^2·s^{-1})');
    legend(ax, 'Box', 'off', 'FontSize', 6.5, 'Location', 'southwest');
    panel_label(ax, 'd');

    export_clean(fig, 'fig2_q2', 'figures/matlab');
    fprintf('fig2_q2 已输出\n');
end
