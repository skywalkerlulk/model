function fig4_q4()
% fig4_q4  问题4: 收缩干燥 (平台阴影 + Δt 定量标注 + 收缩-失水一致性)
    NM = nm_setup();
    d = load('results/matlab/q4.mat');
    q2 = load('results/matlab/q2.mat');
    t_end4 = 184200 / 3600;                     % h
    t_end3 = 207960 / 3600;

    fig = figure('Units', 'centimeters', 'Position', [2 2 17.5 6.2], 'Color', 'w');
    tiled = tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

    % ---- a) 半径收缩曲线 (附件2) + 收缩平台阴影 ----
    ax = nexttile; hold(ax, 'on');
    patch(ax, [252000/3600 72 72 252000/3600], [1.1 1.1 2.1 2.1], ...
        [0.6 0.6 0.6], 'FaceAlpha', 0.12, 'EdgeColor', 'none');
    plot(ax, d.t_R / 3600, d.R_cm, 'o', 'MarkerSize', 2.5, ...
        'MarkerFaceColor', NM.palette(1, :), 'MarkerEdgeColor', 'none');
    plot(ax, d.t / 3600, d.R * 100, '-', 'Color', [0.75 0.35 0.35], 'LineWidth', 1.0);
    plot(ax, [70 70], [1.1 2.1], ':', 'Color', [0.2 0.2 0.2], 'LineWidth', 0.9);
    plot(ax, [t_end4 t_end4], [1.1 2.1], '--', 'Color', [0.1 0.55 0.3], 'LineWidth', 1.1);
    text(ax, 70.4, 1.16, '收缩平台', 'FontSize', 6, 'Color', [0.4 0.4 0.4]);
    style_ax(ax, 't / h', 'R / cm');
    ylim(ax, [1.1 2.1]);
    legend(ax, {'附件2 实测', '模型 PCHIP 插值', '收缩平台 70 h', ...
        sprintf('t_{end} = %.2f h', t_end4)}, ...
        'Box', 'off', 'FontSize', 6, 'Location', 'northeast');
    panel_label(ax, 'a');

    % ---- b) 问题3 vs 问题4 + Δt 双箭头定量标注 ----
    ax = nexttile; hold(ax, 'on');
    plot(ax, q2.t60 / 3600, q2.C60(:, 1), 'Color', NM.palette(4, :), 'LineWidth', 1.3);
    plot(ax, d.t / 3600, d.C(:, 1), 'Color', NM.palette(1, :), 'LineWidth', 1.3);
    plot(ax, [0 72], [0.15 0.15], '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 0.9);
    plot(ax, [t_end3 t_end3], [0 2.6], ':', 'Color', NM.palette(4, :), 'LineWidth', 1.0);
    plot(ax, [t_end4 t_end4], [0 2.6], ':', 'Color', NM.palette(1, :), 'LineWidth', 1.0);
    style_ax(ax, 't / h', 'C(0,t) / (kg·kg^{-1})');
    xlim(ax, [0 72]); ylim(ax, [0 2.6]);
    % 双箭头标注 (数据坐标 → 图归一化坐标)
    axp = ax.Position;
    xn = @(x) axp(1) + axp(3) * x / 72;
    yn = @(y) axp(2) + axp(4) * y / 2.6;
    annotation(fig, 'doublearrow', [xn(t_end4) yn(2.30) xn(t_end3) yn(2.30)], ...
        'Color', [0.2 0.2 0.2], 'LineWidth', 0.9);
    text(ax, (t_end3 + t_end4) / 2, 2.36, 'Δ = 6.60 h (11.4%)', ...
        'FontSize', 6.5, 'HorizontalAlignment', 'center', 'Color', [0.2 0.2 0.2]);
    legend(ax, {'问题3 固定半径', '问题4 收缩模型', '判据 0.15'}, ...
        'Box', 'off', 'FontSize', 6.5, 'Location', 'northeast');
    panel_label(ax, 'b');

    % ---- c) 收缩-失水一致性检验 ----
    ax = nexttile; hold(ax, 'on');
    xi = d.xi(:)';
    dxi = 1 / 400;
    xif = ((1:400) - 0.5) * dxi;
    w = zeros(1, 401);
    w(1) = xif(1)^2;
    w(2:400) = xif(2:400).^2 - xif(1:399).^2;
    w(401) = 1 - xif(400)^2;
    Mw = d.C .* w;
    wl = 1 - sum(Mw, 2) / sum(Mw(1, :));       % 模型失水分数
    vs = 1 - (d.R / 0.02).^2;                  % 体积收缩分数
    plot(ax, wl, vs, 'Color', NM.palette(1, :), 'LineWidth', 1.4);
    plot(ax, [0 1], [0 1], '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 0.9);
    text(ax, 0.16, 0.82, '早期线性相关', 'FontSize', 6.5, 'Color', [0.2 0.2 0.2]);
    text(ax, 0.72, 0.55, '后期收缩饱和', 'FontSize', 6.5, 'Color', [0.2 0.2 0.2]);
    style_ax(ax, '模型失水分数 1 - M_w / M_{w0}', '体积收缩分数 1 - (R/R_0)^2');
    legend(ax, {'收缩-失水轨迹', '理想收缩线 (1:1)'}, ...
        'Box', 'off', 'FontSize', 6.5, 'Location', 'southeast');
    panel_label(ax, 'c');

    export_clean(fig, 'fig4_q4', 'figures/matlab');
    fprintf('fig4_q4 已输出\n');
end
