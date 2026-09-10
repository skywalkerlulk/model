function panel_label(ax, txt)
% panel_label 在面板左上角标注 Nature 风格小写加粗字母 (a, b, c, ...)
    pos = ax.Position;                       % 归一化坐标
    fig = ax.Parent;
    x = pos(1) - 0.05;
    y = pos(2) + pos(4) + 0.012;
    annotation(fig, 'textbox', [x y 0.06 0.06], 'String', txt, ...
        'EdgeColor', 'none', 'FontWeight', 'bold', 'FontSize', 11, ...
        'FontName', get(ax, 'FontName'), 'FitBoxToText', 'on', ...
        'VerticalAlignment', 'middle');
end
