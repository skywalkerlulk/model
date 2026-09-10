function style_ax(ax, xl, yl)
% style_ax 二维坐标轴 Nature 风格微调: 细网格 / 外向刻度 / 标签
    xlabel(ax, xl); ylabel(ax, yl);
    ax.Box = 'on';
    ax.XGrid = 'on'; ax.YGrid = 'on';
    ax.GridAlpha = 0.12; ax.GridColor = [0.35 0.35 0.35];
    ax.TickDir = 'out'; ax.TickLength = [0.004 0.004];
end
