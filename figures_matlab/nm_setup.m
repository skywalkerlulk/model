function NM = nm_setup()
% nm_setup  Nature 风格图表统一设置 (字体 / 配色 / 线宽)
% 用法: 在每个绘图脚本开头调用  NM = nm_setup();
    if ismac
        NM.font = 'PingFang SC';
    elseif ispc
        NM.font = 'Microsoft YaHei';
    else
        NM.font = 'DejaVu Sans';
    end
    % 克制配色 (Nature 风格: 低饱和, 可区分)
    NM.palette = [ ...
        12  123 220;   % 1 蓝
        214  69  65;   % 2 红
        60  180 120;   % 3 绿
        140 120 200;   % 4 紫
        200 170  60;   % 5 金
        120 120 120] / 255;   % 6 灰
    NM.fs  = 8;        % 轴刻度字号
    NM.fsl = 9;        % 轴标签字号
    NM.lw  = 1.1;      % 数据线宽
    NM.alw = 0.6;      % 轴线宽
    set(0, 'DefaultAxesFontName', NM.font, 'DefaultTextFontName', NM.font, ...
           'DefaultAxesFontSize', NM.fs, 'DefaultAxesLineWidth', NM.alw, ...
           'DefaultLineLineWidth', NM.lw, 'DefaultFigureColor', 'w');
end
