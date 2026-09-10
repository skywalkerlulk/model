function export_clean(fig, name, outdir)
% export_clean 统一导出: 600 dpi PNG + 矢量 PDF
% 兼容 R2019b (print) 与 R2020a+ (exportgraphics)
    if ~exist(outdir, 'dir'), mkdir(outdir); end
    png_path = fullfile(outdir, [name '.png']);
    pdf_path = fullfile(outdir, [name '.pdf']);
    try
        exportgraphics(fig, png_path, 'Resolution', 600);
    catch
        print(fig, png_path, '-dpng', '-r600');
    end
    try
        exportgraphics(fig, pdf_path, 'ContentType', 'vector');
    catch
        print(fig, pdf_path, '-dpdf', '-r300');
    end
end
