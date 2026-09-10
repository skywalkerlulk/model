% run_all_figs.m — 一键生成全部 Nature 风格图 (600 dpi PNG + 矢量 PDF)
% 用法: 在 MATLAB 中直接运行本脚本 (或在命令行 matlab -batch run_all_figs)
% 输出: figures/matlab/ 下 fig1_q1 ~ fig6_env (png + pdf)

prj = fileparts(fileparts(mfilename('fullpath')));   % 项目根
cd(prj);
addpath(genpath(fileparts(mfilename('fullpath'))));  % figures_matlab/ 加入路径

fig1_q1();
fig2_q2();
fig3_q3();
fig4_q4();
fig5_validation();
fig6_env();
close all;
fprintf('全部图片已输出至 figures/matlab/\n');
