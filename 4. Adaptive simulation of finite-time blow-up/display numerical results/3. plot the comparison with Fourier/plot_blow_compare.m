% This program plots the comparison between BioMT+TSM4 and Fourier+CRK43
% regarding the blow-up problem. The data were obtained using programs
% equation2D_blow_gauss.m and Fourier_CRK43_fNLS2d.m, respectively.

clear, close all
figure('Position', [0, 200, 1800, 600])

%% plot the L-inf norm versus time
load('BioMT_TSM4_s=0.7_p=1.5_N=256_Bound=20.5.mat');
U_red = U_max; T_red = TT;
h_MT = plot(U_max, TT, 'r-', 'LineWidth', 3); hold on,

t_marker = 1.3; % "t_marker" represents the marker point used to distinguish the four curves corresponding to the Fourier method

load('Fourier_CRK43_s=0.7_p=1.5_N=256.mat');
idx_marker = find(TT <= t_marker, 1, 'last');
h_F256 = plot(U_max, TT, '-o', ...
    'LineWidth', 3, ...
    'Color', [0.20, 0.45, 0.85], ...
    'MarkerIndices', idx_marker, ...
    'MarkerSize', 9);
hold on,

load('Fourier_CRK43_s=0.7_p=1.5_N=512.mat');
idx_marker = find(TT <= t_marker, 1, 'last');
h_F512 = plot(U_max, TT, '-s', ...
    'LineWidth', 3, ...
    'Color', [0.80, 0.5, 0.05], ...
    'MarkerIndices', idx_marker, ...
    'MarkerSize', 9);
hold on,

load('Fourier_CRK43_s=0.7_p=1.5_N=1024.mat');
idx_marker = find(TT <= t_marker, 1, 'last');
h_F1024 = plot(U_max, TT, '-^', ...
    'LineWidth', 3, ...
    'Color', [0.57, 0.30, 0.63], ...
    'MarkerIndices', idx_marker, ...
    'MarkerSize', 9);
hold on,

load('Fourier_CRK43_s=0.7_p=1.5_N=2048.mat');
idx_marker = find(TT <= t_marker, 1, 'last');
h_F2048 = plot(U_max, TT, '-d', ...
    'LineWidth', 3, ...
    'Color', [0, 0.5, 0.5], ...
    'MarkerIndices', idx_marker, ...
    'MarkerSize', 9);

%% blow-up moment
Tstar = 1.27405;
yline(Tstar, 'k--', 'LineWidth', 3); 
text(16.7, 1.285, '$T^* \approx 1.27405$', 'Interpreter', 'latex', 'FontSize', 26)

%% label and axis
xticks = 2:2:20; set(gca,'XTick',xticks); set(gca,'XTickLabel',num2str(xticks'),'FontSize',18);
yticks = 1.2:0.01:1.3; set(gca,'YTick',yticks); set(gca,'YTickLabel',num2str(yticks'),'FontSize',18);
xlabel('$\Vert \psi(\cdot,t) \Vert_{L^\infty(\mathrm{R}^2)}$','Interpreter','latex','FontSize',26);
ylabel('$t$','Interpreter','latex','FontSize',26);
xlim([1.5, 20.5]), ylim([1.20, 1.30])

%% legend
% h_legend = legend('MT + TS4, $N=2^8$', 'Fourier + CRK43, $N=2^8$', 'Fourier + CRK43, $N=2^9$', 'Fourier + CRK43, $N=2^{10}$', 'Fourier + CRK43, $N=2^{11}$', 'Interpreter', 'latex');
h_legend = legend('MT--TS4, $N=256$', 'Fourier--CRK43, $N=256$', 'Fourier--CRK43, $N=512$', 'Fourier--CRK43, $N=1024$', 'Fourier--CRK43, $N=2048$', 'Interpreter', 'latex');
set(h_legend, 'FontSize', 24, 'FontName', 'Times New Roman', ... 
    'FontWeight', 'normal', 'TextColor', [0, 0, 0]);
current_pos = h_legend.Position;
new_pos = current_pos + [-0.34, -0.32, 0, 0]; % position of legend
h_legend.Position = new_pos;

%% plot a rectangle representing the magnified area
x0 = 20; y0 = Tstar; % the center point of the rectangle
a = 0.3; b = 0.003; % the length and width of the rectangle
rectangle('Position', [x0-a/2, y0-b/2, a, b],'EdgeColor', 'k', ...
    'LineStyle', '-','LineWidth', 0.6);

%% plot an arrow to indicate magnification
ax = gca;
drawnow;

% Arrow start and end points: data coordinates (箭头起点和终点：数据坐标)
x_start = x0; y_start = y0 - b/2;
% x_end = 16.86; y_end = 1.26;
x_end = 19.50; y_end = 1.26;

% Get axes position in figure (获取坐标轴在 figure 中的位置)
oldUnits = ax.Units;
ax.Units = 'normalized';
pos = ax.Position;    
ax.Units = oldUnits;

xl = ax.XLim; yl = ax.YLim;

% Convert data coordinates to figure normalized coordinates (数据坐标转换为 figure normalized 坐标)
x_start_n = pos(1) + (x_start-xl(1))/(xl(2)-xl(1))*pos(3);
y_start_n = pos(2) + (y_start-yl(1))/(yl(2)-yl(1))*pos(4);

x_end_n = pos(1) + (x_end-xl(1))/(xl(2)-xl(1))*pos(3);
y_end_n = pos(2) + (y_end-yl(1))/(yl(2)-yl(1))*pos(4);

% Draw text arrow from start to end, arrowhead at the second point (从起点指向终点，箭头尖端位于第二个点)
annotation(gcf, 'textarrow', [x_start_n, x_end_n], [y_start_n, y_end_n], ... 
    'String', '', 'Color', 'k', 'LineWidth', 0.6, ... 
    'HeadStyle', 'vback2', 'HeadLength', 10, 'HeadWidth', 10);

%% plot a magnified inset below the arrow
% Remove possible repeated U-values and arrange them in ascending order
[U_red_unique, ind_unique] = unique(U_red(:), 'sorted');
T_red_unique = T_red(ind_unique);

% Red curve value at x = x0
T_red_x0 = interp1(U_red_unique, T_red_unique, x0, 'linear');

% Position of the inset axes in normalized figure coordinates [left, bottom, width, height]
ax_zoom = axes('Position', [0.65, 0.262, 0.215, 0.345]);
hold(ax_zoom, 'on'), box(ax_zoom, 'on')

% Plot the local red curve
plot(ax_zoom, U_red, T_red, 'r-', 'LineWidth', 3);

% Plot the blow-up time
yline(ax_zoom, Tstar, 'k--', 'LineWidth', 3);

% Use the rectangle's horizontal range
xlim(ax_zoom, [x0-a/2, x0+a/2]);

% Determine a suitable vertical range automatically
local_ind = U_red >= x0-a/2 & U_red <= x0+a/2;
T_local = T_red(local_ind);

ymin_zoom = min([T_local(:); T_red_x0; Tstar]);
ymax_zoom = max([T_local(:); T_red_x0; Tstar]);

% Add vertical padding; also ensure that the true gap is visible
gap = Tstar - T_red_x0;
ypad = max(0.8*abs(gap), 0.10*(ymax_zoom-ymin_zoom));
ylim(ax_zoom, [ymin_zoom-ypad, ymax_zoom+ypad]);

% Appearance of the inset
set(ax_zoom, 'FontSize', 15, 'LineWidth', 0.6, 'XTick', [x0-a/2, x0, x0+a/2]);
grid(ax_zoom, 'off')

%% indicate the gap by a vertical dimension line
y_lower  = min(T_red_x0, Tstar);
y_upper  = max(T_red_x0, Tstar);
y_middle = (y_lower + y_upper)/2;

% Put the dimension line slightly to the right of x = 20
x_gap = x0;

% Vertical line with horizontal caps
h_gap = errorbar(ax_zoom, x_gap, y_middle, ...
    y_middle-y_lower-0.000007, y_upper-y_middle-0.000007, ...
    'k', ...
    'LineStyle', 'none', ...
    'LineWidth', 1, ...
    'CapSize', 12);

% Actual gap label
gap_string = sprintf( ...
    '$\\approx %.2f\\times 10^{-4}$', abs(gap)/1e-4);

text(ax_zoom, x_gap + 0.011, y_middle, gap_string, ...
    'Interpreter', 'latex', ...
    'FontSize', 20, ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'middle');