% This program plots the figure for presentation in paper, based on dataset
% "error_space.mat", which corresponds to initial: Alge + (s,p,gamma,d) =
% (0.6,1,-1,2).

clear
load('error_space.mat'); figure

% Data
Line = 0:512;
loglog(NNN,error_L2_T,'ro','MarkerSize', 8,'LineWidth', 1.5), hold on,
loglog(NNN,error_Linf_T,'bp','MarkerSize',10,'LineWidth', 1), hold on,
loglog(Line,10^(-0.45)*Line.^(-3.29),'r--','LineWidth',2), hold on
loglog(Line,10^(-0.11)*Line.^(-4.80),'b--','LineWidth',2)

% Label
set(gca, 'FontSize', 16); % size of tick
xlabel('$N$','Interpreter','latex','FontSize',20); % xlabel
ylabel('Error','FontSize',20,'FontName','Times New Roman'); % ylabel
set(gca,'XTick',2.^(1:7));
set(gca,'XTickLabel',{'2','4','8','16','32','64','128'});
set(gca,'yTick',10.^(-10:2:-0));
set(gca,'yTickLabel',{'10^{-10}','10^{-8}','10^{-6}','10^{-4}','10^{-2}', ...
    '10^{0}'});
set(gca,'TickLabelInterpreter','tex');
axis([2^2,2^7,10^-10.25,10^0])

% Legend
h_legend = legend('$L^2$-error','$L^\infty$-error','$\propto N^{-3.29}$', ...
    '$\propto N^{-4.80}$','Interpreter', 'latex');
set(h_legend, 'FontSize', 20, 'FontName', 'Times New Roman', ... 
    'FontWeight', 'normal', 'TextColor', [0, 0, 0]);
set(h_legend,'Position', [0.767, 0.766, 0.0, 0.0]);