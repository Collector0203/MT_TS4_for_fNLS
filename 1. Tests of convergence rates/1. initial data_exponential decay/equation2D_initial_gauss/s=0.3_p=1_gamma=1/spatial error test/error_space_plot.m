% This program plots the figure for presentation in paper, based on dataset
% "error_space.mat", which corresponds to initial: Gauss + (s,p,gamma,d) =
% (0.3,1,1,2).

clear
load('error_space.mat'); figure

% Data
Line = 0:512;
loglog(NNN,error_L2_T,'ro','MarkerSize', 8,'LineWidth', 1.5), hold on,
loglog(NNN,error_Linf_T,'bp','MarkerSize',10,'LineWidth', 1), hold on,
loglog(Line,10^(-1.15)*Line.^(-2.41),'r--','LineWidth',2), hold on,
loglog(Line,10^(-1.1)*Line.^(-3.90),'b--','LineWidth',2)

% Label
set(gca,'FontSize',16); % size of tick
xlabel('$N$','Interpreter','latex','FontSize',20); % xlabel
ylabel('Error','FontSize',20,'FontName','Times New Roman'); % ylabel
set(gca,'XTick',2.^(1:7));
set(gca,'XTickLabel',{'2','4','8','16','32','64','128'});
set(gca,'yTick',10.^(-9:2:-1));
set(gca,'yTickLabel',{'10^{-9}','10^{-7}','10^{-5}','10^{-3}','10^{-1}'});
set(gca,'TickLabelInterpreter','tex');
axis([2^2,2^7,10^-9.375,10^-0.5])

% Legend
h_legend = legend('$L^2$-error','$L^\infty$-error','$\propto N^{-2.41}$','$\propto N^{-3.90}$','Interpreter', 'latex');
set(h_legend, 'FontSize', 20, 'FontName', 'Times New Roman', ... 
    'FontWeight', 'normal', 'TextColor', [0, 0, 0]);
set(h_legend,'Position', [0.767, 0.766, 0.0, 0.0]);