% This program plots the figure for presentation in paper, based on dataset
% "error_space.mat", which corresponds to initial: Alge + (s,p,gamma,d) =
% (0.3,1,1,1).

clear
load('error_space.mat'); figure

% Data
Line = 0:512;
loglog(NNN,error_L2_T,'ro','MarkerSize', 8,'LineWidth', 1.5), hold on,
loglog(NNN,error_Linf_T,'bp','MarkerSize',10,'LineWidth', 1), hold on,
loglog(Line,10^(-1.05)*Line.^(-1.68),'r--','LineWidth',2), hold on
loglog(Line,10^(-1.15)*Line.^(-2.46),'b--','LineWidth',2)

% Label
set(gca, 'FontSize', 16); % size of tick
xlabel('$N$','Interpreter','latex','FontSize',20); % xlabel
ylabel('Error','FontSize',20,'FontName','Times New Roman'); % ylabel
set(gca,'XTick',2.^(1:8));
set(gca,'XTickLabel',{'2','4','8','16','32','64','128','256'});
set(gca,'yTick',10.^(-8:0));
set(gca,'yTickLabel',{'10^{-8}','10^{-7}','10^{-6}','10^{-5}','10^{-4}', ...
    '10^{-3}','10^{-2}','10^{-1}','10^0'});
set(gca,'TickLabelInterpreter','tex');
axis([2^1,2^8,10^-7.1,10^0])

% Legend
h_legend = legend('$L^2$-error','$L^\infty$-error','$\propto N^{-1.68}$', ...
    '$\propto N^{-2.46}$','Interpreter', 'latex');
set(h_legend, 'FontSize', 20, 'FontName', 'Times New Roman', ... 
    'FontWeight', 'normal', 'TextColor', [0, 0, 0]);
set(h_legend,'Position', [0.767, 0.766, 0.0, 0.0]);