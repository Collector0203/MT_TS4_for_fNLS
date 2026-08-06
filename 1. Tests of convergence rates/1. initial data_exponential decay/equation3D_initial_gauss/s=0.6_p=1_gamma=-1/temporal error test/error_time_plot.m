% This program plots the figure for presentation in paper, based on dataset
% "error_time.mat", which corresponds to initial: Gauss + (s,p,gamma,d) =
% (0.6,1,-1,3).

clear
load('error_time.mat');

% Data
loglog(Dt,error_L2_T,'ro','MarkerSize', 8,'LineWidth', 1.5), hold on,
loglog(Dt,error_Linf_T,'bp','MarkerSize',10,'LineWidth', 1), hold on, 
loglog(Dt,Dt.^4, '--','Color',[0.3,0.3,0.3],'LineWidth',2);

% Label
set(gca, 'FontSize', 16); % size of tick
xlabel('$\tau$','Interpreter','latex','FontSize',20); % xlabel
ylabel('Error','FontSize',20,'FontName','Times New Roman'); % ylabel
set(gca,'XTick',flip(Dt));
set(gca,'XTickLabel',{'2^{-6}','2^{-5.5}','2^{-5}','2^{-4.5}','2^{-4}','2^{-3.5}','2^{-3}'});
set(gca,'yTick',10.^(-8:-2));
set(gca,'yTickLabel',{'10^{-8}','10^{-7}','10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}'});
set(gca,'TickLabelInterpreter','tex');
axis([2^-6.3,2^-2.7,10^-8,10^-2])

% Legend
h_legend = legend('$L^2$-error','$L^\infty$-error','$\tau^4$','Interpreter', 'latex');
set(h_legend, 'FontSize', 20, 'FontName', 'Times New Roman', ... 
    'FontWeight', 'normal', 'TextColor', [0, 0, 0]);
current_pos = h_legend.Position;
new_pos = current_pos + [-0.427, -0.018, 0, 0]; % position of legend
h_legend.Position = new_pos;