% This program plots the discrete energy deviation at T = 25 as the time
% step decreases for equation parameters (s,p,gamma) = (0.3,1,1).

clear
load('energy_order_(0.3,1,1).mat');

% Data
loglog(Dt,energy_drift_1,'rd','MarkerSize', 8,'LineWidth', 1.5), hold on,
loglog(Dt,energy_drift_2,'bs','MarkerSize',10,'LineWidth', 1.5), hold on, 
loglog(Dt,10^-2.5*Dt.^4, '--','Color',[0.3,0.3,0.3],'LineWidth',2);

% Label
set(gca, 'FontSize', 16); % size of tick
xlabel('$\tau$','Interpreter','latex','FontSize',20); % xlabel
ylabel('$\mathcal{E}_E(T)$','Interpreter','latex','FontSize',20); % ylabel
set(gca,'XTick',flip(Dt));
set(gca,'XTickLabel',{'2^{-6.5}','2^{-6}','2^{-5.5}','2^{-5}','2^{-4.5}','2^{-4}','2^{-3.5}'});
set(gca,'yTick',10.^(-11:1:-6));
set(gca,'TickLabelInterpreter','tex');
axis([2^-6.8,2^-3.8,10^-11,10^-6.5])

% Legend
h_legend = legend('initial data 1','initial data 2','$10^{-2.5}\tau^4$','Interpreter', 'latex');
set(h_legend, 'FontSize', 20, 'FontName', 'Times New Roman', ... 
    'FontWeight', 'normal', 'TextColor', [0, 0, 0]);
current_pos = h_legend.Position;
new_pos = current_pos + [-0.327, -0.015, 0, 0]; % position of legend
h_legend.Position = new_pos;