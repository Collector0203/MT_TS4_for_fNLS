% This program displays the data recorded in the adaptive process
% equation2D_blow_gauss.m, including: mass, energy, time step, scaling
% factor, L^inf-norm, H^1-seminorm.

clear
load("BioMT_TSM4_s=0.7_p=1.5_N=256_Bound=20.5.mat");
num_step = length(TT);

%% subfigure: Relative mass/energy deviation
figure
drift_mass = abs((Mass-Mass(1))/Mass(1));
drift_energy = abs((Energy-Energy(1))/Energy(1));

semilogy(drift_mass(2:num_step),'-','LineWidth', 1.5,'Color', [0, 0.8, 0],'Marker', 'd', ...
    'MarkerIndices', [1:124:2605, 2727],'MarkerSize', 8), hold on,
semilogy(drift_energy(2:num_step),'b-','LineWidth', 1.5,'Marker', 's', ...
    'MarkerIndices', [1:124:2605, 2727],'MarkerSize', 8)

% Label
set(gca,'FontSize',16); % size of tick
xlabel('Time-step index','FontSize',20,'FontName','Times New Roman') % xlabel
ylabel('Relative deviation','FontSize',20,'FontName','Times New Roman') % ylabel
set(gca,'xTick',200:500:2700);
set(gca,'yTick',10.^(-15:2:-1));
set(gca,'yTickLabel',{'10^{-15}','10^{-13}','10^{-11}', ...
    '10^{-9}','10^{-7}','10^{-5}','10^{-3}','10^{-1}'});
axis([1,num_step,10^-15,10^-1])

% Legend
h_legend = legend('mass deviation $\mathcal{E}_M $', ...
    'energy deviation $\mathcal{E}_E $','Interpreter', 'latex');
set(h_legend, 'FontSize', 20, 'FontName', 'Times New Roman', ...
    'FontWeight', 'normal', 'TextColor', [0, 0, 0]);
set(h_legend,'Position', [0.445, 0.805, 0.0, 0.0]);

%% subfigure: Time step & Scaling factor
figure
dt = diff(TT); % dt(m) denotes the time step used at the m-th time step
mu  = Mu(1:num_step-1); % mu(m) denotes the scaling factor used at the m-th time step

blue   = 0.9*[0.0000 0.4470 0.7410];
orange = 0.9*[0.8500 0.3250 0.0980];

% time step
yyaxis left
h1 = semilogy(dt, '-', 'Color', blue, 'LineWidth', 1.5, ...
    'Marker', '+', 'MarkerIndices', [1:124:2605, 2727],'MarkerSize', 10);
ylabel('Time step', 'Interpreter', 'latex', 'FontSize', 20,'FontName', 'Times New Roman')
set(gca,'yTick',10.^(-7:1:-1));
set(gca,'yTickLabel',{'10^{-7}','10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}'});
ylim([10^(-7), 10^(-2)])

% scaling factor
yyaxis right
h2 = plot(mu, '-', 'Color', orange, 'LineWidth', 1.5,...
    'Marker', '*', 'MarkerIndices', [1:124:2605, 2727],'MarkerSize', 10);
ylabel('Scaling factor', 'Interpreter', 'latex','FontSize', 20, 'FontName', 'Times New Roman')
set(gca,'yTick',0:0.5:4.5);
set(gca,'yTickLabel',{'0','0.5','1','1.5','2','2.5','3','3.5','4','4.5'});
ylim([0, 4.5])

% label
set(gca,'FontSize',16); % size of tick
xlabel('Time-step index','FontSize',20,'FontName','Times New Roman') % xlabel
set(gca,'xTick',200:500:2700);
xlim([1,num_step])

% legend
h_legend = legend('time step $\tau$', 'scaling factor $\mu$', 'Interpreter', 'latex');
set(h_legend, 'FontSize', 20, 'FontName', 'Times New Roman', ... 
    'FontWeight', 'normal', 'TextColor', [0, 0, 0]);
set(h_legend,'Position', [0.645, 0.805, 0.0, 0.0]);

%% subfigure: L^inf-norm and H^1-seminorm
figure

% select the last "Num_fit" data for fitting
Num_fit = 500;
TT_fit = TT(end-Num_fit+1:end);
U_max_fit = U_max(end-Num_fit+1:end);
Ugrad_L2_fit = Ugrad_L2(end-Num_fit+1:end);

% fitting L^inf-norm: log(U_max) ~ -k1*log(Tstar-t) + k2
Asy_1 = fitting_Tstar_k1(TT_fit, U_max_fit);
fprintf('RMSE(log-domain) = %.6e\n', Asy_1.RMSE_log); % root mean square error
fprintf('fitting Tstar for L^inf-norm = %.10f\n', Asy_1.Tstar); % the fitting value of T*
fprintf('fitting k1 for L^inf-norm = %.10f\n', Asy_1.k1); % the fitting value of k1
fprintf('fitting k2 for L^inf-norm = %.10f\n', Asy_1.k2); % the fitting value of k2

% fitting H^1-seminorm: log(Ugrad_L2) ~ -k1*log(Tstar-t) + k2
Asy_2 = fitting_Tstar_k1(TT_fit, Ugrad_L2_fit);
fprintf('RMSE(log-domain) = %.6e\n', Asy_2.RMSE_log); % root mean square error
fprintf('fitting Tstar for H^1-seminorm = %.10f\n', Asy_2.Tstar); % the fitting value of T*
fprintf('fitting k1 for H^1-seminorm = %.10f\n', Asy_2.k1); % the fitting value of k1
fprintf('fitting k2 for H^1-seminorm = %.10f\n', Asy_2.k2); % the fitting value of k2

% the average of two fitting blow-up time
Tstar = (Asy_2.Tstar + Asy_1.Tstar)/2;

% plot the numerical results for H^1-seminorm and L^inf-norm
x_axis = (Tstar - TT).^(-1);
loglog(x_axis, U_max, '-', 'LineWidth', 1.5, 'Color',  [1, 0, 1], ...
    'Marker', 's', 'MarkerIndices', [floor(1:2728/20:2728),2728], 'MarkerSize', 8); hold on % L^inf-norm
loglog(x_axis, Ugrad_L2, '-', 'LineWidth', 1.5, 'Color', [0.1034, 0.3967, 0.6000], ...
    'Marker', 'd', 'MarkerIndices', [floor(1:2728/20:2728),2728],'MarkerSize', 8); hold on % H^1-seminorm
xlim([0, 10^4]), ylim([10^0, 10^1.65])

% xlabel
set(gca,'FontSize',16); % size of tick
set(gca,'XTick',10.^(0:0.8:4));
set(gca,'XTickLabel',{'10^{0}','10^{0.8}','10^{1.6}','10^{2.4}','10^{3.2}','10^{4}'});
xlabel('$(T^*-t)^{-1}$','Interpreter','latex','FontSize',20);

% ylabel
set(gca,'YTick',10.^(0:0.3:1.8));
set(gca,'YTickLabel',{'10^{0}','10^{0.3}','10^{0.6}','10^{0.9}', ...
    '10^{1.2}', '10^{1.5}', '10^{1.8}'});
ylabel('$\mathcal{X}$','Interpreter','latex','FontSize',20);

% legend
h_legend = legend('$\mathcal{X} = \Vert \psi(\cdot,t) \Vert_{L^\infty(\mathrm{R}^2)}$', ...
    '$\mathcal{X} = | \psi(\cdot,t) |_{H^1(\mathrm{R}^2)}$', 'Interpreter', 'latex');
set(h_legend, 'FontSize', 20, 'FontName', 'Times New Roman');
set(h_legend,'Position', [0.425, 0.815, 0.0, 0.0]);