% This program is used to plot the asymptotic decay of 6 reference
% solutions (d=1,2,3) initialized with the Gaussian function u(x) =
% exp(-|x|^2).

clear

%% Plot the asymptotic decay for 1D
figure; pt = 10.^(1:0.1:2.5); % plot points

% parameters (s,p,gamma,d) = (0.3,1,1,1)
load('ref_solt_gauss_(0.3,1,1,1).mat'); % import the data of reference solution
unum_coeff = unum_coeff_ref(:,1); % the MT-coefficients for T = 1
[NN,pp] = meshgrid(-N_ref:N_ref-1,pt); MTmatrix = MT(NN,mu_ref*pp);
unum_nodept = MTmatrix*unum_coeff; % the values of the numerical solution on pt
thdrate = pt.^(-1-2*s); % theoretical decay rate

loglog(pt,abs(unum_nodept),'bs-', 'MarkerSize', 8, 'LineWidth', 1.5); hold on % plot the numerical solution
loglog(pt, 10^-0.95*thdrate,'b--','LineWidth', 2), hold on % plot the theoretical decay rate

% parameters (s,p,gamma,d) = (0.6,1,-1,1)
load('ref_solt_gauss_(0.6,1,-1,1).mat'); % import the data of reference solution
unum_coeff = unum_coeff_ref(:,1); % the MT-coefficients for T = 1
[NN,pp] = meshgrid(-N_ref:N_ref-1,pt); MTmatrix = MT(NN,mu_ref*pp);
unum_nodept = MTmatrix*unum_coeff; % the values of the numerical solution on pt
thdrate = pt.^(-1-2*s); % theoretical decay rate

loglog(pt,abs(unum_nodept),'gd-', 'MarkerSize', 8, 'LineWidth', 1.5); hold on % plot the numerical solution
loglog(pt, 10^-0.79*thdrate,'g--','LineWidth', 2), hold on % plot the theoretical decay rate

% plotting mark
xlabel('$x$','Interpreter','latex','FontSize',20); % xlabel
ylabel('$|\psi(x,T)|$','Interpreter','latex','FontSize',20); % ylabel
set(gca,'FontSize',16); % size of tick
set(gca,'XTick',10.^(1:0.3:2.5));
set(gca,'XTickLabel',{'10^1','10^{1.3}','10^{1.6}','10^{1.9}','10^{2.2}','10^{2.5}'});
set(gca,'yTick',10.^(-8:0));
set(gca,'yTickLabel',{'10^{-8}','10^{-7}','10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}','10^0'});
set(gca,'TickLabelInterpreter','tex');
axis([10^0.85,10^2.65,10^-6.5,10^-1.8])

h_legend = legend('$s = 0.3$','$\propto x^{-1.6}$','$s = 0.6$', ...
    '$\propto x^{-2.2}$','Interpreter', 'latex');
set(h_legend, 'FontSize', 20, 'FontName', 'Times New Roman', ... 
    'FontWeight', 'normal', 'TextColor', [0, 0, 0]);
set(h_legend,'Position', [0.785, 0.766, 0.0, 0.0]);

%% Plot the asymptotic decay for 2D
figure; pt = 10.^(1:0.1:2.5); % plot points (pt,0)

% parameters (s,p,gamma,d) = (0.3,1,1,2)
load('ref_solt_gauss_(0.3,1,1,2).mat'); % import the data of reference solution
unum_coeff = unum_coeff_ref(:,:,2); % the MT-coefficients for T = 1
[NN,pp] = meshgrid(-N_ref:N_ref-1,pt);
MTmatrix1 = MT(NN,mu_ref*pp); MTmatrix2 = MT(-N_ref:N_ref-1,0).';
unum_nodept = MTmatrix1*unum_coeff*MTmatrix2; % the values of the numerical solution on (pt,0)
thdrate = pt.^(-2-2*s); % theoretical decay rate

loglog(pt,abs(unum_nodept),'bs-', 'MarkerSize', 8, 'LineWidth', 1.5); hold on % plot the numerical solution
loglog(pt, 10^-0.49*thdrate,'b--','LineWidth', 2), hold on % plot the theoretical decay rate

% parameters (s,p,gamma,d) = (0.6,1,-1,2)
load('ref_solt_gauss_(0.6,1,-1,2).mat'); % import the data of reference solution
unum_coeff = unum_coeff_ref(:,:,2); % the MT-coefficients for T = 1
[NN,pp] = meshgrid(-N_ref:N_ref-1,pt);
MTmatrix1 = MT(NN,mu_ref*pp); MTmatrix2 = MT(-N_ref:N_ref-1,0).';
unum_nodept = MTmatrix1*unum_coeff*MTmatrix2; % the values of the numerical solution on (pt,0)
thdrate = pt.^(-2-2*s); % theoretical decay rate

loglog(pt,abs(unum_nodept),'gd-', 'MarkerSize', 8, 'LineWidth', 1.5); hold on % plot the numerical solution
loglog(pt, 10^-0.95*thdrate,'g--','LineWidth', 2), hold on % plot the theoretical decay rate

% plotting mark
xlabel('$x_1$','Interpreter','latex','FontSize',20); % xlabel
ylabel('$|\psi(x_1,0,T)|$','Interpreter','latex','FontSize',20); % ylabel
set(gca,'FontSize',16); % size of tick
set(gca,'XTick',10.^(1:0.3:2.5));
set(gca,'XTickLabel',{'10^1','10^{1.3}','10^{1.6}','10^{1.9}','10^{2.2}','10^{2.5}'});
set(gca,'yTick',10.^(-10:0));
set(gca,'yTickLabel',{'10^{-10}','10^{-9}','10^{-8}','10^{-7}','10^{-6}','10^{-5}','10^{-4}', ...
    '10^{-3}','10^{-2}','10^{-1}','10^0'});
set(gca,'TickLabelInterpreter','tex');
axis([10^0.85,10^2.65,10^-9.3,10^-2.7])

h_legend = legend('$s = 0.3$','$\propto x_1^{-2.6}$','$s = 0.6$', ...
    '$\propto x_1^{-3.2}$','Interpreter', 'latex');
set(h_legend, 'FontSize', 20, 'FontName', 'Times New Roman', ... 
    'FontWeight', 'normal', 'TextColor', [0, 0, 0]);
set(h_legend,'Position', [0.785, 0.766, 0.0, 0.0]);

%% Plot the asymptotic decay for 3D
figure; pt = 10.^(1:0.1:2.5); % plot points (pt,0,0)

% parameters (s,p,gamma,d) = (0.3,1,1,3)
load('ref_solt_gauss_(0.3,1,1,3).mat'); % import the data of reference solution
unum_coeff = unum_coeff_ref(:,:,:,2); % the MT-coefficients for T = 1
[NN,pp] = meshgrid(-N_ref:N_ref-1,pt); 
MTmatrix1 = MT(NN,mu_ref*pp); MTmatrix2 = MT(-N_ref:N_ref-1,0).';
unum_nodept_temp = pagemtimes(unum_coeff,MTmatrix2); 
unum_nodept_temp = squeeze(unum_nodept_temp);
unum_nodept = MTmatrix1*unum_nodept_temp*MTmatrix2; % the values of the numerical solution on (pt,0,0)
thdrate = pt.^(-3-2*s); % theoretical decay rate

loglog(pt,abs(unum_nodept),'bs-', 'MarkerSize', 8, 'LineWidth', 1.5); hold on % plot the numerical solution
loglog(pt, 10^-0.35*thdrate,'b--','LineWidth', 2), hold on % plot the theoretical decay rate

% parameters (s,p,gamma,d) = (0.6,1,-1,3)
load('ref_solt_gauss_(0.6,1,-1,3).mat'); % import the data of reference solution
unum_coeff = unum_coeff_ref(:,:,:,2); % the MT-coefficients for T = 1
[NN,pp] = meshgrid(-N_ref:N_ref-1,pt);
MTmatrix1 = MT(NN,mu_ref*pp); MTmatrix2 = MT(-N_ref:N_ref-1,0).';
unum_nodept_temp = pagemtimes(unum_coeff,MTmatrix2);
unum_nodept_temp = squeeze(unum_nodept_temp);
unum_nodept = MTmatrix1*unum_nodept_temp*MTmatrix2; % the values of the numerical solution on (pt,0,0)
thdrate = pt.^(-3-2*s); % theoretical decay rate

loglog(pt,abs(unum_nodept),'gd-', 'MarkerSize', 8, 'LineWidth', 1.5); hold on % plot the numerical solution
loglog(pt, 10^-1.0*thdrate,'g--','LineWidth', 2), hold on % plot the theoretical decay rate

% plotting mark
xlabel('$x_1$','Interpreter','latex','FontSize',20); % xlabel
ylabel('$|\psi(x_1,0,0,T)|$','Interpreter','latex','FontSize',20); % ylabel
set(gca,'FontSize',16); % size of tick
set(gca,'XTick',10.^(1:0.3:2.5));
set(gca,'XTickLabel',{'10^1','10^{1.3}','10^{1.6}','10^{1.9}','10^{2.2}','10^{2.5}'});
set(gca,'yTick',10.^(-12:2:-4));
set(gca,'yTickLabel',{'10^{-12}','10^{-10}','10^{-8}','10^{-6}','10^{-4}'});
set(gca,'TickLabelInterpreter','tex');
axis([10^0.85,10^2.65,10^-12,10^-3.5])

h_legend = legend('$s = 0.3$','$\propto x_1^{-3.6}$','$s = 0.6$', ...
    '$\propto x_1^{-4.2}$','Interpreter', 'latex');
set(h_legend, 'FontSize', 20, 'FontName', 'Times New Roman', ... 
    'FontWeight', 'normal', 'TextColor', [0, 0, 0]);
set(h_legend,'Position', [0.785, 0.766, 0.0, 0.0]);