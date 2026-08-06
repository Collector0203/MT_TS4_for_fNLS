% This program plots the time evolution of the discrete mass deviation and
% discrete energy deviation for equation parameters (s,p,gamma) =
% (0.6,1,-1).

clear
load('conservation_(0.6,1,-1).mat');

%% Plot the deviation of the discrete mass
figure
drift_mass1 = abs((Mass1-Mass1(1))./Mass1(1));
drift_mass2 = abs((Mass2-Mass2(1))./Mass2(1));

plot((0:Tstep)*dt,drift_mass1,'r-','LineWidth',1.5,'Marker', 'd', ...
    'MarkerIndices', 1:25:(Tstep+1),'MarkerSize', 8), hold on,
plot((0:Tstep)*dt,drift_mass2,'b-','LineWidth',1.5,'Marker', 's', ...
    'MarkerIndices', 1:25:(Tstep+1),'MarkerSize', 8), hold on,

set(gca,'FontSize',16); % size of tick
xlabel('$t_n$','Interpreter','latex','FontSize',20); % xlabel
set(gca,'XTick',[0,5,10,15,20,25]);
set(gca,'XTickLabel',{'0','5','10','15','20','25'});
ylabel('$\mathcal{E}_M(t_n)$','Interpreter','latex','FontSize',20); % ylabel
axis([0,25,10^-30,2.1*10^-12])

h_legend = legend('initial data 1','initial data 2','Interpreter', 'latex');
set(h_legend, 'FontSize', 20, 'FontName', 'Times New Roman', ... 
    'FontWeight', 'normal', 'TextColor', [0, 0, 0]);

%% Plot the deviation of the discrete energy
figure
drift_energy1 = abs((Energy1-Energy1(1))./Energy1(1));
drift_energy2 = abs((Energy2-Energy2(1))./Energy2(1));

plot((0:Tstep)*dt,drift_energy1,'r-','LineWidth',1.5,'Marker', 'd', ...
    'MarkerIndices', 1:25:(Tstep+1),'MarkerSize', 8), hold on,
plot((0:Tstep)*dt,drift_energy2,'b-','LineWidth',1.5,'Marker', 's', ...
    'MarkerIndices', 1:25:(Tstep+1),'MarkerSize', 8), hold on,

set(gca,'FontSize',16); % size of tick
xlabel('$t_n$','Interpreter','latex','FontSize',20); % xlabel
set(gca,'XTick',[0,5,10,15,20,25]);
set(gca,'XTickLabel',{'0','5','10','15','20','25'});
ylabel('$\mathcal{E}_E(t_n)$','Interpreter','latex','FontSize',20); % ylabel
axis([0,25,10^-30,2.5*10^-6])

h_legend = legend('initial data 1','initial data 2','Interpreter', 'latex');
set(h_legend, 'FontSize', 20, 'FontName', 'Times New Roman', ... 
    'FontWeight', 'normal', 'TextColor', [0, 0, 0]);