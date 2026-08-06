% This program uses the data computed by equation2D_blow_gauss.m to plot
% the corresponding numerical solutions.

clear
load('BioMT_TSM4_s=0.7_p=1.5_N=256_Bound=20.5.mat');

%% Plot the numerical solutions: 2D figure
M = 2; % plotting area [-M,M]^2
pt = -M:M/2500:M; [pt2,pt1] = meshgrid(pt); % plotting points [pt1,pt2]
[NN,ppt] = meshgrid(-N:N-1,pt); % usd in computing the numerical values of \psi(x,t)

figure('Position', [0, 100, 1800, 400])
plot_par = tight_subplot(1, 4, [0.05 0.05], [0.1 0.1], [0.04 0.015]); 
xticks = -2:1:2;
yticks = xticks;
zticks = 0:5:20;

for k = 1:4
    axes(plot_par(k))
    %% read data
    flag = coeff_store{k,1};
    t = TT(flag); mu = Mu(flag); % read the scaling factor and time corresponding to "flag"
    coeff = MTcoeff2_origin(coeff_store{k,2},N); % read the MT-coefficients corresponding to "flag"

    %% plot numerical solutions
    unum_nodept = MT(NN,mu*ppt)*coeff*MT(NN,mu*ppt).'; % numerical values of \psi(x,t) on [pt1,pt2]
    mesh(pt1, pt2, abs(unum_nodept)), grid off,
    colormap(gca, 'turbo');
    xlim([pt1(1), pt1(end)]);
    ylim([pt2(1), pt2(end)]);
    zlim([0, 20]), daspect([1 1 4.5])

    %% view
    view(50,15)

    %% label and axis
    hx = xlabel('$x_1$','Interpreter','latex','FontSize',14);
    hy = ylabel('$x_2$','Interpreter','latex','FontSize',14);
    hz = zlabel('$|\psi|$','Interpreter','latex','FontSize',14);
    set(gca,'XTick',xticks); set(gca,'XTickLabel',num2str(xticks'),'FontSize',14);
    set(gca,'YTick',yticks); set(gca,'YTickLabel',num2str(yticks'),'FontSize',14);
    set(gca,'ZTick',zticks); set(gca,'ZTickLabel',num2str(zticks'),'FontSize',14);

    drawnow;
    hx.Units = 'data';
    pos = hx.Position;
    pos(1) = pos(1)-0.7;
    pos(2) = pos(2)+0.5;
    hx.Position = pos; % modify the position of xlabel

    drawnow;
    hy.Units = 'data';
    pos = hy.Position;
    pos(1) = pos(1)-0.4;
    pos(2) = pos(2)+0.8;
    hy.Position = pos; % modify the position of ylabel

    %% title
    title(sprintf(['$\\Vert \\psi(\\cdot,t = %.',num2str(k),'f) \\Vert_{L^\\infty(\\mathrm{R}^2)} = %.2f$'], t, U_max(flag)),'Interpreter','latex','FontSize',21);
end

%% Plot the numerical solutions: Far-field behavior
figure('Position', [0, 100, 1800, 360])
plot_par = tight_subplot(1, 4, [0.05 0.05], [0.18 0.1], [0.04 0.015]);

for k = 1:4
    axes(plot_par(k))
    %% compute and plot the numerical solution in the x1-direction
    pt = 10.^(1:0.1:2.5); % plot points (pt,0)
    flag = coeff_store{k,1}; coeff = MTcoeff2_origin(coeff_store{k,2},N); mu = Mu(flag); % read data
    [NN,pp] = meshgrid(-N:N-1,pt); MTmatrix1 = MT(NN,mu*pp); MTmatrix2 = MT(-N:N-1,0).';
    unum_nodept = MTmatrix1*coeff*MTmatrix2; % compute the numerical solution in the x1-direction
    loglog(pt,abs(unum_nodept),'-s','MarkerSize', 8, 'LineWidth', 1.5); hold on % plot the numerical solution in the x1-direction

    %% plot the theoretical decay rate
    thdrate = pt.^(-2-2*s); % theoretical decay rate
    loglog(pt, 10^-0.49*thdrate,'b--','LineWidth', 2), hold on

    %% label
    set(gca,'FontSize',14.5); % size of tick
    xlabel('$x_1$','Interpreter','latex','FontSize',16); % xlabel
    ylabel('$|\psi(x_1,0,t)|$','Interpreter','latex','FontSize',16); % ylabel
    set(gca,'XTick',10.^(1:0.3:2.5));
    set(gca,'XTickLabel',{'10^1','10^{1.3}','10^{1.6}','10^{1.9}','10^{2.2}','10^{2.5}'});
    set(gca,'yTick',10.^(-9:0));
    set(gca,'yTickLabel',{'10^{-9}','10^{-8}','10^{-7}','10^{-6}','10^{-5}','10^{-4}', ...
        '10^{-3}','10^{-2}','10^{-1}','10^0'});
    set(gca,'TickLabelInterpreter','tex');
    axis([10^1,10^2.5,10^-9,10^-3])

    %% legend
    t = TT(flag);
    fac = 10^k; t_show = floor(t*fac)/fac; % truncate t, not round
    h_legend = legend(sprintf(['$t =  %.',num2str(k),'f$'], t_show),'$\propto x_1^{-3.4}$','Interpreter', 'latex');
    set(h_legend, 'FontSize', 17, 'FontName', 'Times New Roman', 'FontWeight', 'normal', 'TextColor', [0, 0, 0]);
end