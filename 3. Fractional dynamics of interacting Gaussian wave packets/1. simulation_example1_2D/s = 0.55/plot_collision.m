% This program uses the data computed by equation2D_collision.m to plot the
% corresponding numerical evolution process.

clear, close all;
load('simulation2D_collision_s=0.55.mat');

%% Compute the values of numerical solution ​​at plotting points
M = 15; % plotting area [-M,M]^2
pt = -M:M/1500:M; [pt2,pt1] = meshgrid(pt); % plotting points [pt1,pt2]

[NN,ppt] = meshgrid(-N:N-1,pt); MTmatrix = MT(NN,mu*ppt); % common matrix, which is used to compute the value of sum(m=-N:N-1)sum(n=-N:N-1) U(m,n)*MT(m,mu*x)*MT(n,mu*y) on plotting points [pt1,pt2] based on U.
unum_coeff = zeros(2*N,2*N,inter+1); % unum_coeff(:,:,q) denotes the MT-coefficients of \psi(x,y,(q-1)*(Tstep/inter)*dt) under {MT(m,mu*x)*MT(n,mu*y)}
for m = 1:inter+1
    unum_coeff(:,:,m) = MTcoeff2(unum_node(:,:,m),N,opt_Sign_SignM,opt_iTheta);
end
unum_nodept = pagemtimes(pagemtimes(MTmatrix,unum_coeff),MTmatrix.'); % unum_nodept(:,:,q) denotes the values of \psi(x,y,(q-1)*(Tstep/inter)*dt) on points [pt1,pt2]

%% Plot the numerical solution
figure('Position', [0, 200, 1800, 400])
plot_par = tight_subplot(1, 4, [0.05 0.05], [0.1 0.1], [0.04 0.015]); 
% parameters of "plot_par": number of rows, number of columns, spacing of subgraphs (子图间距), margins of top and bottom (上下边距) , margins of left and right (左右边距)
for m = 1:4
    axes(plot_par(m))
    mesh(pt1, pt2, abs(unum_nodept(:,:,m))), colorbar
    xlabel('$x_1$','Interpreter','latex','FontSize',14); % xlabel
    ylabel('$x_2$','Interpreter','latex','FontSize',14); % ylabel
    set(gca,'FontSize',14); 
    set(gca,'xTick',-12:6:12); % xTick
    set(gca,'yTick',-12:6:12); % yTick
    title(sprintf('$t = %d$', (m-1)*(Tstep/inter)*dt), 'Interpreter','latex', 'FontSize', 20) % title
    axis([-M, M, -M, M]), daspect([1, 1, 1]), view(2);
end